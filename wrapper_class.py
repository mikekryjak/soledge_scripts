from files.load_soledge_mesh_file				import load_soledge_mesh_file
from files.load_plasma_files						import load_plasma_files
from files.load_ions_list								import load_ions_list
from files.load_refpar_file				import load_refpar_file
from files.load_eirene_triangles					import load_eirene_triangles
from mesh.get_rz_core_sep						import get_rz_core_sep
from mesh.get_rho_in_out_core_sep		import get_rho_in_out_core_sep
from mesh.find_zones_intersections			import find_zones_intersections
from mesh.compute_mesh_intersections	import compute_mesh_intersections
from routines.intersect_contour				import intersect_2contours
from routines.utils_walls							import get_in_out_walls, plot2d_walls, get_dmax_points_walls
from routines.h5_routines							import h5_read
from eirene.get_wall_triangle					import get_wall_triangle
from math													import sqrt, exp
import numpy											as np
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import contextlib
import h5py
import scipy


class SOLEDGEcase():
    def __init__(self, 
                 path,
                 verbose = False):
        
        with HiddenPrints() if verbose is False else contextlib.nullcontext():  # Suppress all prints
            
            self.path = path 
            self.verbose = verbose
            
            #	Read mesh
            self.Config = load_soledge_mesh_file(os.path.join(path,"mesh.h5"))
            
            # Load species names
            self.ions = load_ions_list(path)   # Includes electrons and ions (same nomenclature as original)

            # Load species parameters from both the plasma and EIRENE grid.
            # The first index of Plasmas is ions/neutrals and second is electrons
            self.Plasmas = load_plasma_files(path, nZones=len(self.Config.Zones), Evolution=0, iPlasmas=range(len(self.ions)))

            # Load parameter names
            self.iparams = self.Plasmas[0][0].VNames
            self.eparams = self.Plasmas[1][0].VNames
            self.params = self.iparams + self.eparams
            
            # Load reference parameters
            self.RefPar = load_refpar_file(os.path.join(path, "Results"))
            
            # Load triangle grid
            if_tri	 = h5py.File(os.path.join(path,"triangles.h5"), "r")
            self.TriKnots = h5_read(if_tri,"triangles/tri_knots", messages = 0)   # Nodes of all EIRENE triangles
            self.TriKnots = self.TriKnots - 1 										#Matlab/Fortan to python indexes
            self.R		 = h5_read(if_tri,"knots/R", messages = 0)*0.01
            self.Z		 = h5_read(if_tri,"knots/Z", messages = 0)*0.01
            if_tri.close()
            self.TripTriang = mpl.tri.Triangulation(self.R, self.Z, triangles=self.TriKnots)
            self.Eirene = load_eirene_triangles(os.path.join(path, "triangles.h5"))
            
            ### Plasma grid extents and wall
            self.Rlim = [self.R.min(), self.R.max()]
            self.Zlim = [self.Z.min(), self.Z.max()]
            if_mesh = h5py.File(os.path.join(path, "mesh.h5"), "r")

            try:
                self.Rwall	= h5_read(if_mesh, "walls/wall1/R")
                self.Zwall	= h5_read(if_mesh, "walls/wall1/Z")
            except:
                self.Rwall	= h5_read(if_mesh, "wall/R")
                self.Zwall	= h5_read(if_mesh, "wall/Z")

            if_mesh.close()
            
            # EIRENE wall extents
            ZeroTriangle, ZeroSide, Ri, Zi = get_wall_triangle(self.Eirene, rz0_line=[0,0], theta_line=0, no_plot=1, no_print=1)
            
    
    def plot_2d(
        self,
        fig,
        ax,
        param = None, 
        data = None,
        norm = None,
        vmin = None, 
        vmax = None, 
        logscale = True, 
        separatrix = True, 
        sep_width = 2,
        cmap = "Spectral_r",
        cbar = False,
        cbar_label = "",
        linewidth = 0,
        linecolor = "k",
        grid_only = False,
        verbose = False):
        """
        Plots a 2D tripcolor plot on tri grid
        Tri grid is useful because it contains the actual wall unlike the hex grid
        
        Inputs
        -----
        - param = str: Parameter name in SOLEDGE2D nomenclature
        - ax = matplotlib axis object: Axis to plot on
        - norm = matplotlib norm: If None, one is created
        - vmin, vmax: Minimum and maximum for colormap if no norm is passed
        - logscale = bool: For creating norm if none is passed. True creates symlog.
        - separatrix = bool: plot separatrix? This is currently done from eqb instead of grid.
        - sep_width = int: line width of the separatrix plot
        - cmap = str: Matplotlib colormap choice
        """
        
        # if param == None and data == []:
        #     raise Exception("Must provide param or data")
        # if param != None and data.isinstance(np.array):
        #     raise Exception("Must provide param or data")
        
        if param != None:
            iPlasma, iPar = self.get_param_indices(param, triangles = True)
            data = self.Plasmas[iPlasma][0].Triangles.Values[iPar]
        
        
        # Calculate range if not provided
        vmin = min(data) if vmin is None else vmin
        vmax = max(data) if vmax is None else vmax
        
        print(vmin,vmax)
        
        if norm is None:
            norm = create_norm(logscale, norm, vmin, vmax)
        if grid_only is True:
            cmap = mpl.colors.ListedColormap(["white"])
        
        tp = ax.tripcolor(self.TripTriang, data, norm = norm, cmap = cmap,  linewidth=linewidth, edgecolors = linecolor)
        ax.set_aspect("equal")
        if separatrix is True:
            lw = sep_width
            c = "w"
            lhs, rhs = self._get_rz_sep()
            ax.plot(lhs["R"], lhs["Z"], lw = lw, c = c)
            ax.plot(rhs["R"], rhs["Z"], lw = lw, c = c)
            
        if cbar is True:
            fig.colorbar(tp, norm=norm, label = cbar_label)

            
        # self.params = self.Plasmas[iPlasma][0].Triangles.VNames + self.Plasmas[iPlasma][1].Triangles.VNames
          
    def _get_2d_on_tri(self,param):
        """
        Return plottable 2D data on triangle mesh
        """
        
        iPlasma, iPar = self.get_param_indices(param, triangles = True)
        data = self.Plasmas[iPlasma][0].Triangles.Values[iPar]
        
        return data
            
    def get_1d_radial_data(
        self, 
        params, 
        rz0_line = [], 
        theta_line = 0, 
        verbose = False,
        debug_line = False):
        
        """ 
        Prepares radial 1D plot
        
        Inputs
        ------
        - param = string: Parameter name in SOLEDGE2D nomenclature
        - rz0_line = [R,Z]: RZ coordinates of the origin of the line. if 0, will calculate eqb centre coords
        - theta_line = 0: Angle of the line, where 0 degrees is the OMP and 180 is the IMP
        - debug_line = bool: produce a plot showing the location of the line compared to the separatrix
        
        Returns
        ------
        Dataframe with poloidal/parallel length, RZ coordinates and parameter data.
        The data column is named after the parameter of choice
        """
        
        with HiddenPrints() if verbose is False else contextlib.nullcontext():  # Suppress all prints
        
            # Find mesh cells along line
            if(len(rz0_line) == 0):
                
                Rcore, Zcore, CoreMegazone = get_rz_core_sep(self.Config, core_and_sep = False)
                rz0_line = [0.5*(Rcore.min() + Rcore.max()), Zcore[np.argmax(Rcore)]]
                if verbose is True: print(f"rz0_line not set. Calculated as {rz0_line}")
                
            rMax		= 6*get_dmax_points_walls(self.Config, rz0_line[0], rz0_line[1], plasma_wall=True, eirene_wall=False, extra_wall=False)
            theta_line	= theta_line*np.pi/180.
            RZLine		= np.array([[rz0_line[0],						  rz0_line[1]], \
                                    [rz0_line[0]+rMax*np.cos(theta_line), rz0_line[1]+rMax*np.sin(theta_line)]])


            # Calculate the OMP line
            Cut = find_zones_intersections(self.Config, RZLine)
            Lengths, IntRZ, IntCEll = compute_mesh_intersections(self.Config, Cut, also_pos=True, use_mag_zones=False)
            in_wall, out_wall = get_in_out_walls(self.Config, IntRZ[:,0], IntRZ[:,1])
            Lengths = Lengths[in_wall,:]
            IntRZ	= IntRZ[in_wall,:]    # Selected cell RZ coordinates
            IntCEll	= IntCEll[in_wall,:]   # Selected cell indices in [pol, rad, corner]
            
            #    Calculate separatrix distances for the X axis of the plot
            dist = Lengths[:,0] - Lengths[0,0]
            
            #    dist = distance from separatrix
            #    In_Sep, Out_Sep = indices of the above for inside and outside separatrix
            #    RZcore, RZsep = RZ coordinates for inside and outside separatrix
            #    xName can be Psi_pol or Rho_pol... need to look into it
            Rho, In_Sep, Out_Sep, RZcore, RZsep = get_rho_in_out_core_sep(self.Config, IntRZ[:,0], IntRZ[:,1], rho_type = "Psi_pol")

            if(len(In_Sep) > 0):
                if(Out_Sep[-1] < In_Sep[0]):
                    Out_Sep = np.append(Out_Sep, In_Sep[0])
                else:
                    Out_Sep = np.append(In_Sep[-1], Out_Sep)

            xName = "DIST"
            Ri, Zi, is1, is2  = intersect_2contours(RZsep[:,0], RZsep[:,1], IntRZ[:,0], IntRZ[:,1])
            if(len(Ri)==0):
                if verbose is True: print(f"Contour intersection failed!")
                plt.plot(RZsep[:,0], RZsep[:,1],'k-')
                plt.plot(IntRZ[:,0], IntRZ[:,1],'r-')
                plt.show()

            dsep = sqrt((IntRZ[0,0] - Ri[0])**2 + (IntRZ[0,1] - Zi[0])**2)
            dist -= dsep
            xSep = 0.
            
            df = pd.DataFrame()
            
            df["dist"] = dist
            df["R"] = IntRZ[:,0]
            df["Z"] = IntRZ[:,1]
            
            for param in params:
                iPlasma, iPar = self.get_param_indices(param)
                data = get_plasma_parameter_on_mesh(self.Plasmas[iPlasma], iPar, IntCEll)
                df[param] = data
            
            if debug_line is True:
                fig, ax = plt.subplots()
                ax.plot(RZsep[:,0], RZsep[:,1])
                ax.plot(IntRZ[:,0], IntRZ[:,1])
                ax.set_aspect("equal")
                ax.grid()
                ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
        
        return df
    
    
    def get_1d_poloidal_data(
        self, 
        params, 
        rz0_line = [], 
        theta_line = 0, 
        parallel_length = False,
        d_from_sep = 0.0001,
        verbose = False):
        
        """
        Prepares data for a 1D poloidal or parallel plot. Based on original SOLEDGE2D routines
        
        Inputs
        ------
        - param = string: Parameter name in SOLEDGE2D nomenclature
        - rz0_line = [R,Z]: Origin of radial line to use as basis. If default, code will find domain centre
        - theta_line = 0 : Bearing of basis radial line in degrees from 0 at the OMP (i.e. IMP is 180)
        - parallel_length = bool: Use parallel instead of poloidal length for X axis
        - d_from_sep = 0.01: Ideal distance of poloidal cut from the separatrix
        - verbose = bool: Print more info
        
        Returns
        ------
        Dataframe with poloidal/parallel length, RZ coordinates and parameter data.
        The data column is named after the parameter of choice
        """
        
        with HiddenPrints() if verbose is False else contextlib.nullcontext():  # Suppress all prints
        
            Zones = self.Config.Zones

            # Read metric coefficients if parallel length requested

            l_pol = 1 if parallel_length is False else 0 

            if(l_pol == 0):
                if_metric = h5py.File(os.path.join(path, "Results/metric"), "r")
                Gmet = []
                for k in range(len(self.Config.Zones)):
                    zone = "zone{:d}".format(k+1)
                    Gmet.append(h5_read(if_metric, zone+ '/G', order = 'F'))									#[Nx+2,Nz+2]

                if_metric.close()
                
            # Generate radial cut to form basis of poloidal cut

            if(len(rz0_line) == 0):
                Rcore, Zcore, CoreMegazone = get_rz_core_sep(self.Config, core_and_sep = False)
                rz0_line = [0.5*(Rcore.min() + Rcore.max()), 0.]

            rMax		= 6*get_dmax_points_walls(self.Config, rz0_line[0], rz0_line[1], plasma_wall=True, eirene_wall=False, extra_wall=False)
            theta_line	= theta_line*np.pi/180.
            RZLine		= np.array([[rz0_line[0],						  rz0_line[1]], \
                                    [rz0_line[0]+rMax*np.cos(theta_line), rz0_line[1]+rMax*np.sin(theta_line)]])


            Cut = find_zones_intersections(self.Config, RZLine)
            Lengths, IntRZ, IntCEll = compute_mesh_intersections(self.Config, Cut, also_pos=True, use_mag_zones=False)
            in_wall, out_wall = get_in_out_walls(self.Config, IntRZ[:,0], IntRZ[:,1])
            Lengths = Lengths[in_wall,:]
            IntRZ	= IntRZ[in_wall,:]
            IntCEll	= IntCEll[in_wall,:]

            dist = Lengths[:,0] - Lengths[0,0]
            Rho, In_Sep, Out_Sep, RZcore, RZsep = get_rho_in_out_core_sep(self.Config, IntRZ[:,0], IntRZ[:,1])

            if(len(In_Sep) > 0):
                if(Out_Sep[-1] < In_Sep[0]):
                    Out_Sep = np.append(Out_Sep, In_Sep[0])
                else:
                    Out_Sep = np.append(In_Sep[-1], Out_Sep)

            Ri, Zi, is1, is2  = intersect_2contours(RZsep[:,0], RZsep[:,1], IntRZ[:,0], IntRZ[:,1])
            if(len(Ri)==0):
                plt.plot(RZsep[:,0], RZsep[:,1],'k-')
                plt.plot(IntRZ[:,0], IntRZ[:,1],'r-')
                plt.show()

            dsep = np.sqrt((IntRZ[0,0] - Ri[0])**2 + (IntRZ[0,1] - Zi[0])**2)
            dist -= dsep   # Distance from separatrix


            # Select cell

            iCell = np.argmin(np.abs(dist-d_from_sep))
            iZone 	= IntCEll[iCell,0]
            ix		= IntCEll[iCell,1]
            iTheta	= IntCEll[iCell,2]

            # Find zones along poloidal coordinate

            if(Zones[iZone].Chi[ix,-1] == 1):
                iThEast = np.array([np.min(np.where(Zones[iZone].Chi[ix,iTheta:] == 1)[0])+iTheta])
                East = -1
            else:
                iThEast = np.array([Zones[iZone].Chi.shape[1]])
                East = Zones[iZone].Neighbour.east

            if(Zones[iZone].Chi[ix,0] == 1):
                iThWest = np.array([np.max(np.where(Zones[iZone].Chi[ix,:iTheta] == 1)[0])])
                West = -1
            else:
                iThWest = np.array([0])
                West = Zones[iZone].Neighbour.west

            iThetaOff  = iTheta - iThWest[0]
            nThetaPts  = iThEast[0] - iThWest[0]
            iZones	   = np.array([iZone])

            #   Look East

            while (East > -1):
                iZones = np.append(iZones, East)
                iThWest = np.append(iThWest,0)
                if(Zones[East].Chi[ix,-1] == 1):
                    iThEast = np.append(iThEast, np.min(np.where(Zones[East].Chi[ix,:] == 1)[0]))
                    East = -1
                else:
                    iThEast = np.append(iThEast, Zones[East].Chi.shape[1])
                    East 	 = Zones[East].Neighbour.east
                nThetaPts += iThEast[-1]

            #   Look West

            while (West > -1):
                iZones = np.append(West, iZones)
                iThEast = np.append(Zones[West].Chi.shape[1], iThEast)
                if(Zones[West].Chi[ix,0] == 1):
                    iThWest = np.append(np.max(np.where(Zones[West].Chi[ix,:] == 1)[0])+1, iThWest)
                    West = -1
                else:
                    iThWest = np.append(0, iThWest)
                    West = Zones[West].Neighbour.west
                iThetaOff += iThEast[0] - iThWest[0]
                nThetaPts += iThEast[0] - iThWest[0]

            Rpol = np.empty((nThetaPts), dtype = 'f8')
            Zpol = np.empty((nThetaPts), dtype = 'f8')
            jOff = 0
            for k in range(len(iZones)):
                Rpol[jOff: jOff + iThEast[k] - iThWest[k]] = Zones[iZones[k]].gridRc[ix, iThWest[k]:iThEast[k]]
                Zpol[jOff: jOff + iThEast[k] - iThWest[k]] = Zones[iZones[k]].gridZc[ix, iThWest[k]:iThEast[k]]
                jOff += iThEast[k] - iThWest[k]

            if(l_pol == 0):
                dl = np.empty((nThetaPts), dtype = 'f8')
                jOff = 0
                for k in range(len(iZones)):
                    dtheta = (Zones[iZones[k]].zb[ix, 1:] - Zones[iZones[k]].zb[ix, :-1])*2.*np.pi
                    dlZone	   = 2.*dtheta/Gmet[iZones[k]][ix+1, 1:-1]
                    dl[jOff: jOff + iThEast[k] - iThWest[k]] = dlZone[iThWest[k]:iThEast[k]]
                    jOff += iThEast[k] - iThWest[k]
                dlZone = 0
                Lpara = np.cumsum(dl)
                dl	  = 0
            else:
                Lpara = np.sqrt((Rpol[1:]-Rpol[:-1])**2 + (Zpol[1:]-Zpol[:-1])**2)
                Lpara = np.append(0., np.cumsum(np.sqrt((Rpol[1:]-Rpol[:-1])**2 + (Zpol[1:]-Zpol[:-1])**2)))

            Lpara = Lpara - Lpara[iThetaOff]   # Parallel or poloidal distance starting at midplane, depending on setting l_pol
            
            
            # Grab the data, pack in dataframe and return
            df = pd.DataFrame()
            for param in params:
                iPlasma, iPar = self.get_param_indices(param)
                data = get_plasma_parameter_on_pol(self.Plasmas[iPlasma], iPar, ix, iZones, iThWest, iThEast, nThetaPts)
                df[param] = data
            
           
            
            
            df["dist"] = Lpara
            df["R"] = Rpol
            df["Z"] = Zpol
            df[param] = data
        
        return df
    
    def get_param_indices(self, param, verbose = False, triangles = False):
        """
        Parse a parameter name to get the species and parameter indices
        """
        
        with HiddenPrints() if verbose is False else contextlib.nullcontext():  # Suppress all prints
            
            # Parse parameter name
            #    Get species index
            if param in self.iparams:   # heavy species
                iPlasma = 0
            elif param in self.eparams:   # electrons
                iPlasma = 1
            else:
                raise Exception(f"Parameter {param} not found")
            
            # Get parameter index
            #    Not all data is interpolated onto triangles so it's a different list of params
            if triangles is True:
                iPar = self.Plasmas[iPlasma][0].Triangles.VNames.index(param) 
            else:
                iPar = self.Plasmas[iPlasma][0].VNames.index(param) 
            
        return iPlasma, iPar
    
    def get_wall_fluxes(self, verbose = True, split = True, grid = "tight"):
        """
        Return wall particle and heat fluxes for the entire wall. If split = true,
        will split the fluxes between  the walls, PFR and targets. This relies 
        on hardcoded grid indices and is only implemented for the "tight" grid right now.
        Verbose will give you a print of the results as well as a debug plot showing the regions.
        """
        
        ### Get wall coordinates and lengths
        ZeroTriangle, ZeroSide, Ri, Zi = get_wall_triangle(self.Eirene, rz0_line=[0,0], theta_line=0, no_plot=1, no_print=1)

        R = self.R
        Z = self.Z
        TriKnots = self.TriKnots
        
        try:
            AFluxF	 = np.loadtxt(os.path.join(self.path, "soledge2D.atoms_fluxes_wall_1"), dtype='f8')
        #            		[0:iTri, 1:iSide, 2:iProp, 3:dlWall, 
        #					 4:Flux_Atom_Incident, 			5:Flux_Mol_Incident,
        #					 6:Flux_Atom_Emitted_from_Ion, 	7:Flux_Mol_Emitted_from_Ion,
        #					 8:Flux_Atom_Emitted_from_Atom, 9:Flux_Mol_Emitted_from_Atom]
        except:
            print("\tError: Not found "+ self.path+"soledge2D.atoms_fluxes_wall_1")
            exit()

        iTri	= AFluxF[:,0]-1 										#Matlab/Fortan to python indexes
        iSide	= AFluxF[:,1]-1 										#Matlab/Fortan to python indexes
        WallL	= AFluxF[:,3]


        Vertex	= np.array([[0,1],[1,2],[2,0]], dtype='i4')
        iTri	= iTri.astype(int)
        iSide	= iSide.astype(int)

        iZero = np.where(iTri == ZeroTriangle); iZero = iZero[0]
        if(len(iZero) != 1):
            print("\tError: not found wall triangle intersection with line")
            exit()

        iZero  = iZero[0]
        zOrder = np.arange(len(iTri))
        zOrder = np.append(zOrder[iZero:], zOrder[:iZero]) 
        iTri   = iTri[zOrder]

        iSide  = iSide[zOrder]
        WallL[iZero:]  = WallL[iZero:] - WallL[iZero]
        WallL[:iZero]  = WallL[:iZero] - WallL[0] + WallL[-1]
        WallL		   = WallL[zOrder]

        WallR	= 0.5*(R[TriKnots[iTri,Vertex[iSide,0]]]+R[TriKnots[iTri,Vertex[iSide,1]]])
        WallZ	= 0.5*(Z[TriKnots[iTri,Vertex[iSide,0]]]+Z[TriKnots[iTri,Vertex[iSide,1]]])
        WalldL  = np.sqrt((R[TriKnots[iTri,Vertex[iSide,0]]]-R[TriKnots[iTri,Vertex[iSide,1]]])**2 + \
                            (Z[TriKnots[iTri,Vertex[iSide,0]]]-Z[TriKnots[iTri,Vertex[iSide,1]]])**2)
        
        df = pd.DataFrame(columns = ["iTri", "iSide", "iProp", "dlSurf", "dlWall", "Area"])
        df["R"] = WallR
        df["Z"] = WallZ

        ## --------------- Read ion fluxes
        try:
            FluxiF	 = np.loadtxt(os.path.join(self.path,"soledge2D.ion_fluxes_wall_1"), dtype='f8')
        except:
            print("\tError: Not found "+ self.path+"soledge2D.ion_fluxes_wall_1")
            exit()
        Fluxi =  FluxiF[zOrder,4]*1e-22

        for i, name in enumerate(["iTri", "iSide", "iProp", "dlSurf", 
                        "F_Ion_Incident", "F_Atom_Emitted_from_Ion"]):
            
            df[name] = FluxiF[zOrder, i]
            
        ## --------------- Read neutral fluxes (again for some reason)
        # FluxiN is a list of arrays, each a different neutral species
        FluxiN=[]
        try:
            FluxiN.append(np.loadtxt(os.path.join(self.path, "soledge2D.atoms_fluxes_wall_1"), dtype='f8'))
            i=2
            while(os.path.isfile('./soledge2D.atoms_fluxes_wall_'+str(i))):
                FluxiN.append(np.loadtxt(self.path+"soledge2D.atoms_fluxes_wall_"+str(i), dtype='f8'))
                i=i+1
        except:
            print("\tError: Not found "+ self.path+"soledge2D.atoms_fluxes_wall_1")
            exit()
            
        for i, name in enumerate(["iTri", "iSide", "iProp", "dlWall", 
                    "F_Atom_Incident", 			"F_Mol_Incident",
                    "F_Atom_Emitted_from_Ion", 	"F_Mol_Emitted_from_Ion",
                    "F_Atom_Emitted_from_Atom", "F_Mol_Emitted_from_Atom"]):
            
            df[name] = FluxiN[0][zOrder, i]

        ## --------------- Read neutral energy fluxes
        try:

            EFluxF	 = np.loadtxt(os.path.join(self.path, "soledge2D.energy_fluxes_details_1"), dtype='f8', comments="%")
            
            for i, name in enumerate(
                ["dlWall", "E_incident_Electron", 	"E_incident_Ions",
                "E_Total",				"E_incident_Atoms",
                "E_Radiation_Atoms",		"E_Rad_Recombination",
                "E_Recombination_in_Wall", "E_incident_Molecules"]):
                
                if name not in ["dlWall"]:
                    df[name] = EFluxF[zOrder,i]

            EFluxF	 = np.loadtxt(os.path.join(self.path, "soledge2D.energy_fluxes_1"), dtype='f8', comments="%")
            
            for i, name in enumerate(
                [
                "iTri", "iSide",                   "iProp", "dlWall", 
                "E_incident_Ions_EIRENE", "E_incident_Electron",
                "E_incident_Ions",		   "E_Total",				
                "E_incident_Atoms",	   "E_Radiation",		
                "Area"
                ]
            ):
                if not any([x in name for x in ["iTri", "iSide", "iProp", "dlWall"]]):
                    df[name] = EFluxF[zOrder, i]
                
        except:
            raise Exception("\tError: Not found "+ self.path+"soledge2D.energy_fluxes_details_1")
            
        ## Not sure what's happening here
        # fluxN=[]	#Neutral flux to the wall	[[fD][fD2][fI1]...[fIn]]
        # fluxN.append(FluxiN[0][zOrder,4])
        # fluxN.append(FluxiN[0][zOrder,5])
        # for i in range(len(FluxiN)-1):
        #     fluxN.append(FluxiN[i+1][zOrder,4])
        # for i in range(len(fluxN)):
        #     fluxN[i]=fluxN[i]*1e-22
        
        ### Get the total integrals
        if verbose is True:
            print(f"\n### TOTAL WALL INTEGRALS------------")
            print(f"\n### Heat flows:")
            for col in df.columns:
                if "E_" in col:
                    print(f"{col}: ---- {(df[col]*df['Area'] * 1e-6).sum():.2f} [MW]")
            
            print(f"\n### Particle flows:")      
            for col in df.columns:
                if "F_" in col:
                    print(f"{col}: ---- {(df[col]*df['Area']).sum():.3e} [s-1]")
        
        wfluxes_all = df.copy()
                    
        if split is True:
            if grid != "tight": raise Exception(f"Grid type {grid} not implemented yet")
            
            wfluxes = dict()
            df = wfluxes_all.copy()
            wfluxes["inner_lower_target"] = df.iloc[slice(22,44),:]
            wfluxes["lower_pfr"] = df.iloc[slice(44,74),:]
            wfluxes["outer_lower_target"] = df.iloc[slice(74,97),:]
            wfluxes["outer_wall"] = df.iloc[slice(97,156),:]
            wfluxes["outer_upper_target"] = df.iloc[slice(156,176),:]
            wfluxes["upper_pfr"] = df.iloc[slice(176,209),:]
            wfluxes["inner_upper_target"] = df.iloc[slice(209,231),:]
            wfluxes["inner_wall"] = df.iloc[np.r_[slice(0,22), slice(231,253)],:]
            
            if verbose is True:
                fig, ax = plt.subplots(dpi = 170)

                for region in wfluxes.keys():

                    ax.plot(df["R"], df["Z"], c = "k", alpha = 1, lw = 0, marker = "o", markersize = 2, markeredgewidth=0.1, markerfacecolor="None")
                    if region != "wall":
                        ax.plot(wfluxes[region]["R"], wfluxes[region]["Z"],  alpha = 0.5, lw = 0,  label = region, markersize = 4, marker = "o")
                        
                # Plot duplicates if any
                allregions = pd.concat(wfluxes.values())
                dupl = allregions[allregions.duplicated(subset="R")]
                if len(dupl) > 0:
                    ax.scatter(dupl["R"], dupl["Z"], c = "r", edgecolors="yellow", s = 50, label = "DUPLICATES", zorder = 200)
                else:
                    print("No duplicates found")
                    
                ax.set_aspect("equal")
                fig.legend(loc="upper left", bbox_to_anchor=(0.75, 0.7))
                ax.set_xlabel("R [m]")
                ax.set_ylabel("Z [m]")
                ax.set_ylim(-0.9, 0.9)
                
            # Return wall fluxes split by region
            return wfluxes

        else:
            # Return all wall fluxes
            return df
    
    
    def get_wall_ntmpi(self):
        """
        Get density, temperature, pressure of plasma and neutrals on wall
        return df
        """
        with HiddenPrints() if self.verbose is False else contextlib.nullcontext():  # Suppress all prints
            ## --------------- Get N, T, M, P on wall
            Eirene = self.Eirene
            RefPar = self.RefPar
            ions = self.ions
            
            # EIRENE results are here on the full tri mesh
            eirene_neutrals = h5py.File(os.path.join(self.path, "Results", "eirene_neutrals"), 'r')
            
            # Get wall ids. iWallKnots are the indices for the plasma parameters interpolated to tri
            # iWallTriangles are the indices for EIRENE parameters on full tri grid.
            ZeroTriangle, ZeroSide, Ri, Zi, RWallTriangles, ZWallTriangles, iWallTriangles, iWallSide, iWallKnots  = \
                    get_wall_triangle(Eirene, rz0_line=[0,0], theta_line=0, no_plot=1, no_print=1, no_triangles=0)

            RWallTriangles = np.append(RWallTriangles, RWallTriangles[0])
            ZWallTriangles = np.append(ZWallTriangles, ZWallTriangles[0])
            iWallKnots	   = np.append(iWallKnots, iWallKnots[0])
            # iWallTriangles	   = np.append(iWallTriangles, iWallTriangles[0])

            WalldL  	   = np.sqrt((RWallTriangles[1:]-RWallTriangles[:-1])**2 + (ZWallTriangles[1:]-ZWallTriangles[:-1])**2)
            DistKnots	   = np.cumsum(np.append(0., WalldL))
            DistTriangles  = 0.5*(DistKnots[:-1]+DistKnots[1:])
                
            evolution = 0

            if(evolution == 0):
                base_plasma_name = os.path.join(self.path, "Results") 
            else:
                base_plasma_name = os.path.join(self.path, "/Evolution/{:d}_".format(evolution))

            if_plasma	= h5py.File(os.path.join(base_plasma_name, "plasma_0"), "r")
            Te			= h5_read(if_plasma,"triangles/temperature")*RefPar.T0eV
            if_plasma.close()

            if_plasma	= h5py.File(os.path.join(base_plasma_name,"plasma_1"), "r")
            Ti			= h5_read(if_plasma,"triangles/temperature")*RefPar.T0eV
            if_plasma.close()

            Te 	= Te[iWallKnots]
            Ti 	= Ti[iWallKnots]
            
            for iPlasma in range(len(ions)):
                try:
                    if_plasma	= h5py.File(os.path.join(base_plasma_name, "plasma_"+str(iPlasma)), "r")
                except:
                    raise Exception("Cannot read plasma file")

                temperature	= h5_read(if_plasma,"triangles/temperature")*RefPar.T0eV
                density		= h5_read(if_plasma,"triangles/density")*RefPar.n0
                velocity	= h5_read(if_plasma,"triangles/velocity")*RefPar.c0

                temperature = temperature[iWallKnots]
                density 	= density[iWallKnots]
                velocity 	= velocity[iWallKnots]

                Jsat		= np.abs(1.6022e-19*velocity*density)*1e-3											#eletronic charge
                M			= velocity/np.sqrt((Te+Ti)/2)*(np.sqrt(RefPar.T0eV)/RefPar.c0)	
                
                if(iPlasma > 0):
                    
                    try:
                        Sn_tri	= h5_read(if_plasma,"triangles/Sn")*RefPar.n0/RefPar.tau0
                        if(Sn_tri.max() - Sn_tri.min() > 0.):
                            Sn_tri = set_min_positive(Sn_tri)
                            Sn	   = Sn_tri[iWallTriangles]
                            Sn	   = 0.5*np.append(np.append(Sn[0]+Sn[-1], Sn[:-1]+Sn[1:]), Sn[0]+Sn[-1])
                        else: Sn = 0.
                        Sn_tri = 0
            #				Sn	= -h5_read(if_plasma,"triangles/Sn")*RefPar.n0/RefPar.tau0									#recombination
                    except:
                        raise Exception("Cannot read Sn")
                        Sn  = 0.

                    try:
                        Nn = eirene_neutrals["atomic_species"]["dens_1"][:]
                        Tn = eirene_neutrals["atomic_species"]["T_1"][:]
                        Nm = eirene_neutrals["molecular_species"]["dens_1"][:]
                        Tm = eirene_neutrals["molecular_species"]["T_1"][:]

                        Nn = Nn[iWallTriangles]
                        Nm = Nm[iWallTriangles]
                        Tn = Tn[iWallTriangles]
                        Tm = Tm[iWallTriangles]

                        Tn = set_min_positive(Tn)
                        Tm = set_min_positive(Tm)
                        Nn	= 0.5*np.append(np.append(Nn[0]+Nn[-1], Nn[:-1]+Nn[1:]), Nn[0]+Nn[-1])*RefPar.n0
                        Nm	= 0.5*np.append(np.append(Nm[0]+Nm[-1], Nm[:-1]+Nm[1:]), Nm[0]+Nm[-1])*RefPar.n0
                        Tn	= 0.5*np.append(np.append(Tn[0]+Tn[-1], Tn[:-1]+Tn[1:]), Tn[0]+Tn[-1])
                        Tm	= 0.5*np.append(np.append(Tm[0]+Tm[-1], Tm[:-1]+Tm[1:]), Tm[0]+Tm[-1])

                        Pn	= (Nn*Tn+Nm*Tm)*1.6e-19
                        if(Pn.max() - Pn.min() < 0.): Pn = 0.
                    except:
                        raise Exception("Cannot read neutral data")
                else:
                    Sn	= 0.; Pn	= 0.

            def do_some_append(x):
                return x[:-1]
            
            df = pd.DataFrame()
            
            # TODO resolve the length issue here
            # print("walldL", len(WalldL))
            # print("Rtriangles", len(RWallTriangles))
            # print("Nn", len(Nn))
            # print("iWallKnots", len(iWallKnots))
            
            df["Nn"] = Nn
            df["Nm"] = Nm
            df["Tn"] = Tn
            df["Tm"] = Tm
            df["Sn"] = Sn
            df["M"] = M
            df["Te"] = Te
            df["Ti"] = Ti
            df["Ne"] = density
            
            df["walldL"] = np.append(WalldL, WalldL[-1])
            df["L"] = np.cumsum(df["walldL"])
            # No idea what's happening with the lengths
            df["R"] = RWallTriangles #do_some_append(RWallTriangles)
            df["Z"] = ZWallTriangles # do_some_append(ZWallTriangles)
            df["iWallKnots"] = iWallKnots
            
        
        return df
    
    def get_wall_data_on_target(self, df, target):
        """
        Align data with Te peaks for each divertor
        """
        out = scipy.signal.find_peaks(df["Te"], prominence = 1)
        peakid = out[0]
        indices = {
            "inner_lower":peakid[0],
            "outer_lower":peakid[1],
            "outer_upper":peakid[2],
            "inner_upper":peakid[3]
        }
        target = "inner_lower"
        df["dist"] = df["L"] - df["L"][indices[target]]

        return df
    
    def _get_rz_sep(self):

        """
        Queries grid regions to obtain separatrix R and Z coordinates
        Extends to the equilibium boundary
        Requires a Config (grid/eqb) file with MagZones computed
        Returns lhs and rhs for each of the sides of the separatrix for easy contiguous line plots
        Those are dicts with keys R and Z for the respective coords
        """
        Config = self.Config
        sep = dict()
        sep["UOT"] = Config.MagZones[2].north
        sep["UIT"] = Config.MagZones[3].north
        sep["LIT"] = Config.MagZones[4].north
        sep["LOT"] = Config.MagZones[5].north

        sep["ISEP"] = Config.MagZones[0].north
        sep["OSEP"] = Config.MagZones[1].north

        lhs = dict()
        rhs = dict()
        for coord in ["R", "Z"]:
            # lhs[coord] = np.
            lhs[coord] = sep["LIT"].__dict__[coord]
            lhs[coord] = np.concatenate([lhs[coord], sep["ISEP"].__dict__[coord]])
            lhs[coord] = np.concatenate([lhs[coord], sep["UIT"].__dict__[coord]])
            rhs[coord] = sep["UOT"].__dict__[coord]
            rhs[coord] = np.concatenate([rhs[coord], sep["OSEP"].__dict__[coord]])
            rhs[coord] = np.concatenate([rhs[coord], sep["LOT"].__dict__[coord]])
            
        return lhs, rhs
    
    
    
class HiddenPrints:
    """
    Suppress printing by:
    with HiddenPrints():
        ...
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        


#===================================================================

def get_transp_parameter_on_mesh(Config, IntCell, iBall):

	Zones = Config.Zones

	Transp = np.empty((IntCell.shape[0]), dtype='f8')

	iZones = np.unique(IntCell[:,0])
	for iZone in iZones:
		index = np.where(IntCell[:,0] == iZone)[0]
		Transp[index] = Zones[iZone].Ballooning[iBall][IntCell[index,1],IntCell[index,2]]

	return Transp


#===================================================================

def get_plasma_parameter_on_mesh(Plasma, iPar, IntCell):

	Par = np.empty(IntCell.shape[0], dtype='f8')

	if(Plasma[0].Nx == Plasma[0].Values[iPar].shape[0]):	ijOff = 0
	else:													ijOff = 1
	iZones = np.unique(IntCell[:,0])
	for iZone in iZones:
		index = np.where(IntCell[:,0] == iZone)[0]
		Par[index] = Plasma[iZone].Values[iPar][IntCell[index,1]+ijOff,IntCell[index,2]+ijOff]

	return Par


#===================================================================

def get_flux_parameter_on_mesh(Fluxes, Config, iPar, IntCell):

	Zones		= Config.Zones

	Par = np.empty((IntCell.shape[0], 5), dtype='f8')

	iZones = np.unique(IntCell[:,0])
	for iZone in iZones:
		index = np.where(IntCell[:,0] == iZone)[0]
		ii = IntCell[index,1]
		jj = IntCell[index,2]
		Par[index,:4] = Fluxes[iZone].Values[iPar][ii,jj,:4]

		ii1			  = ii + 1
		jj1			  = jj + 1
		SideLen		  = np.sqrt((Zones[iZone].gridR[ii1,jj1] - Zones[iZone].gridR[ii1,jj])**2 + (Zones[iZone].gridZ[ii1,jj1] - Zones[iZone].gridZ[ii1,jj])**2)			#North
		TotLen		  = np.copy(SideLen)
		Par[index,4]  = -Par[index,0]*SideLen
		SideLen		  = np.sqrt((Zones[iZone].gridR[ii,jj1]  - Zones[iZone].gridR[ii,jj])**2  + (Zones[iZone].gridZ[ii,jj1]  - Zones[iZone].gridZ[ii,jj])**2)			#South
		TotLen		 += SideLen
		Par[index,4] += Par[index,1]*SideLen
		SideLen		  = np.sqrt((Zones[iZone].gridR[ii1,jj1] - Zones[iZone].gridR[ii,jj1])**2 + (Zones[iZone].gridZ[ii1,jj1] - Zones[iZone].gridZ[ii,jj1])**2)			#East
		TotLen		 += SideLen
		Par[index,4] -= Par[index,2]*SideLen
		SideLen		  = np.sqrt((Zones[iZone].gridR[ii1,jj]  - Zones[iZone].gridR[ii,jj])**2  + (Zones[iZone].gridZ[ii1,jj]  - Zones[iZone].gridZ[ii,jj])**2)			#West
		TotLen		 += SideLen
		Par[index,4] += Par[index,3]*SideLen
		Par[index,4] /= TotLen

	return Par


def get_flux_parameter_on_flux_surfaces(Fluxes, Config, iPar, IntCell, RefPar):

	Zones		= Config.Zones
	Megazones	= Config.Megazones

	Par = np.zeros((IntCell.shape[0], 3), dtype='f8')

	iZones = np.unique(IntCell[:,0])
	for iZone in iZones:

		index = np.where(IntCell[:,0] == iZone)[0]
		ii = IntCell[index,1]
		ii1 = ii + 1

		for kZone in Megazones[Zones[iZone].mz].list[:]:
			LenSouth2R	  = np.sqrt((Zones[kZone].gridR[ii ,1:] - Zones[kZone].gridR[ii, :-1])**2 + (Zones[kZone].gridZ[ii ,1:] - Zones[kZone].gridZ[ii, :-1])**2)* \
								    (Zones[kZone].gridR[ii ,1:] + Zones[kZone].gridR[ii, :-1])*(1.-Zones[kZone].Chi[ii,:])
			LenNorth2R	  = np.sqrt((Zones[kZone].gridR[ii1,1:] - Zones[kZone].gridR[ii1,:-1])**2 + (Zones[kZone].gridZ[ii1,1:] - Zones[kZone].gridZ[ii1,:-1])**2)* \
								    (Zones[kZone].gridR[ii1,1:] + Zones[kZone].gridR[ii1,:-1])*(1.-Zones[kZone].Chi[ii,:])

			Par[index,0] += np.sum(Fluxes[kZone].Values[iPar][ii,:,0]*LenNorth2R,axis=1)				#North
			Par[index,1] += np.sum(Fluxes[kZone].Values[iPar][ii,:,1]*LenSouth2R,axis=1)				#South

			InOut	 = Zones[kZone].Chi[ii,1:] - Zones[kZone].Chi[ii,:-1]
			InOutMax = np.max(InOut, axis=1)															#East jump up
			InOutMin = np.min(InOut, axis=1)															#West jump down
			
			if(InOutMax.max() > 0):
				iigt, jgt = np.where(InOut > 0)															#Index jumps plasma to wall at East
				igt		 = ii[iigt]
				igt1 	 = igt + 1
				jgt1	 = jgt + 1
				LenEast	 = np.sqrt((Zones[kZone].gridR[igt1,jgt1] - Zones[kZone].gridR[igt,jgt1])**2 + (Zones[kZone].gridZ[igt1,jgt1] - Zones[kZone].gridZ[igt,jgt1])**2)* \
							       (Zones[kZone].gridR[igt1,jgt1] + Zones[kZone].gridR[igt,jgt1])
				Par[index[iigt],2] -= Fluxes[kZone].Values[iPar][igt,jgt,2]*LenEast

			if(InOutMin.min() < 0):
				iilt, jlt = np.where(InOut < 0)															#Index jumps  dowm wall to plasma West
				ilt		  = ii[iilt]
				ilt1 	  = ilt + 1
				jlt		  = jlt + 1																			#move to first cell with chi = 0
				LenEast	  = np.sqrt((Zones[iZone].gridR[ilt1,jlt] - Zones[iZone].gridR[ilt,jlt])**2 + (Zones[iZone].gridZ[ilt1,jlt] - Zones[iZone].gridZ[ilt,jlt])**2)* \
								    (Zones[iZone].gridR[ilt1,jlt] + Zones[iZone].gridR[ilt,jlt])
				Par[index[iilt],2] += Fluxes[kZone].Values[iPar][ilt,jlt,3]*LenEast
		Par[index,2] += Par[index,1] - Par[index,0]


	Par *= np.pi		#2*pi*R (alredy multiply by 2*R)

	return Par

#===================================================================

def get_plasma_parameter_on_pol(Plasma, iPar, ix, iZones, iThWest, iThEast, nThetaPts):

	Par = np.zeros(nThetaPts, dtype='f8')

	if(iPar > -1):
		jOff = 0
		for k in range(len(iZones)):
			if(Plasma[iZones[k]].Nz < Plasma[iZones[k]].Values[iPar].shape[1]):
				Par[jOff: jOff + iThEast[k] - iThWest[k]] = Plasma[iZones[k]].Values[iPar][ix, iThWest[k]+1:iThEast[k]+1]		#[Nx+2,Nz+2] parameter
			else:
				Par[jOff: jOff + iThEast[k] - iThWest[k]] = Plasma[iZones[k]].Values[iPar][ix, iThWest[k]:iThEast[k]]				#[Nx,Nz] parameter

			jOff += iThEast[k] - iThWest[k]

	return Par

#
# This routine set to the minimum positive value all negative values
#

def set_min_positive(Values):
	index_n = np.where(Values < 0.)
	if(len(index_n[0]) > 0):
		index_p = np.where(Values >= 0.)
		min_p   = np.min(Values[index_p])
		Values[index_n] = min_p

	return Values


def indexes(list_strs, str):
	return  [i for i, j in enumerate(list_strs) if j == str]


def create_norm(logscale, norm, vmin, vmax):
    if logscale:
        if norm is not None:
            raise ValueError(
                "norm and logscale cannot both be passed at the same time."
            )
        if vmin * vmax > 0:
            # vmin and vmax have the same sign, so can use standard log-scale
            norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            # vmin and vmax have opposite signs, so use symmetrical logarithmic scale
            if not isinstance(logscale, bool):
                linear_scale = logscale
            else:
                linear_scale = 1.0e-5
            linear_threshold = min(abs(vmin), abs(vmax)) * linear_scale
            if linear_threshold == 0:
                linear_threshold = 1e-4 * vmax   # prevents crash on "Linthresh must be positive"
            norm = mpl.colors.SymLogNorm(linear_threshold, vmin=vmin, vmax=vmax)
    elif norm is None:
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    return norm