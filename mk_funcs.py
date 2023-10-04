import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys, pathlib
import platform
import traceback
import xarray
import xbout
import scipy
import re
import netCDF4 as nc

onedrive_path = onedrive_path = str(os.getcwd()).split("OneDrive")[0] + "OneDrive"
sys.path.append(os.path.join(onedrive_path, r"Project\python-packages\sdtools"))
sys.path.append(os.path.join(onedrive_path, r"Project\python-packages\soledge"))
sys.path.append(os.path.join(onedrive_path, r"Project\python-packages"))


from gridtools.b2_tools import *
from gridtools.utils import *

from hermes3.fluxes import *
from hermes3.case_db import *
from hermes3.load import *
from hermes3.named_selections import *
from hermes3.plotting import *
from hermes3.grid_fields import *
from hermes3.accessors import *
from hermes3.utils import *
from code_comparison.viewer_2d import *
from code_comparison.code_comparison import *

from gridtools.solps_python_scripts.read_b2fgmtry import *

from soledge.plot1d_wall_fluxes import plot1d_wall_fluxes



def get_soledge_wall_fluxes(
    path_solps,
    path_soledge,
    plot = False
    ):
    """ 
    Compute soledge wall fluxes on a per region basis, where the regions are obtained by 
    interrogating a SOLPS case
    """

    # Open SOLPS results and set the boundaries of each target

    with nc.Dataset(os.path.join(path_solps, "balance.nc")) as d:
        crx = d["crx"]
        cry = d["cry"]

    g = read_b2fgmtry(where=path_solps)
    p = SOLPSplot(path_solps, "ti")

    if plot is True:
        plt.close("all")
        fig, ax = plt.subplots(dpi = 200, figsize = (8,4))
        p.plot(fig = fig, ax = ax, antialias = True, linewidth = 0.1)
        ax.set_ylim(-1,1)

    upper_break = int(np.mean([g["rightcut"][1], g["leftcut"][0]]))

    bounds = {
        # It' [corner, radial(logical y), poloidal(logical x)]
        "inner_lower_inner_corner" : [1, 0, 0],
        "inner_lower_outer_corner" : [1, -1, 0],
        "inner_upper_outer_corner" : [1, -1, upper_break + 10],
        "inner_upper_inner_corner" : [1, 0, upper_break + 10],
        "outer_upper_inner_corner" : [1, 0, upper_break + 11],
        "outer_upper_outer_corner" : [1, -1, upper_break + 11],
        "outer_lower_outer_corner" : [1, -1, -1],
        "outer_lower_inner_corner" : [1, 0, -1],
    }
    rz_extents = dict()

    # Plot to check - some of the extents were fudged manually
    for name in bounds:
        p1 = bounds[name]
        rz_extents[name] = dict()
        rz_extents[name]["R"] = crx[p1[0], p1[1], p1[2]]
        rz_extents[name]["Z"] = cry[p1[0], p1[1], p1[2]]
        
        if plot is True:
            ax.scatter(crx[p1[0], p1[1], p1[2]], cry[p1[0], p1[1], p1[2]], label = name)
            
    # Collect SOLEDGE wall fluxes. I modified the routine to output a dataframe with the results
    
    wfluxes_all = dict()
    wfluxes_all = plot1d_wall_fluxes(path_soledge, no_plot = 1)

    wfluxes = dict()
    df = wfluxes_all.copy()

    # Filtering the wall fluxes by each region
    wfluxes["inner_lower"] = df.query(f"Z < {rz_extents['inner_lower_outer_corner']['Z']} & R < {rz_extents['inner_lower_inner_corner']['R']}")
    wfluxes["inner_wall"] = df.query(f"Z > {rz_extents['inner_lower_outer_corner']['Z']} & R < {rz_extents['inner_upper_outer_corner']['R']} & Z < {rz_extents['inner_upper_outer_corner']['Z']}")
    wfluxes["inner_upper"] = df.query(f"Z > {rz_extents['inner_upper_outer_corner']['Z']} & R < {rz_extents['inner_upper_inner_corner']['R']}")
    wfluxes["upper_pfr"] = df.query(f"R > {rz_extents['inner_upper_inner_corner']['R']} & R < {rz_extents['outer_upper_inner_corner']['R']} & Z > 0")
    wfluxes["outer_upper"] = df.query(f"R > {rz_extents['outer_upper_inner_corner']['R']} & R < {rz_extents['outer_upper_outer_corner']['R']} & Z > {rz_extents['outer_upper_inner_corner']['Z']}")
    wfluxes["lower_pfr"] = df.query(f"R > {rz_extents['inner_upper_inner_corner']['R']} & R < {rz_extents['outer_upper_inner_corner']['R']} & Z < 0")
    wfluxes["outer_lower"] = df.query(f"R > {rz_extents['outer_lower_inner_corner']['R']} & R < {rz_extents['outer_lower_outer_corner']['R']} & Z < {rz_extents['inner_lower_inner_corner']['Z']}")

    # Now going by index. Note that SOLEDGE goes anticlockwise, but we go clockwise...
    wfluxes["upper_chamber"] = df.query(f"index < {wfluxes['outer_upper'].index[0]} & index > {wfluxes['outer_upper'].index[0]} - 45")
    wfluxes["upper_baffle"] = df.query(f"index < {wfluxes['upper_chamber'].index[0]} & index > {wfluxes['upper_chamber'].index[0]} - 16")

    wfluxes["lower_chamber"] = df.query(f"index > {wfluxes['outer_lower'].index[-1]} & index < {wfluxes['outer_lower'].index[-1]} + 45 + 2")    # Asymmetric offset due to oddities in SOLEDGE grid
    wfluxes["lower_baffle"] = df.query(f"index > {wfluxes['lower_chamber'].index[-1]} & index < {wfluxes['lower_chamber'].index[-1]} + 16 - 2")

    # Index starts at OMP, so there is index discontinuity in the middle
    wfluxes["outer_lower_wall"] = df.query(f"index > {wfluxes['lower_baffle'].index[-1]}")# & index < {wfluxes['upper_baffle'].index[0]}")
    wfluxes["outer_upper_wall"] = df.query(f"index < {wfluxes['upper_baffle'].index[0]}")
        
    # Now plot the SOLEDGE region splits to check that they make sense
    if plot is True:
        fig, ax = plt.subplots(dpi = 150)
        for region in wfluxes.keys():
            ax.plot(df["R"], df["Z"], c = "k", alpha = 1, lw = 0, marker = "o", markersize = 2, markeredgewidth=0.1, markerfacecolor="None")
            if region != "wall":
                ax.plot(wfluxes[region]["R"], wfluxes[region]["Z"],  alpha = 0.5, lw = 0,  label = region, markersize = 4, marker = "o")
        for i, name in enumerate(rz_extents.keys()):
            label = "SOLPS target boundary" if i == 0 else ""
            ax.scatter(rz_extents[name]["R"], rz_extents[name]["Z"], marker = "x", c = "deeppink", label = label, zorder = 100)
            
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
        
    # Now prepare wall region integrals
    import copy

    df = pd.DataFrame()
    wfluxes_integral = copy.deepcopy(wfluxes)

    params = [x for x in wfluxes_all.columns if x not in ["R", "Z", "L"]]

    for region in wfluxes_integral.keys():
        for param in params:
            wfluxes_integral[region][param] = wfluxes_integral[region][param] * 2*np.pi * wfluxes_integral[region]["R"] * wfluxes[region]["dL"]
            # if "hflux_par_tot" not in param: # This is a scalar
            df.loc[region, param] = wfluxes_integral[region][param].sum()

    df = wfluxes_integral

    return wfluxes, wfluxes_integral