import matplotlib.pyplot as plt
import csv
import numpy as np
import matplotlib as mpl


def read_target_C_MOD(filename):
		Dict = {}
		with open(filename, "r") as f:
		  for i, line in enumerate(f):
		        if i>0:
		            DATA = line.split(" ")
		            for g,el in enumerate(NameDATA):
		                  Dict[el] .append(float(DATA[g]))
		        else:
		            DATA = line.split(" ");lenlist = len(DATA);NameDATA = DATA;print(NameDATA)
		            for g,el in enumerate(NameDATA):
		                  list = [];Dict[el] = list
		return  Dict


def read_exp_C_mod(filename):
	DictI = {}
	with open(filename, "r") as f:
		for line in f:
			if "psinorm" in line:
				List0 = [];Dict = {}
				Ar0 = line.split(" ")[1:];	Arnumb = int(line.split(" ")[0])
				print(Arnumb)
				Arg1 = Ar0[0];Arg2 = Ar0[1];Arg3 = Ar0[2];
				for el in Ar0:
					Dict[el] = []
				DictI[Ar0[1]] = Dict
				i = 0
			if i<=Arnumb:
				if i>0:
				  Ar = line.split(" ");Ar= [enn for enn in Ar if enn != ""] #''.join(line.split())
				  Dirr = DictI[Ar0[1]] 
				  pp = list(Dirr.keys())
				  Dirr[Arg1].append(float(Ar[0]))
				  Dirr[Arg2].append(float(Ar[1]))
				  Dirr[Arg3].append(float(Ar[2]))
				i+=1
	return DictI

def read_file(Name):
        i = 0;Dict = {}
        with open(Name,"r") as f:
                        for line in f:
                          if i>0:
                                 A = line.split("\t");B = [float(el) for el in A if el !=""]
                                 for imk,er in enumerate(List_name):
                                         Dict[er].append(B[imk])
                                 i+=1
                          else:
                                 A = line.split("\t");List_name = [el for el in A if el !=""]
                                 for er in List_name:
                                         Listt =[];Dict[er] = Listt
                                 i+=1
        return Dict


def read_csv(Name,delimiter = ",",flag = False,flag_count = False):
    with open(Name) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=delimiter)
        i = 0; Dict = {};List_name2 =[]
        for column in readCSV:
            if i  == 0:
              List_name = column
#              print(List_name)
              for er in List_name:
                if er in Dict.keys():
                    Listt =[];Dict[str(er+'_0'+str(i))] = Listt
                    List_name2.append(str(er+'_0'+str(i)))
                else:
                    List_name2.append(str(er))
                    Listt =[];Dict[str(er)] = Listt
            if i>0:
                for imk,er in enumerate(List_name2):
                      if flag == True:
                          if imk< len(List_name2)-1:
                             Dict[er].append(float(column[imk]))
                      else:
                          Dict[er].append(float(column[imk]))
            i+=1
    #for elem in Dict.keys():
    #    Listi = Dict[elem];
    #    mylist = list( dict.fromkeys(Listi));Dict[elem] = mylist
    return Dict








           
        

