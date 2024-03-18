#!/usr/bin/env python

# Author: Niamh Smith, e-mail: niamh.smith.17 [at] ucl.ac.uk
# Date:   13-03-2024; 18-03-2024 [updated]

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

class adiabatic_potentials:
    """ 
        Plotting and deriving relaxation energies from an adiabatic potential curve (total energy Vs reaction coordinate) of a defect in two charge states q1 and q2. Reaction coordinate, Q, is taken to be the vector multiplier of the xyz column vectors which describe the difference in atomic positions between relaxed geometries of q1 and q2. 
        Method follows that described by "6. Nonradiative multiphonon transition rates, T. Grasser, Microelectronics Reliability 52 (2012) 39–70"
    
        Class definitions:
            energy(str): string within CP2K log file in which final total energy is found. 
            at_structure(str): string within CP2K log file which signals that printed atom coordinates used in calculation follow.  
            Bandgap(float): global variables E_VBM and E_CBM (defined outside of class) must be the eigenvalue energies (in eV) of the bulk defect free material's highest occupied and lowest unoccupied molecular orbits, respectively. 
    
        Inputs:
            Ename(str): Name of band edge energy defect is closest to. 
            Ebm(float): energy of band edge defect is closest to. 
            E_q2(float): final energy of relaxed charge state q2 structure
            x_q2(list of floats): x-coordinates of relaxed charge state q2 structure
            y_q2(list of floats): y-coordinates of relaxed charge state q2 structure
            z_q2(list of floats): z-coordinates of relaxed charge state q2 structure
            E_q1(float): final energy of relaxed charge state q1 structure
            x_q1(list of floats): x-coordinates of relaxed charge state q1 structure
            y_q1(list of floats): y-coordinates of relaxed charge state q1 structure
            z_q1(list of floats): z-coordinates of relaxed charge state q1 structure
            CTL(float): charge transition level calculated between charge states q1 and q2 wrt the VBM
            q2_rsk_log(str): file path of log file for ENERGY calc of q1 geometry in charge state q2.
            q1_rsk_log(str): file path of log file for ENERGY calc of q2 geometry in q1 charge state.
            
        Attributes [self. variables]:
            name(str): assigned as name of the defect within function self.renew file path. 
            Ename(str): Ename input saved as class attribute
            Ebm(float): Ebm input saves as class attribute
            CTL(float): CTL input saved as class attribute
            E_q2min_pos(None -> int), E_q1min_pos(None -> int): Q value minimum of charge state parabola, assigned by function self.comparing
            M(None -> float): associated modal mass assigned by function self.find_each_atom_vector.
            fitting_cst1(None -> list of floats), fitting_cst2(None -> list of floats): curve fit [scipy.optimize] coefficient to fit self.objective eqn to each charge state results. 
            E_q2_xyz(array of floats), E_q1_xyz(array of floats): 3 x n ndarray, where n is number of atoms, of x, y, z relaxed coord inputs for charge state. 
            ε12(None -> float), ε21(None -> float): theory SRH levels comparible with experimental extracted trap/interface levels, assigned by self.CTL_2_exp_comparison()
            mass_1(list of float), mass_2(list of float): atomic masses of all atoms in system. 
            elements(list of str): atomic species names for all atom in system. 
            E_q2(float): E_q2 input plus Ebm input (energy adjustment wrt closest band egde) saved as class attribute
            E_q1(float): E_q1 input saved as class attribute. 
            E_q2_rsk(float): energy returned from function self.variables_from_logs + E_bm input (energy adjustment wrt closest band egde) saved as class attribute
            E_q1_rsk(float): energy returned from function self.variables_from_logs saved as class attribute
            Es_q2(list), Es_q1(list): energy points to be plotted for charge state 
            Qs_q2(list), Qs_q1(list): reaction coordinate points to be plotted for charge state. 
    """
    
    energy = 'ENERGY| Total FORCE_EVAL'
    at_structure = 'MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom'
    
    def __init__(self, Ename, Ebm, E_q2, x_q2, y_q2, z_q2, E_q1, x_q1, y_q1, z_q1, q2_rsk_log, q1_rsk_log):
        # check that both q1 & q2 have coordinates for the same number of atoms.
        
        if len(x_q2) != len(x_q1): 
            raise ValueError("length of the coordinates of q2 don't match length of coordinates of q1")
        
        # Set up class attributes from input variables 
        self.name = "" 
        self.Ename, self.Ebm, self.CTL, self.E_q2min_pos, self.E_q1min_pos, self.M, self.fitting_cst1, self.fitting_cst2, self.ε12, self.ε21 = Ename, Ebm, CTL, None, None, None, None, None, None, None
        self.E_q2_xyz, self.E_q1_xyz = np.array([x_q2,y_q2,z_q2]), np.array([x_q1,y_q1,z_q1])
        self.Es_q2, self.Qs_q2, self.Es_q1, self.Qs_q1  = [],[], [], []
        self.E_q2, self.E_q1 = E_q2+Ebm, E_q1
        
        # search the CP2K log files for energies, atom names, and atom masses. 
        E_q2_rsk, element, self.mass_1 = self.variables_from_logs("q2*", "q1", q2_rsk_log, self.arrays_(x_q1,y_q1,z_q1))
        E_q1_rsk, self.element, self.mass_2 = self.variables_from_logs("q1*", "q2", q1_rsk_log, self.arrays_(x_q2, y_q2, z_q2))
        self.E_q2_rsk, self.E_q1_rsk =  E_q2_rsk+Ebm, E_q1_rsk
        
        # setting up coordinate files for different Q values and extracting energies from the ENERGY calc log files of these geometries
        self.deriving_Q_coord_xyzs(self.E_q1_xyz, self.E_q2_xyz, self.renew_file_path(q2_rsk_log))
        
        self.plotting(self.Es_q2, self.Qs_q2, self.Es_q1, self.Qs_q1)
        
    def arrays_(self, x, y, z):
        # each coordinate element within the x,y,z lists rounded to 6 decimal places using function roundlistelements (outside class). 
        for j in range(len(x)):
            x, y, z = roundlistelements(x,j,6), roundlistelements(y,j,6), roundlistelements(z,j,6) # Lists updated for each value of j
        xyz = np.array([x,y,z]) 
        
        return xyz
        
    def variables_from_logs(self, strg1, strg2, rsk_log, compare):
        """
            Search the CP2K log files for energies, atom names, and atom masses.
            
            inputs:
                strg1(str): name of excited charge state which is in the relaxed geometry of strg2
                strg2(str): name of charge state from relaxed coordinates were derived.
                rsk_log(str): path of log file of CP2K ENERGY calc for excited charge state strg1 with geometry of strg2
                compare(array): array of strg2 coordinates 
                
            outputs:
                E(float): total energy calculated from ENERGY calculation of excited charge extracted from log file
                ELEMENT(list of str): element names of all atoms within the system extracted from log file
                MASS(list of floats): atomic masses of all atom within the system extracted from log file
        """
        
        log_rsk = open(rsk_log, 'r') # open up log file and set variables for data extraction including lists for extracted data storage     
        lines, index, index2, X, Y, Z, ELEMENT, MASS = log_rsk.readlines(), 0, len(log_rsk.readlines()) + 1, [], [], [], [], []
        lines2 = lines
        vars()["%sfound"%'at_structure'] = False # boolean set up to identify when certain strings have been found in log file
        
        for line in lines: # starting at top of log file
            index += 1
            if adiabatic_potentials.at_structure in line and vars()["%sfound"%'at_structure'] is False:
                for p in range(0,[int(var[1]) for var in enumerate(lines[index-20].split()) if var[0] in [3]][0]): # for p from 0 to total amount of atoms reached
                    element, x, y, z, mass = [var[1] for var in enumerate(lines[index+2+p].split()) if var[0] in [2,4,5,6,8]] # data extraction
                    for f, F in zip([float(x),float(y),float(z),element, mass],[X,Y,Z,ELEMENT, MASS]):
                        F.append(f) # saving extracted data witihn predefined lists.
                array = np.array([X,Y,Z])
                vars()["%sfound"%'at_structure'] = None # changing boolean to None to break out of if loop.
                
        self.comparing(array, compare, strg1, strg2)
        
        E = self.log_energy(lines2, index2)
                
        log_rsk.close()
        return E, ELEMENT, MASS
    
     def comparing(self, array, compare, strg1, strg2):
        
        if type(compare) is np.ndarray:
            """
                comparing to check coordinates from strg1 excited charge state log file are same as coordinates of relaxed geometry coordinates of strg2 charge state. 
            """
            for d in range(len(array)): 
                for i in range(len(array[d])):
                    if array[d][i] != compare[d][i]: # error if array of excited charge state geometry doesn't match array of relaxed state geometry.
                        raise ValueError("for coordinate",d,i,"atomic positions array created from log file of",strg1,"does not match the positions of relaxed",strg2,"charge state, with",strg1,"coordinate of",array[d][i],"and",strg2,"coordinate of",compare[d][i])
        else:
            """ when input compare is given as either 0 or 1. 
            """
            # boolean variables to be change from None to False when a following condition is met. 
            q2_ = None
            q1_ = None
            # putting total energies, extracted from log files, in to context wrt to reaction coordinates (either 0 or 1).  
            for d in range(len(array)):
                for i in range(len(array[d])):
                    if array[d][i] != self.E_q2_xyz[d][i]:
                        q2_ = False
                    elif array[d][i] != self.E_q1_xyz[d][i]:
                        q1_ = False
                    if i == len(array[d]) - 1 and d == len(array) - 1 : # defining particular Q position of minima
                        if q2_ == None and q1_ == False:
                            self.E_q2min_pos = compare 
                        elif q1_ == None and q2_ == False:
                            self.E_q1min_pos = compare
                        elif q2_ == None and q1_ == None or q2_ == False and q1_ == False:
                            raise ValueError("something is wrong")
        
    def log_energy(self, lines2, index2):
        vars()["%sfound"%'En'] = False # boolean set up to identify when certain strings have been found in log file
        for line in reversed(lines2): # starting from bottom of log file 
            index2 -= 1 
            if adiabatic_potentials.energy in line and vars()["%sfound"%'En'] is False:
                E = [(round(float(var[1]),10) * 27.211 ) for var in enumerate(line.split()) if var[0] in [8]][0]
                vars()["%sfound"%'En'] = None
                
        if vars()["%sfound"%'En'] == False: # changing boolean to None to break out of if loop.
            raise ValueError("Final total energy not found")
        
        return E
   
    def deriving_Q_coord_xyzs(self, xyz0, xyz1, path):
        """
            inputs:
                xyz0(array): 3 x n ndarray, where n is number of atoms, of x, y, z coordinates of relaxed charge state q1
                xyz1(array): 3 x n ndarray, where n is number of atoms, of x, y, z coordinates of relaxed charge state q2
                path(str): path of directory where geometry and log CP2K files for different reaction coordinates are expected to be found. 
        """
        # finding minimums of each charge state parabola 
        self.comparing(xyz0,0,"xyz0","")
        self.comparing(xyz1,1,"xyz1","")
        
        # defining column vector difference between geometries of charge states.
        vectors = self.find_each_atom_vector(xyz0, xyz1)
        
        num = len(xyz0[0]) 
        
        for Q in [round(Q_,1) for Q_ in np.arange(-2,3.2,0.20)]: # rounded to 1 dp so you don't get weird values of Q such as 2.1333335
            flname = str("Q_{}_{}_pos.xyz".format(abs(Q),self.name)) if Q <= 0 else str("Q{}_{}_pos.xyz".format(Q,self.name)) # defining geometry file name for each Q, for negative Q values the '-' is changed to a '_'.
            
            if indirectory(path, flname) == False: # if Q geometry file name not in directory
                file = str("{}Q_{}_{}_pos.xyz".format(path,abs(Q),self.name)) if Q <= 0 else str("{}Q{}_{}_pos.xyz".format(path,Q,self.name)) # including path with file name
                print(file) # to show this file is being produced. 
                first = str("     {}\n \n".format(num)) # first line of CP2K xyz file needs number of atoms stated
                with open(file, 'w') as  f:
                    f.write(first) # write first line to file. 
                    
                for l in range(num): # for each atom. 
                    x_coord, y_coord, z_coord = "{0:.8f}".format(xyz0[0][l] + (Q * vectors[l][0])), "{0:.8f}".format(xyz0[1][l] + (Q * vectors[l][1])), "{0:.8f}".format(xyz0[2][l] + (Q * vectors[l][2])) # setting up strings of each coordinate where the resultant coordinate is the q1 coordinate + (column vector component multipled by Q) 
                    
                    if round(Q,2)== 1.00: # for when Q is 1, check method and column vector with coordinates of xyz1. 
                        if xyz0[0][l] + (Q * vectors[l][0]) != xyz1[0][l] or xyz0[1][l] + (Q * vectors[l][1]) != xyz1[1][l] or xyz0[2][l] + (Q * vectors[l][2]) != xyz1[2][l]:
                            raise ValueError("something has gone wrong with vectors 1")
                        if xyz1[0][l] - (Q * vectors[l][0]) != xyz0[0][l] or xyz1[1][l] - (Q * vectors[l][1]) != xyz0[1][l] or xyz1[2][l] - (Q * vectors[l][2]) != xyz0[2][l]:
                            raise ValueError("something has gone wrong with vectors 2")
                            

                    line = str(self.element[l]).rjust(2)+str(x_coord).rjust(19)+str(y_coord).rjust(19)+str(z_coord).rjust(19)+"\n" # atom coordinate line
                    with open(file, 'a') as  f:
                        f.write(line) # write line to file. 
                        
            else: # if Q geometry file in directory. 
                if Q == -0.0 or Q == 1.0: # energies already extracted for Q = 0 and Q = 1.
                    if abs(Q) == self.E_q2min_pos: # if energy minimum of q2 is at Q then so is the energy of excited q1 
                        self.Es_q2.append(self.E_q2)
                        self.Qs_q2.append(abs(Q)) # abs used here because np.arange returns -0.0 not 0.0 
                        self.Es_q1.append(self.E_q1_rsk)
                        self.Qs_q1.append(abs(Q)) 
                    elif abs(Q) == self.E_q1min_pos: # if energy minimum of q1 is at Q then so is the energy of excited q2
                        self.Es_q2.append(self.E_q2_rsk)
                        self.Qs_q2.append(abs(Q)) 
                        self.Es_q1.append(self.E_q1)
                        self.Qs_q1.append(abs(Q)) 

                elif Q > 2.0: # Energy calculations for Q values greater than 2 were only performed on the neutral (q2) charge state. 
                    flname = str("Neu_Q{}_ENERGY.log".format(Q))
                    if indirectory(path, flname) == True:
                        file = str("{}{}".format(path,flname))
                        log = open(file, 'r')        
                        lines, index = log.readlines(),len(log.readlines()) + 1
                        E = self.log_energy(lines, index) # extract energy value
                        log.close()
                        E = E + self.Ebm # add energy adjustment wrt closest band egde. 
                        self.Es_q2.append(E)
                        self.Qs_q2.append(Q)
                elif Q in [-0.6,-0.2,0.2,0.6,1.4,1.8]: # don't want energy values from these values of Q. 
                    w = 1
                else:
                    for N, D in zip(['Neu', 'Neg'],['q2','q1']):
                        flname = str("{}_Q_{}_ENERGY.log".format(N,abs(Q))) if Q <= 0 else str("{}_Q{}_ENERGY.log".format(N,Q)) # name of CP2K log file for Q ENERGY calculations, for negative Q values the '-' is changed to a '_'
                        if indirectory(path, flname) == True: 
                            file = str("{}{}".format(path,flname))
                            log = open(file, 'r')        
                            lines, index = log.readlines(),len(log.readlines()) + 1
                            E = self.log_energy(lines, index) # extract energy value
                            log.close()
                            if D == 'q2': # if ENERGY value calculated for q2, add energy adjustment wrt closest band egde. 
                                E = E + self.Ebm
                                self.Es_q2.append(E)
                                self.Qs_q2.append(Q)
                            else:
                                self.Es_q1.append(E)
                                self.Qs_q1.append(Q)
                                
    def find_each_atom_vector(self, xyz0, xyz1):   
        """
            inputs:
                xyz0(np.ndarray): 3 x n ndarray, where n is number of atoms, of x, y, z coordinates of relaxed charge state q1
                xyz1(np.ndarray): 3 x n ndarray, where n is number of atoms, of x, y, z coordinates of relaxed charge state q2
                
            outputs:
                atom_vectors(np.ndarray): 3 x n ndarray, where n is number of atoms, of column vector  / x \ components for the displacement of atoms in
                                                                                                      |  y  |
                                                                                                       \ z /
                position xyz1 from xyz0
        """
        # create new array for each atoms column vector.
        atom_vectors = np.ndarray(shape = (len(xyz0[0]), 3))
        
        for l in range(len(xyz0[0])): # calculating the difference in each dirn between charge state geometries for each atom
            dx, dy, dz = float(xyz1[0][l]) - float(xyz0[0][l]), float(xyz1[1][l]) - float(xyz0[1][l]), float(xyz1[2][l]) - float(xyz0[2][l])

            atom_vectors[l] = dx, dy, dz # 
            
            r = np.sqrt(dx**2 + dy**2 + dz**2) # magnitude of vector between charge state geometry positions for an atom
            
            if self.mass_2[l] != self.mass_1[l]: # same defect, different charge states - arrays of masses for each charge state must be the same 
                raise ValueError("You can't do it this way")
                
            M = float(self.mass_2[l]) * r**2 if l == 0 else M + float(self.mass_2[l]) * r**2 # calculation of the associated modal mass
            
        self.M = M 
                
        return atom_vectors
    
    def renew_file_path(self, path):
        """
            inputs:
                path(str): directory path of log file for ENERGY calc of q1 geometry in charge state q2
                
            outputs:
                Path(str): path of directory where geometry and log CP2K files for different reaction coordinates are expected to be found. 
        """
        Path = path.split('/') # make list of each directory in path from splitting path by the '/' within it. 
        self.name = Path[-3] # name is the 3rd to last item in list
        
        # remove and replace last item of list with directory where geometry and log CP2K files for different reaction coordinates are expected to be found
        Path.pop(-1)
        Path.append("adiabatic_potentials/")
        Path = '/'.join(Path) # join list into one big string and input '/' between list items. 
        return Path 
    
    def plotting(self,Es_q2, Qs_q2, Es_q1, Qs_q1): 
        """
            inputs:
                Es_q2(list), Es_q1(list): energy points to be plotted for charge state 
                Qs_q2(list), Qs_q1(list): reaction coordinate points to be plotted for charge state.
        """
        
        fig, ax = plt.subplots(figsize =(12,8))
        
        εR_12, εR_21, E_21, E_12 = self.energies()
        
        # big for loop to save lines related to line and point plotting for each charge state. 
        for q, q_, E, Eq, c1, c2, c3, c4, c5, c6, lab1, lab2, lab3, txt1, num, num2, n1, n2, n3, n4, n5, n6 in zip(['q2','q1'],['q1','q2'],['E2','E1'],[self.E_q2,self.E_q1], ['brown','seagreen'],['tomato','mediumaquamarine'],['red','cyan'],['crimson','blue'],['deeppink','violet'], ['orangered','teal'], [str('D0 + {}'.format(self.Ename)), 'Dn1'],["$E_{D0}$"+str(' + {}'.format(self.Ename)),'$E_{Dn1}$'],["$E_{D0*}$"+str(' + {}'.format(self.Ename)), '$E_{Dn1*}$'],["$\epsilon_{R_{12}}$ = %.3f eV"%εR_12, "$\epsilon_{R_{21}}$ = %.3f eV"%εR_21], [1,2], [2,1], [0.6,-2], [0.58,-0.6], [self.E_q2min_pos,-2.1], [3.05,0], [-0.8, -0.2], [0.01,0.02]):
            
            exec(f'ax.scatter(Qs_{q}, Es_{q}, marker="o", color=c1, label=lab1)') # scatter of energies and reaction coordinate points
            exec(f'ax.scatter(self.E_{q}min_pos, self.E_{q}, marker="*", s=100, color=c2,label=lab2)') # make minimum of charge state parabola stand out
            exec(f'ax.scatter(self.E_{q_}min_pos, self.E_{q}_rsk, marker="o", color=c2,label=lab3)') # make excited point of charge state stand out
            
            if len(eval("Qs_{}".format(q))) > 2: 
                fitting_csts = self.curvefitting(np.array(eval("Qs_{}".format(q))), np.array(eval("Es_{}".format(q))), eval("self.E_{}".format(q)), False)# fitting data to parabola defined in self.objective, with boolean input argument as false. 
                ax.plot(np.arange(-2.4, 3.4, 0.1), self.objective(np.arange(-2.4, 3.4, 0.1), fitting_csts[0], fitting_csts[1],fitting_csts[2]),'-.', color=c4,label='fit')
                exec(f'self.fitting_csts{num2} = fitting_csts')
                
            else: # for when using more than 2 points to create a parabola results in an avoided crossing. 
                exec(f'self.fitting_csts{num2} = False')

            if self.Ename == "$E_{VBM}$ ($α_{_{HOMO_{perf}}}$)": # defining how to plot arrows of ε˚12 and ε˚21 based on where the closest band edge to the defect state is the valence band. 
                op2 = "+"
                op1 = "-"
            else: # or the conduction band.
                op2 = "-"
                op1 = "+"
            
            exec(f'ax.arrow(self.E_{q}min_pos,self.E_{q},0,(self.E_{q_}_rsk-self.E_{q}{eval("op{}".format(num))}0.04),head_width=0.04,head_length=0.04,width=0.01, linewidth= 0.0001, edgecolor=c5,facecolor=c5,alpha =0.5)')
            
            # pointing out and stating energy value of charge state minimum. 
            ax.text(eval("self.E_{}min_pos".format(q))+n1,Eq-0.01,str("{} = ".format(E))+"%.3f eV"% Eq,color=c2,fontsize=14)
            exec(f'ax.plot([self.E_{q}min_pos,self.E_{q}min_pos+n2],[Eq,Eq],ls="-",color=c2,linewidth=1)')
            ax.plot([n3,-2.5],[Eq,Eq],ls='--',color=c2,linewidth=1)
            ax.plot([n4,3.4],[Eq,Eq],ls='--',color=c2,linewidth=1)
            
            ax.text(eval("self.E_{}min_pos".format(q))+n5,eval("self.E_{}_rsk".format(q))+n6,txt1,color=c6,fontsize=14) # Text corresponds to vertical energies of R21 and R12, place on top of dotted line parallel to x-axis, as no space for it next to the dash-dotted line parallel to y-axis. 
            exec(f'ax.plot([self.E_{q}min_pos,self.E_{q}min_pos],[self.E_{q}_rsk,Eq],ls="-.",color=c6,linewidth=1.5)') # energy difference between excited and ground state - R21/R12. 
            exec(f'ax.plot([self.E_{q_}min_pos,self.E_{q}min_pos],[self.E_{q}_rsk,self.E_{q}_rsk],ls=":",color=c6,linewidth=1,alpha=0.8)') # reaction coordinate difference between excited and ground state.
            
            
        for q, c4, num, num2 in zip(['q2','q1'],['crimson','blue'],[1,2], [2,1]):
            if eval("self.fitting_csts{}".format(num2)) is False and eval("self.fitting_csts{}".format(num)) is not False: # for when using more than 2 points to create a parabola results in an avoided crossing. 
                fitting_csts = self.curvefitting(np.array(eval("Qs_{}".format(q))), np.array(eval("Es_{}".format(q))), eval("self.E_{}".format(q)), True) # fitting data to parabola defined in self.objective, with boolean input argument true. 
                ax.plot(np.arange(-2.4, 3.4, 0.1), self.objective(np.arange(-2.4, 3.4, 0.1), fitting_csts[0], fitting_csts[1],fitting_csts[2]),'-.', color=c4,label='fit')
                exec(f'self.fitting_csts{num2} = fitting_csts')
        
        Qx, TEx, R_cap, R_emis, lwrE, opp, limit = self.relax_E()
        
        ax.plot([Qx,Qx],[TEx,TEx+R_cap],ls='--',color="orange") # vertical line from crossing point to lower energy minimum 
        ax.text((Qx+opp[0]), ((TEx+R_cap)+0.01), "$R_{cap}$ = %.3f eV"%(R_cap),color="orange",fontsize = 14) # stating relaxation energy for capture 
        ax.plot([Qx+0.02,Qx+0.02],[TEx,TEx+R_emis],ls='--',color="green") # vertical line from crossing point to higher energy minimum
        ax.text((Qx+opp[1]), ((TEx+R_emis)+0.02), "$R_{emis}$ = %.3f eV"%(R_emis),color="green",fontsize = 14) # stating relaxation energy for emission. 
        
        ax.scatter(Qx,TEx,marker="X",s=200,color="yellow", edgecolors="k")
        
        ax.text(2.5,((self.E_q2+self.E_q1)/2)-0.01,"$E_{21}$ = %.3f eV"%(E_21),color='sienna',fontsize=14)
        ax.text(-2.5,((self.E_q2+self.E_q1)/2)-0.01,"$E_{12}$ = %.3f eV"%(E_12),color='sienna',fontsize=14)
        
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, ncol=4, fontsize =14)  
        ax.set_ylim(lwrE-0.05,limit+0.35)
        plt.ylabel('Total Energy (eV)', fontsize =16) 
        plt.xlabel('Q (Å)', fontsize =16)  
        plt.savefig("./Relaxation_energies/{}.png".format(self.name))
        plt.show()

        self.CTL_2_exp_comparison(R_cap,R_emis)
        
    def energies(self):
        """ E21 = E2 - E1 
            E12 = E1 - E2
            
            Franck–Condon principle energies - Fig. 50  
                ε˚12 = E*2 - E1
                εR12 =  ε˚12 - E21 [Eqn 78]
                     = E*2 - E1 - (E2 - E1) 
                     = E*2 - E2
                     
                ε˚21 = E2 - E*1
                εR21 = E21 - ε˚21 [rearrangement of Eqn 81]
                     = E2 - E1 - (E2 - E*1)
                     = - E1 + E*1
        """
        E_21 = self.E_q2 - self.E_q1
        E_12 = self.E_q1 - self.E_q2
        
        ε_R21 = self.E_q1_rsk - self.E_q1
        ε_R12 = self.E_q2_rsk - self.E_q1 - (self.E_q2 - self.E_q1)
        
        return ε_R12, ε_R21, E_21, E_12
    
    def objective(self, Q, w, qi, Ei):
        """ Eqn 77 
        
            ET_i = 0.5Mω_i^2 (Q-Q_i)^2 + E_i 
            
            ET_i is the total energies of charge state i, M is the associated modal mass, ω_i is the angular frequency of charge state i, Q is reaction coordinates which act as x, Q_i is the reaction coordinate of the charge state minimum, E_i is the minimum total energy of the charge state.  
        """
        return (0.5 * self.M * (w**2) * ((Q-qi)**2)) + Ei
    
    def curvefitting(self, Q, E, E_Di, boolean):
        """ fitting to date via use of parabolic eqn in self.objective. 
        
            inputs:
                E(list):energy points to be plotted for charge state 
                Q(list): reaction coordinate points to be plotted for charge state.
                E_Di(float): minimum total energy of the charge state
                boolean(boolean): whether there are only two entries in E and Q  
                
            outputs:
                popt(list): list of fitted coefficients 
        """
        if boolean == False:
            popt, _ = curve_fit(self.objective,Q,E,bounds=([-np.inf, 0, E_Di-0.00000001],[np.inf, 1, E_Di+0.00000001]))# popt[0]=ω, popt[1]=qi, popt[2]=Ei
        elif boolean == True:
            popt, _ = curve_fit(self.objective,Q,E,bounds=([((4/5)*self.fitting_csts1[0]), 0.88, E_Di-0.00000001],[self.fitting_csts1[0], 1.1, E_Di+0.00000001]))
        ω = popt[0] * np.sqrt(((1.6022 * 10 ** -19) /((1.66054 * (10 ** -27)) * (10 ** -10)**2))) # ω = sqrt(2 ET/M Q_rel^2)
        print(str('TE = 0.5Mω^2(Q-qi)^2 + Ei -> M = %.5f amu; ω = %.5f eV(1/2)•amu-(1/2)•Å-1  = %.4E s-1; qi = %.1f Å; Ei = %.5f eV' % (self.M,popt[0],Decimal(ω),popt[1],popt[2])))
        
        return popt
     
    def relax_E(self): 
        """ Finding the barrier and crossing between the two charge states 
            (1) 0.5Mω_1^2 (Q-q_1)^2 + E_1 = 0.5Mω_2^2 (Q-q_2)^2 + E_2 
            (2)     a      (x-e)^2  +  c  =    b      (x-f)^2   +  d
            
            (3) x= [ea - bf- \sqrt{ abf^2 - 2eabf + ad + e^2ab + cb - ac - bd} ]/{a-b}
            
            outputs:
                Q(float): reaction coordinate of crossing point
                y(float): energy of crossing point 
                R_cap(float): relaxation energy of capture
                R_emis(float): relaxation energy of emission
                col1(str): color of R_cap notation on graph based on where R_cap is from crossing to minimum of q1 curve or q2 curve 
                col2(str): color of R_emis notation on graph based on where R_emis is from crossing to minimum of q1 curve or q2 curve
                lwrE(float): lower charge state minimum energy
                opp(list of float/int): constant to add to x positioning of text for R_cap and R_emis on graph. 
                hgr_rsk(float): higher excited charge state energy
        """
        
        if self.E_q2 < self.E_q1:
            lwrE, opp, hgr_rsk = self.E_q2, [-2,-2], self.E_q1_rsk
            
        else:
            lwrE, opp, hgr_rsk = self.E_q1, [+ 0.02, + 0.06], self.E_q2_rsk
            
        #defining variables corresponding to eqn (2) in docstring
        a, b, c, d, e, f = 0.5*self.M*self.fitting_csts2[0]*self.fitting_csts2[0], 0.5*self.M*self.fitting_csts1[0]*self.fitting_csts1[0], self.fitting_csts2[2], self.fitting_csts1[2], self.fitting_csts2[1], self.fitting_csts1[1]
            
        # defining variables corresponding to eqn (3) in docstring 
        ea, bf, abf2, eabf, ad, e2ab, cb, ac, bd = e*a, b*f, a*b*f*f, e*a*b*f, a*d, e*e*a*b, c*b, a*c, b*d
        Q = (ea - bf - np.sqrt(abf2 - 2*eabf + ad + e2ab + cb - ac - bd))/(a - b)

        y = (0.5 * self.M * (self.fitting_csts1[0]*self.fitting_csts1[0])*((Q - self.fitting_csts1[1])*(Q - self.fitting_csts1[1]))) + self.fitting_csts1[2]
        
        R_cap = self.E_q2 - y
        R_emis = self.E_q1 - y
        
        return Q, y, R_cap, R_emis, lwrE, opp, hgr_rsk
        
    def CTL_2_exp_comparison(self, R_cap, R_emis):
        """ Converting defect CT level into a SRH (ε12 & ε21) levels comparible with experimental extracted trap/interface (e- capture & e- emission) levels collected from TDRC and Id-DLTS.
            Both experimental techniques measure electron traps via either emission (TDRC) or capture (Id-DLTS), so CTL value must be converted from being wrt VBM [self.CTL] to being wrt CBM [CTL_].
            
            ε12 = (εR + E21)^2/4εR [Eqn 89]
                = (R_cap + CTL_)^2/4R_cap 
            ε21 = (εR - E21)^2/4εR [Eqn 90]
                = (R_emis - CTL_)^2/4R_emis
            
            inputs:
                R_cap(float): relaxation energy of capture
                R_emis(float): relaxation energy of emission
        """
        if self.CTL > adiabatic_potentials.bandgap/2: # defect lvl close to CB
            CTL_ = adiabatic_potentials.bandgap - self.CTL # CT lvl wrt to CBM is bandgap - lvl wrt to EBM

            self.ε12 = ((R_cap + CTL_)**2) / (4 * R_cap)
            self.ε21 = ((R_emis - CTL_)**2) / (4 * R_emis)
        elif self.CTL < adiabatic_potentials.bandgap /2: # defect level close to VB
            # must calculate ε12 and ε21 wrt to VBM first.
            ε12 = ((R_cap + self.CTL)**2) / (4 * R_cap)
            ε21 = ((R_emis - self.CTL)**2) / (4 * R_emis)
            # then convert to ε12 and ε21 wrt to CBM via bandgap -  ε12; bandgap - ε21 
            self.ε12 = adiabatic_potentials.bandgap - ε12 
            self.ε21 = adiabatic_potentials.bandgap - ε21 
        
    def returnCTL2exp(self):
        return self.ε12, self.ε21

def indirectory(directory, flname):
    """
        Check if a certain file is already present within a directory. 
        
        inputs:
            directory(str): path of directory within which the file is being searched for 
            flname(str): name of file being searched for
        output: 
            in_(boolean): answer to where file is in directory 
    """
    in_ = None # boolean set up to identify where a certain file is found within a directory 
    
    if len(os.listdir(directory)) == 0: # if directory is empty
        in_ = False
    else:
        for entry, item in zip(os.scandir(directory),range(len(os.listdir(directory)))):
            if entry.is_file() and entry.name.startswith(str(flname)): # file found in directory
                in_ = True 
            if item == len(os.listdir(directory)) - 1 and in_ == None: # if at the end of list of items in the directory and file has not been found.
                in_ = False
    return in_

def roundlistelements(list,pos,num):
    """
        Function to found float elements within a list to a certain number of decimal places. 
        
        Inputs:
            list: the whole list containing the elements to be rounded. 
            pos(int): position of element wanting to round
            num(int): number of decimal places which element is rounded to. 
            
        Output:
            list: the whole list containing the element which has now been rounded. 
    """
    strg = str(list[pos])
    if strg[strg.index('.')+num+1:] == str(5): # for when the digit the rounding is based upon is a 5 (as python inbuilt round function has werid behaviour when rounding based on the digit 5), complete this test. 
        L_pos = str(list[pos]) # str version of list element before being rounded 
        num1 = num + L_pos.index('.') + 1
        posn = round(list[pos],num)
        if float(L_pos[:num1]) == posn: # check if python round function has rounded list element down to the number of decimal places. 
            list[pos] = list[pos] + (1 * (10 ** (-num-1))) # changing the digit of the decimal place before the cut of number of decimals places from 5 to 6. 
            list[pos] = round(list[pos],num) 
            
        else: 
            list[pos] = posn
        
    else:
        list[pos] = round(list[pos],num) 
    return list
