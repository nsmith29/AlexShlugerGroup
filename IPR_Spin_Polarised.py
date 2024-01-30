#!/usr/bin/env python

# Full IPR analysis script  
# Adaption of original script by Jack Strand
#---------------------------------------------------
# Usage: 
# Make sure pdos files obtained from IPR calc are in directory with log file of calc,
# python file get-smearing-pdos.py and pdos.py must also be in directory
#
# run via:
#        ./IPR_Spin_Polarised.py IPR.log 
#
# Output: 
#         IPR_output-Alpha.dat: alpha ipr spectrum
#         IPR_output-Beta.dat: beta ipr spectrum
#         smearedAlpha_kn-X.dat: smeared alpha DOS of atomic kind X from ALPHA_kn-1.pdos file
#         smearedBeta_kn-X.dat: smeared beta DOS of atomic kind X from BETA_kn-1.pdos file 
#         MO_Alpha_projection.dat: atom projection for each alpha MO with an IPR value above 0.01
#         MO_Beta_projection.dat: atom projection for each beta MO with an IPR value above 0.01
#         ipr_pdos.png: figure of system's ipr and pdos data, with title which gives the average 
#                       ipr values of the VB and CB, the ipr value of the defect state, and the
#                       defect state's Fermi level above the single electron KS HOMO energy.
#---------------------------------------------------
# Author: Niamh Smith [orignal: Jack Strand] e-mail: niamh.smith.17 [at] ucl.ac.uk
# Date:   22-09-2023; 30-01-2024 [updated]
#---------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import subprocess


class IPR_Spin_Polarised:
    def __init__(self, input_file):
        self.input_file=input_file
        
        ### calculation project name
        self.project_name = self.find_project_name()
        
        ### Beta Calculations
        self.IPR_Beta=self.Read_data('Beta')
        
        ### Alpha Calculations
        self.IPR_Alpha=self.Read_data('Alpha')
        
    ############################# finding name of project from inputted log file
    def find_project_name(self):
        name_found=False
        logfile_project_name =' GLOBAL| Project name'
        
        with open(self.input_file, 'r') as file:
            lines = file.readlines()
            index = 0
            
            for line in lines:
                index +=1
                if logfile_project_name in line and name_found is False:
                    name_found = True
                    string = line 
                    
        string_arr = string.split()
        project_name = string_arr[3]
        return project_name
    
    ############################# reading MO data from inputted log file
    def Read_data(self, Spin):
        record_alpha=False    
        record=False
        ipr_full_info=[] # MO,Energy,atom number, element type ,a^2, a^4
        current_mos=[]
        temp_array=[]
        temp_array_2=[]
        basis_number=0
        last_array=[]
        first_MO_number=1
        index = len(open(self.input_file).readlines()) + 1
        for line in reversed(open(self.input_file).readlines()):
            index-=1
            if record == True:
                data_line=line.split()
                
                if basis_number == 1:       #For when we get to the final basis number (labelled '1', since we are going backwards)
                    # print(index,"if basis_number == 1; script read final basis number correctly")
                    # print(index, "if basis_number == 1; data_line", data_line)
                    if len(data_line) != 1:
                        last_array.append(data_line)
                        print(index, "if basis_number == 1; if len(data_line) != 1; last_array", last_array)
                    if len(last_array)==3:
                        last_array[2].remove('MO|')
                        last_array[1].remove('MO|')
                        print(index,"if basis_number == 1; if len(data_line) != 3; last_array",'[2]:',last_array[2],'[1]:',last_array[1])
                        
                        for i in range(len(last_array[2])):
                            q=len(last_array[2])-i-1
                            ipr_full_info.append([int(last_array[2][q]),last_array[1][q],last_array[0][q],temp_array_2[q]]) 
                            
                        print(last_array[2][0])
                        
                        if int(last_array[2][0])==first_MO_number:
                            break
                        temp_array=[]
                        temp_array_2=[]
                        last_array=[]  
                        basis_number=0
                        current_mos=[]
                        
                elif len(data_line) != 1:
                    # print(data_line)
                    basis_number=int(data_line[1])
                    if basis_number == 1:
                        print(index, 'len(data_line) != 1; basis_number == 1','data_line:',data_line)
                    atom_number=int(data_line[2])
                    atom_type=data_line[3]
                    num_of_mos=len(data_line) - 5
                    
                    if len(temp_array) == 0:
                        temp_array.append(data_line[1])
                    #Calculate how many MOs we are currently dealing with (will usually be 4, someimes the final batch will be less than 4)
                        for i in range(num_of_mos):
                            temp_array.append([0,0])
                            temp_array_2.append([]) 
                            
                    #Calculate and store the sum of a^2 and a^4 for each MO in the current batch
                    for i in range(num_of_mos):
                        temp_array[1+i][0]=temp_array[1+i][0]+float(data_line[5+i])**2
                        temp_array[1+i][1]=temp_array[1+i][1]+float(data_line[5+i])**4   
                        
                   #Store in temp array with format [Atom number, Atom type,  sum of a^2, sum of a^4]
                if len(data_line) == 1 and len(temp_array) != 0:
                    
                    for i in range(num_of_mos):
                        temp_array_2[i].append([atom_number,atom_type,temp_array[1+i][0],temp_array[1+i][1]])   
                    
                    temp_array=[]     
                    
            if len(line.split()) == 6 and line.split()[1] == 'E(Fermi):' and Spin=='Beta':
                record=True 
                
            if len(line.split()) == 6 and line.split()[1] == 'E(Fermi):' and record_alpha==True :
                record=True
                
            if len(line.split()) == 6 and line.split()[1] == 'E(Fermi):' and Spin=='Alpha':
                record_alpha=True
                
        return ipr_full_info 

    ############################# Calculate IPR Spectrum
    def IPR_Spectrum(self,Spin):
        if Spin == 'Alpha':
            s=1
            ipr_full_info = self.IPR_Alpha
            
        else:
            s=-1
            ipr_full_info = self.IPR_Beta
            
        a_2=0
        a_4=0
        ipr=0
        MO_ipr_values=[]
        
        for MO in reversed(ipr_full_info):
            
            for atom in MO[3]:
                a_2=a_2+atom[2]
                a_4=a_4+atom[3]
                
            ipr=s*a_4/(a_2**2)
            a_2=0
            a_4=0
            
            MO_ipr_values.append([MO[0],MO[1],ipr])
            
        output_file=open("IPR_output-"+Spin+".dat","w")
        
        for MO in MO_ipr_values:
            output_file.write(str(MO[1])+"   "+str(MO[2])+"\n")
            
        output_file.close()
        
    ########################### defining average ipr value of the VB and CB
    def bandIPR_ave(self, homo_MO, lumo_MO, IPR):
        # VB_ave, CB_ave = None, None
        for band, start, finish in zip(['VB','CB'],[1,int(lumo_MO)],[int(homo_MO),int(len(IPR))]):
            adding = 0
            counter = 0
            for n in range(start,finish):
                counter += 1
                adding = adding + IPR[n]
            _ave = round(adding / counter,6)
            if band == 'VB':
                VB_ave = _ave
            else:
                CB_ave = _ave
        return VB_ave, CB_ave
        
    ############################# processing pdos files to get pdos.dat files, adjusting IPR eigenvalue for fermi level, finding localised MOs in IPR analysis. 
    def pdos_ipr_processing(self, cwd):
        for S, s in zip(['ALPHA','BETA'],['Alpha','Beta']):
            MOs_written = False
            # unpacking data from IPR output files 
            IPR_file = str("IPR_output-{}.dat".format(s))
            exec(f'eigenvalue0_{S}, IPR_{s} = np.loadtxt(IPR_file, unpack=True)')

            if S == 'ALPHA':
                IPR_ALPHA = eval("IPR_{}".format(s))
                Alpha_dats, MOs_Alpha, is_Alpha_occ, is_Alpha_unocc = [], [], [], []
                
            elif S == 'BETA':
                IPR_BETA = eval("IPR_{}".format(s))
                Beta_dats, MOs_Beta, is_Beta_occ, is_Beta_unocc = [], [], [], []
                
            eigenvalue = eval("eigenvalue0_{}".format(S))
            # file name basis of pdos files
            pdos_stem = str("{}-{}".format(self.project_name,S))
            
            for entry in os.scandir(cwd):
                
                ## finding pdos files within current working directory
                if entry.is_file() and entry.name.startswith(pdos_stem):
                    pdos_file = entry.name
                    
                    # getting atomic kind and fermi level from pdos 
                    with open(pdos_file, 'r') as file:
                        firstline = file.readline().strip().split()
                        secondline = file.readline().strip().split()
                        kind = firstline[6]
                        efermi = float(firstline[15])
                        var_num = len(secondline) # number of headings in the pdos file
                        
                    eigen_energy = [(float(ex) - efermi)*27.211384523 for ex in eigenvalue]
                    
                    if S == 'ALPHA':
                        eigenvalue_ALPHA = eigen_energy
                        
                    elif S == 'BETA':
                        eigenvalue_BETA = eigen_energy
                        
                    # getting name to change 'smeared.dat' to for each pdos
                    pdos_list = []
                    
                    for l in range(len(self.project_name),len(pdos_file)-6): 
                        L = pdos_file[l] 
                        pdos_list.append(L) 
                        
                    name = str('smeared' + "".join(pdos_list) + kind +'.dat') 
                    exec(f'{s}_dats.append(name)')
                    # executing get-smearing-pdos.py and renaming 'smeared.dat'
                    subprocess.run(["python3", "get-smearing-pdos.py",str("{}".format(pdos_file))])
                    subprocess.run(["mv","smeared.dat",str("{}".format(name))])
                    
                    # unpacking the pdos files
                    if var_num == 7:
                        Mo, Eigenvalue, Occupation, S_occ, P_occ = np.loadtxt(pdos_file, unpack=True)
                        
                    elif var_num == 8:
                        Mo, Eigenvalue, Occupation, S_occ, P_occ, D_occ = np.loadtxt(pdos_file, unpack=True)
                        
                    elif var_num == 9:
                        Mo, Eigenvalue, Occupation, S_occ, P_occ, D_occ, F_occ = np.loadtxt(pdos_file, unpack=True)

                    # finding the MOs of the homo and lumo states
                    for occ, val, band, indx in zip(['occ','unocc'], [1.0,0.0], ['homo','lumo'], [-1,0]):
                        for j in range(len(Occupation)): 
                            if round(Occupation[j],0) == val: 
                                exec(f'is_{s}_{occ}.append(j)')
                        exec(f'{band}_{S}_MO = int(Mo[is_{s}_{occ}[indx]])')

                    # calculate average IPR values of the valence and conduction bands
                    VB_ave, CB_ave = self.bandIPR_ave(eval("homo_{}_MO".format(S)), eval("lumo_{}_MO".format(S)), eval("IPR_{}".format(S)))

                    if MOs_written is False:
                        for n, m in zip(eval("IPR_{}".format(S)),eigenvalue):
                            
                            if abs(float(n)) >= 0.01:
                                exec(f'MOs_{s}.append(Mo[np.where(Eigenvalue == m)[0][0]])')
                                
                        MOs_written = True

                    defectstatefound = False
                    for i,j in zip(range(len(eval("eigenvalue_{}".format(S)))),range(len(eval("IPR_{}".format(S))))):
                        # findng the eigenvalue and ipr values of the defect state
                        if eval("eigenvalue_{}[i]".format(S)) > 0 and abs(eval("IPR_{}[j]".format(S))) > abs(VB_ave) and abs(eval("IPR_{}[j]".format(S))) > abs(CB_ave) and defectstatefound is False:
                            def_state_eigen, def_state_IPR = eval("eigenvalue_{}[i]".format(S)), round(eval("IPR_{}[j]".format(S)),6)
                            VB_Edge  = Eigenvalue[np.where(Mo == eval("homo_{}_MO".format(S)))[0][0]]
                            # calculation of defect state fermi level of homo 
                            DS_Ef = round(((efermi - VB_Edge)*27.211) + def_state_eigen,6)
                            
                            defectstatefound = None
                        
        Alpha_dats = np.sort(Alpha_dats) 
        Beta_dats = np.sort(Beta_dats) 
                    
        return eigenvalue_ALPHA, eigenvalue_BETA, IPR_ALPHA, IPR_BETA, Alpha_dats, Beta_dats, MOs_Alpha, MOs_Beta, VB_ave, CB_ave, DS_Ef, def_state_IPR
    
    ############################# Project MOs onto atoms
    def MO_projection_atoms(self, MO_number,number_of_atoms,spin):
        if spin == 'alpha':
            ipr_full_info = self.IPR_Alpha
            
        else:
            ipr_full_info = self.IPR_Beta

        atoms=[]
        a_2_values=[]
        
        for i in range(number_of_atoms):
            atoms.append('N')
            a_2_values.append(0)
            
        for entry in ipr_full_info:
            
            if entry[0] == MO_number:
                
                for atom in entry[3]:
                    
                    for i in range(len(a_2_values)):
                        
                        if atom[2] > a_2_values[i]:
                            atoms[i]=atom[0]
                            a_2_values[i]=atom[2]
                            break
                            
        output_list=[]
        
        for i in range(len(a_2_values)):
            output_list.append([atoms[i],a_2_values[i]])
            
        return output_list

    ########################### Function to project several MOs
    def MO_projection_scan(self, list_of_MOs, spin):
        output_list=[]
        
        for element in list_of_MOs:
            output_list.append([element,self.MO_projection_atoms(element,3,spin)])
            
        filename = str("MO_{}_projection.dat".format(spin))
        output_file=open(filename, "w")

        for element in output_list:
            output_file.write(str(element[0]))
            output_file.write("\n")
            
            for atom in element[1]:
                output_file.write(str(atom[0])+"   "+str(atom[1]))
                output_file.write("\n")
                
            output_file.write("\n")    
            
        output_file.close()
        
    ########################### deriving colours to be used within plot
    def get_cmap(self, n, name='hsv'):
    #Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    # RGB color; the keyword argument name must be a standard mpl colormap name.
        return plt.cm.get_cmap(name, n)
    
    ########################### function to create figure showing pdos and ipr 
    def plotting_pdos_IPR(self, eigenvalue_ALPHA, eigenvalue_BETA, IPR_ALPHA, IPR_BETA, Alpha_dats, Beta_dats, VB_ave, CB_ave, DS_Ef, def_state_IPR):
        fig = plt.figure(figsize=(13,8))
        ax = plt.subplot()
        ax.set_xlabel('Energy - Fermi level energy (eV)',fontsize =22)
        ## first y-axis for pdos data
        ax.set_ylabel('Density of states',fontsize =22)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        ## second y-axis for ipr data 
        ax2 = ax.twinx() 
        ax2.set_ylabel('IPR',fontsize =22)
        ax.plot([-6,6],[0,0],ls = '--', color = 'grey', lw=1)
        # to show the different spin directions on the graph
        ax.text(3,-175,'BETA',fontsize=20,alpha = 0.75)
        ax.text(3,175,'ALPHA',fontsize=20,alpha = 0.75)
        cmap = self.get_cmap(len(Alpha_dats)+1)
        
        ## plotting alpha energies and densities of each pdos kind
        for i, dat in zip(range(len(Alpha_dats)),list(Alpha_dats)):
            dat_stem = 'smeared-Alpha'
            energy, density = np.loadtxt(dat, unpack=True)
            dat_letters = []
            
            # getting label to put in legend
            for l in range(len(dat_stem)+4,len(dat)-4):
                L = dat[l]
                dat_letters.append(L)
                
            ax.plot(energy,density,ls='-',c=cmap(i),label="".join(dat_letters))
            
        ## plotting beta energies and densities of each pdos kind
        for i, dat in zip(range(len(Beta_dats)),list(Beta_dats)):
            dat_stem = 'smearedBeta'
            energy, density = np.loadtxt(dat, unpack=True)
            density = - density 
            ax.plot(energy,density,ls='-',c=cmap(i))

        
        ## plotting the ipr energies and ratios
        for S in 'ALPHA','BETA':
            
            for i,j in zip(range(len(eval("eigenvalue_{}".format(S)))),range(len(eval("IPR_{}".format(S))))):
                eigen_energy = eval("eigenvalue_{}[i]".format(S))
                IPR =  eval("IPR_{}[j]".format(S))
                ax2.plot([eigen_energy,eigen_energy],[0,IPR],ls ='-',color='k',lw=1.75)

        eqn1, eqn2, eqn3, eqn4 = str("$\overline{IPR(\psi_{VB})}$"), str("$\overline{IPR(\psi_{CB})}$"), str("$\delta\epsilon_{f}[\psi_{DS}]$"), str("$\overline{IPR(\psi_{DS})}$")
        string = str("PDOS & IPR: {} = {}, {} = {}, \n {} = {}, {} = {} ".format(eqn1,VB_ave,eqn2,CB_ave,eqn3,DS_Ef,eqn4,def_state_IPR))
        
        ax.set_title(string,fontsize =22)
        
        hands, labs = ax.get_legend_handles_labels()
        ax.legend(hands, labs, bbox_to_anchor=(0.25, 0.7), ncol=2, fontsize=22) 
        
        ax.set_ylim(-200,200) # y-axis limits for pdos y-axis
        ax2.set_ylim(-0.09,0.09) # y-axis limits for ipr y-axis 
        ax.set_xlim(-5,5) # x-axis limit for energy range wanting to be shown
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.savefig('ipr_pdos.png')
        plt.show()

input_file=sys.argv[1]
ipr = IPR_Spin_Polarised(input_file)
# producing the IPR output dat files for both spins
ipr.IPR_Spectrum('Beta')
ipr.IPR_Spectrum('Alpha')

## pdos file processing
cwd = os.getcwd()
eigenvalue_ALPHA, eigenvalue_BETA, IPR_ALPHA, IPR_BETA, Alpha_dats, Beta_dats, MOs_Alpha, MOs_Beta, VB_ave, CB_ave, DS_Ef, def_state_IPR =ipr.pdos_ipr_processing(cwd)
for s in 'Alpha','Beta':
    ipr.MO_projection_scan(eval("MOs_{}".format(s)),s)
ipr.plotting_pdos_IPR(eigenvalue_ALPHA, eigenvalue_BETA, IPR_ALPHA, IPR_BETA, Alpha_dats, Beta_dats, VB_ave, CB_ave, DS_Ef, def_state_IPR)           
            
