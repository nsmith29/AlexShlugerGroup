""" _Python_automation_of_CP2K_calculation_results_collection_ 

    READ BEFORE USE
    
    Currently code allows automatic extraction of:
    - Total Energies, Band gap width, Mulliken and Hirshfeld charges, spins, and populations from .log files
    - HOMO and LUMO eigenvalues converted into eV, smeared DOS .dat files, and plot-ready energy and density arrays from .pdos files
    - atomic co-ordinates obtained in last geometry optimisation step of calculation, plot-ready displacements of atoms in last geometry step from perfect defect free geometry
      and distances of atoms after geometry optimisation from defect site from .xyz files

    Code originally written for being run in a jupyter notebook/lab environment. 
    Indication of where code should be inputted into a new jupyter code cell is given by:
# ====================================================================================
    
    Code variables which need to be individually set to reflect your specific calculations are specified as {...} followed by a comment '#' which identifies what the variable is 
    
    General comments about what the code is doing will be specified as '##'. 
    
    Author: Niamh Smith [niamh.smith.17@ucl.ac.uk]
    v1.0 (02/06/2022)
    
    Please get in touch if having trouble adapting or adding to the code 
""" 
# ====================================================================================
import numpy as np
from math import pi, sqrt	
import numpy as np
import sys

# ====================================================================================
class pdos: ## projected electronic density of states from CP2K output files 
    def __init__(self, infilename): ## read a CP2k .pdos file and build a pdos instance
        input_file = open(infilename, 'r')
        firstline = input_file.readline().strip().split()
        secondline = input_file.readline().strip().split()
        
        self.atom = firstline[6] ## atom name of kind where DoS is projected
        self.iterstep = int(firstline[12][:-1]) ## iteration step from CP2K job, [:-1] deletes ","
        self.efermi = float(firstline[15]) ## energy of fermi level [a.u.]
        
        secondline[0:5] = []
        self.orbitals = secondline ## keeps just the orbital names 
        
        lines = input_file.readlines()
        
        eigenvalue = []
        self.occupation = []
        data = []
        self.pdos = []
        for index, line in enumerate(lines):
            data.append(line.split())
            data[index].pop(0)
            eigenvalue.append(float(data[index].pop(0)))
            self.occupation.append(int(float(data[index].pop(0))))
            self.pdos.append([float(i) for i in data[index]])
            
        self.e = [(x-self.efermi)*27.211384523 for x in eigenvalue]
        
        self.tpdos = []
        for i in self.pdos:
            self.tpdos.append(sum(i))
            
    def __add__(self, other): ## Return sum of two PDOS objects
        sumtpdos = [i+j for i,j in zip(self.tpdos,other.tpdos)]
        return sumtpdos 
    
    def delta(self,emin,emax,npts,energy,width): ## return a delta-function centered at energy, emin= min eigenvalue, emax= max eigenvalue, npts= # of points in smeared pdos, energy= energy where gaussian 
                                                 ## is centered, width= dispersion parameter, delta: array of delta function values 
        energies = np.linspace(emin, emax, npts)
        x = -((energies - energy) / width)**2
        return np.exp(x) / (sqrt(pi) * width) 
    
    def smearing(self, npts, width,): ## return a gaussian smeared DOS
        d = np.zeros(npts)
        emin = min(self.e)
        emax = max(self.e)
        for e, pd in zip(self.e,self.tpdos):
            d +=pd*self.delta(emin,emax,npts,e,width)
            
        return d
        
# ====================================================================================
total_atoms_in_system = {...} # number of atoms in material included within calculations

tot_atoms = str("     {}".format(total_atoms_in_system))
total_atoms_plus_2 = total_atoms_in_system + 2
total_atoms_in_calc = []
number_of_total_atoms_array = np.linspace(1, total_atoms_plus_2, total_atoms_plus_2) - 1
for t in number_of_total_atoms_array:
    t = round(t)
    total_atoms_in_calc.append(t)

gap = 'HOMO - LUMO gap'
energy = 'ENERGY| Total FORCE_EVAL'
pop1 = 'Mulliken Population Analysis'      
pop2 = 'Hirshfeld Charges'
# ====================================================================================

number_of_kinds_in_perfect = {...} # number of different atomic kinds included within perfect calculation
kinds_in_perfect = np.linspace(1, number_of_kinds_in_perfect, number_of_kinds_in_perfect) - 1
kinds_in_perfect = np.around(kinds_in_perfect, decimals= 1)
perfect_included_atoms = '{...}' # symbol of elements included in perfect calculation listed in same order as their kind sections in input file 

number_of_kinds_in_defect = {...} # number of different atomic kinds included within each defect calculation 
kinds_in_defect = np.linspace(1, number_of_kinds_in_defect, number_of_kinds_in_defect) - 1
kinds_in_defect = np.around(kinds_in_defect,0)
defect_included_atoms = '{...}' # symbol of elements included in defect calculations listed in same order as their kind sections in input file; if different defects are related to addition of
                              # different impurities comment out line and uncomment first lines after "else" of "if X == perfect" statement.

# ====================================================================================
perfect_project_name = str("{...}") # PROJECT_NAME set within no defect pefect material calculation cp2k input file 
perfect_logfile_name = str("{...}") # name of .log file set within ARCHER2 submission script of no defect calculation, do not include '.log' extension in variable

defect_project_name_stem = str("{...}") # PROJECT_NAME stem set within defect calculation cp2k input files
defect_logfile_name_stem = str("{...}") # name of .log file set within ARCHER2 submission scripts of defect calculations, do not include '.log' extension in variable

output_file_hierarchy_directory = str("{...}") # path to directory where output files of all calculations are store in subdirectories on local machine/server. str needs to end with /

# ====================================================================================
for X in 'perfect','{...}': # list all defects in material being studied. If defects are related to addition of impurity atoms, list the symbols of impurities used in the different calculations rather than defect name
    if X == 'perfect':
        system = str("perf_" + perfect_project_name)
    else:
        # defect_included_atoms = ['{...}',X] # comment out line if not dealing with impurity related defects. 
        system = str(X)
        
    for state in '{...}':# charge state of defects; for only one charge state comment out line and write replacement line state = {...} below. Remember to move code lines included within for state in {...} 
                        # back an indent which can be done by highlighting code and pressing shift tab. 
        if X == 'perfect':
            state = '{...}' # state may need to be reset for perfect structure to stop any iteration errors from occuring for charge states incompatible with perfect structure 
            file = str("{...}") # remaining path to subdirectory containing only perfect calculation results. Should continue on hierachy directory path saved in output_file_hierarchy_directory. Str needs to end with /
            logname = str("{}.log".format(perfect_logfile_name))
            for k in kinds_in_perfect:
                kr = round(k)
                K = kr + 1
                pdosname_alpha = str("{}-ALPHA_k{}-1.pdos".format(perfect_project_name,K))
                exec(f'pdosname_alpha{kr} = pdosname_alpha')
                pdosname_beta = str("{}-BETA_k{}-1.pdos".format(perfect_project_name,K))
                exec(f'pdosname_beta{kr} = pdosname_beta')
        else:
            file = str("{...}") # remaining path to subdirectory comtaining calculation results for a given defect in X [other than perfect]. This section of path needs to be iteratible with X and state. Str needs to end with /
            logname = str("{}_{}.log".format(X,defect_logfile_name_stem)) ## if logfile names for each individual defect calculation is the same, change to logname = str('{}.log'.format(defect_logfile_name_stem)), line at current assumes different 
                                                                          ## logfile name prefix for each individual defect with defect_logfile_name_stem being the suffix
            for k in kinds_in_defect:
                kr = round(k)
                K = kr + 1
                pdosname_alpha = str("{}{}-ALPHA_k{}-1.pdos".format(defect_project_name_stem,X,K))
                exec(f'pdosname_alpha{kr} = pdosname_alpha')
                pdosname_beta = str("{}{}-BETA_k{}-1.pdos".format(defect_project_name_stem,X,K))
                exec(f'pdosname_beta{kr} = pdosname_beta')
        exec(f'{system}_{state}_energy_line = []')
        exec(f'{system}_{state}_gap_line = []')
        exec(f'{system}_{state}_gap= []')
        exec(f'{system}_{state}_energy = []')
        for n in 1, 2:
            for A in '{...}': # list of atom index and chem formula of important atoms (i.e. substituted atom [if impurity defect, not vacancy defect] nearest atomic neighbours) around defects - e.g. C304, C313
                for var in 'charge', 'spin', 'beta_pop', 'alpha_pop':
                    exec(f'{system}_{state}_Pop{n}Impurity{A}_line = []')
                    exec(f'{system}_{A}_{state}_pop{n}_{var} = []')
                    
        logfile = str("{}{}{}".format(output_file_hierarchy_directory,file,logname)) 
        log  = open(logfile,'r') ## open file name saved as logfile
        index = 0
        for line in log: ## iterate over every line in file
            index += 1
            if energy in line:  ## if the string called energy is found within a line
                exec(f'{system}_{state}_energy_index = index - 1') ## first line is treated as line0 so number of line found containing string must be adjusted
                exec(f'{system}_{state}_energy_line.append({system}_{state}_energy_index)') ## the line number of every line in which string is found is appended into an array
            if gap in line: 
                exec(f'{system}_{state}_gap_index = index - 1') 
                exec(f'{system}_{state}_gap_line.append({system}_{state}_gap_index)')
            if pop1 in line: ## if string called pop1 is found within a line
                PopAnalysisStart = index ## record index of line containing pop1 string
                for A in '{...}': #  list of atom index and chem formula of important atoms (i.e. substituted atom [if impurity defect, not vacancy defect] nearest atomic neighbours) around defects - e.g. C304, C313
                    i_{...}_pop1 = PopAnalysisStart + {...} # first {...} needs to = first list item in A ,for second {...}, remove chem. formula and add 1 to index number of atom [first line is treated as line0]
                                                          ## get line number of Mulliken Population analysis results for atom {...} in list A
                    {...}etc # repeat line above for each subsequent list items listed in A
                    exec(f'{system}_{state}_Pop1Impurity{A}_line.append(i_{A}_pop1)') ## specific array created for each atom in list A appending every Mulliken Population analysis line number for the atom
            if pop2 in line: ## record index of line containing pop2 string
                PopAnalysisStart = index 
                for A in '{...}': #  list of atom index and chem formula of important atoms (i.e. substituted atom [if impurity defect, not vacancy defect] nearest atomic neighbours) around defects - e.g. C304, C313
                    i_{...}_pop2 = PopAnalysisStart + {...} # first {...} needs to = first list item in A ,for second {...}, remove chem. formula and add 1 to index number of atom [first line is treated as line0]
                                                          ## get line number of Hirshfeld Population analysis results for atom {...} in list A
                    {...}etc # repeat line above for each subsequent list items listed in A
                    exec(f'{system}_{state}_Pop2Impurity{A}_line.append(i_{A}_pop2)') ## specific array created for each atom in list A appending every Hirshfield Population analysis line number for the atom
        log.close() ## file must be closed after this
        
        for item in 'energy', 'gap':
            read_lines = eval("{}_{}_{}_line".format(system,state,item)) ## read_lines needs to be define for each item so that the following for loop can be carried out automatically
            log = open(logfile,'r')  ## log file must be reasigned 
            for position, line in enumerate(log): ## enumerate used to search file and return the text string found on line numbers specified by position
                if position in read_lines:
                    strg = line
                    exec(f'arr_{system}_{state}_{item} = strg.split()') ## split text string of line into individual string based on column grouping of line
                    exec(f'arr_{system}_{state}_{item} = float(arr_{system}_{state}_{item}[-1])') ## record line's last column string value as a float
                    exec(f'{system}_{state}_{item}.append(arr_{system}_{state}_{item})') ## append each line's last column float value within an array
            log.close()
            
        exec(f'E_{system}_{state} = {system}_{state}_energy[-1] * 27.211') ## To get final energy of calculation and convert it from hartree units into eV
        exec(f'alpha_HOMO_LUMOgap_{system}_{state} = {system}_{state}_gap[-2]') ## Second to last array entry is the seperation in energy (eV) between HOMO and LUMO of alpha spin state
        exec(f'beta_HOMO_LUMOgap_{system}_{state} = {system}_{state}_gap[-1]') ## Last array entry is the seperation in energy (eV) between HOMO and LUMO of beta spin state
        
        for n in 1, 2: ## for two two different population analyses
            item = str("Pop{}Impurity".format(n)) 
            for A in 'Si246', '250', 'Si373', 'Si357', 'Si242': #  list of atom index and chem formula of important atoms (i.e. substituted atom [if impurity defect, not vacancy defect] nearest atomic neighbours) around defects - e.g. C304, C313
                lines_to_read = eval("{}_{}_{}{}_line".format(system,state,item,A)) ## read_lines needs to be define for each population analysis and atom in list A can be carried out automatically
                log = open(logfile,'r')
                for position, line in enumerate(log):
                    if position in lines_to_read:
                        string = line
                        for s in 'beta', 'alpha': ## for spin polarised calculations where beta and alpha spin are treated seperately
                            exec(f'{system}_{state}_{n}_{A}_arr = string.split()')
                            fullstring = str("{}_{}_{}_{}_arr".format(system,state,n,A)) ## save name of above array in a string called fullstring
                            u1 = str("{}_1_{}".format(state,A)) ## string section of array name indicating array is for Mulliken population analysis
                            u2 = str("{}_2_{}".format(state,A)) ## string section of array name indicating array is for Hirshfield population analysis
                            fullstring2 = str("{}_{}_{}_pop{}_{}_pop".format(system,A,state,n,s)) ## save name of current iteration's population variable as a string called fullstring2
                            if fullstring.find(u1) != -1: ## if array contains string section for Mulliken population analysis
                                exec(f'char1 = float({system}_{state}_{n}_{A}_arr[5])') ## Within Mulliken population table charges are given in column 6 (first column is counted as [0])
                                exec(f'{system}_{A}_{state}_pop{n}_charge.append(char1)')
                                exec(f'spin1 = float({system}_{state}_{n}_{A}_arr[6])') ## Within Mulliken population table spins are given in column 7 (first column is counted as [0])
                                exec(f'{system}_{A}_{state}_pop{n}_spin.append(spin1)')
                                if fullstring2.find("beta") != -1: ## if the current iteration's Mulliken population is for the beta spin state
                                    exec(f'b1 = float({system}_{state}_{n}_{A}_arr[4])') ## Within Mulliken population table beta populations are given in column 5 (first column is counted as [0])
                                    exec(f'{system}_{A}_{state}_pop{n}_{s}_pop.append(b1)')
                                elif fullstring2.find("alpha")  != -1: ## if the current iteration's Mulliken population is for the alpha spin state
                                    exec(f'a1 = float({system}_{state}_{n}_{A}_arr[3])') ## Within Mulliken population table alpha populations are given in column 4 (first column is counted as [0])
                                    exec(f'{system}_{A}_{state}_pop{n}_{s}_pop.append(a1)')
                            elif fullstring.find(u2) != -1: ## if array contains string section for Hirshfield population analysis
                                exec(f'char = float({system}_{state}_{n}_{A}_arr[7])') ## Within Hirshfield population table charges are given in column 8 (first column is counted as [0])
                                exec(f'{system}_{A}_{state}_pop{n}_charge.append(char)')
                                exec(f'spin = float({system}_{state}_{n}_{A}_arr[6])') ## Within Hirshfield population table spins are given in column 7 (first column is counted as [0])
                                exec(f'{system}_{A}_{state}_pop{n}_spin.append(spin)')
                                if fullstring2.find("beta") != -1: ## if the current iteration's Hirshfield population is for the beta spin state
                                    exec(f'b2 = float({system}_{state}_{n}_{A}_arr[5])') ## Within Hirshfield population table beta populations are given in column 6 (first column is counted as [0])
                                    exec(f'{system}_{A}_{state}_pop{n}_{s}_pop.append(b2)')
                                elif fullstring2.find("alpha") != -1: ## if the current iteration's Hirshfield population is for the alpha spin state
                                    exec(f'a2 = float({system}_{state}_{n}_{A}_arr[4])') ## Within Hirshfild population table alpha populations are given in column 5 (first column is counted as [0])
                                    exec(f'{system}_{A}_{state}_pop{n}_{s}_pop.append(a2)')
                log.close()
        
        for s in 'alpha', 'beta': ## for spin polarised calculations where beta and alpha spin are treated seperately
            if X == 'perfect':
                for k in kinds_in_perfect:
                    kr = round(k)
                    pdosname = eval("pdosname_{}{}".format(s,kr))
                    infilename = str("{}{}{}".format(output_file_hierarchy_directory,file,pdosname)) 
                    Pdos = pdos(infilename)
                    npts = len(Pdos.e)
                    pdos_smeared = Pdos.smearing(npts,0.1) 
                    eigenvalues = np.linspace(min(Pdos.e), max(Pdos.e),npts)
                    ## above lines in for loop of k taken from get_smearing_pdos.py
                    kind = perfect_included_atoms[kr] ## corresponding symbol of element in perfect represented by k
                    pdos_dat_file = str("{}{}{}_{}.dat".format(output_file_hierarchy_directory,file,kind,s)) ## name and file path of .dat file being wrote
                    g = open(pdos_dat_file,'w')
                    for i,j in zip(eigenvalues, pdos_smeared):
                        t = str(i).ljust(15) + '     ' + str(j).ljust(15) + '\n'
                        g.write(t) 
                        
                    exec(f'{system}_{state}_{kind}_{s}_energy, {system}_{state}_{kind}_{s}_density = np.loadtxt(pdos_dat_file, unpack=True)') 
                    if s == 'beta':
                        exec(f'{system}_{state}_{kind}_{s}_density = - {system}_{state}_{kind}_{s}_density')
            else:
                for k in kinds_in_defect:
                    kr = round(k)
                    pdosname = eval("pdosname_{}{}".format(s,kr))
                    infilename = str("{}{}{}".format(output_file_hierarchy_directory,file,pdosname)) 
                    Pdos = pdos(infilename)
                    npts = len(Pdos.e)
                    pdos_smeared = Pdos.smearing(npts,0.1) 
                    eigenvalues = np.linspace(min(Pdos.e), max(Pdos.e),npts)
                    kind = defect_included_atoms[kr]
                    
                    pdos_dat_file = str("{}{}{}_{}.dat".format(output_file_hierarchy_directory,file,kind,s))
                    g = open(pdos_dat_file, 'w')
                    for i,j in zip(eigenvalues, pdos_smeared):
                        t = str(i).ljust(15) + '     ' + str(j).ljust(15) + '\n'
                        g.write(t) 
                        
                    exec(f'{system}_{state}_{kind}_{s}_energy, {system}_{state}_{kind}_{s}_density = np.loadtxt(pdos_dat_file, unpack=True)')
                    if s == 'beta':
                        exec(f'{system}_{state}_{kind}_{s}_density = - {system}_{state}_{kind}_{s}_density')
                    
        for s in 'alpha', 'beta':
            if X == 'perfect':
                pdosname = eval("pdosname_{}0".format(s)) ## for perfect structure take both pdos files of first atomic kind
                pdosfile = str("{}{}{}".format(output_file_hierarchy_directory,file,pdosname)) 
            else:
                pdos_kind = kinds_in_defect[-1] ## number of last atomic kind in defect calculations
                pdos_kind = round(pdos_kind)
                pdosname = eval("pdosname_{}{}".format(s,pdos_kind)) ## for defect structures take both pdos files of last atomic kind
                pdosfile = str("{}{}{}".format(output_file_hierarchy_directory,file,pdosname)) 
            exec(f'{system}_{state}_{s}_Mo, {system}_{state}_{s}_Eigenvalue, {system}_{state}_{s}_Occupation = np.loadtxt(pdosfile, usecols=(0,1,2), unpack=True)') ## extracting only first 3 columns of pdos file
            exec(f'{system}_{state}_{s}_Eigenvalue = {system}_{state}_{s}_Eigenvalue * 27.211') ## converting eigenvalues into eV
            exec(f'indices_{s}_occ = []')
            exec(f'indices_{s}_unocc = []')
            Occupation = eval("{}_{}_{}_Occupation".format(system,state,s)) ## Occupation variable made iteratible for for loop regardless of X or State
            for j in range(len(Occupation)): 
                if Occupation[j] == 1.00000: ## if value of entry j of occupation equals 1
                    exec(f'indices_{s}_occ.append(j)')
            for i in range(len(Occupation)):
                if Occupation[i] ==0.00000: ## if value of entry j of occupation equals 0
                    exec(f'indices_{s}_unocc.append(i)')
            exec(f'{system}_{state}_occupied_{s} = {system}_{state}_{s}_Eigenvalue[indices_{s}_occ]') ## all occupied oribtal eigenvalues(eV)
            exec(f'{system}_{state}_HOMO_{s} = {system}_{state}_occupied_{s}[-1]') ## last occupied eigenvalue is eigenvalue of HOMO"
            exec(f'{system}_{state}_HOMO_{s} = round({system}_{state}_HOMO_{s},2)')
            exec(f'{system}_{state}_unoccupied_{s} = {system}_{state}_{s}_Eigenvalue[indices_{s}_unocc]') ## all unoccupied oribtal eigenvalues(eV)
            exec(f'{system}_{state}_LUMO_{s} = {system}_{state}_unoccupied_{s}[0]')
            exec(f'{system}_{state}_LUMO_{s} = round({system}_{state}_LUMO_{s},2)') ## first unoccupied eigenvalue is eigenvalue of LUMO
            
        """IMPORTANT!! On first run of script/jupyter notebook cells comment out all lines which have the comment #!*! at then end of them to generate the -L.xyz files of your calculations. The script/jupyter
        notebook cells then should be rerun again with those lines which have the comment #!*! uncommented and line which comment #¡***¡ commented out. This is so the full -L.xyz files of your calculations 
        can be unpacked to {system}_{state}_x, {system}_{state}_y, {system}_{state}_z"""
        itr_start = []
        L_itr_lines = []
        if X == 'perfect':
            xyz_file = str("{}{}{}-pos-1.xyz".format(output_file_hierarchy_directory,file,perfect_project_name)) #¡***¡
            new_xyz_file = str("{}{}{}-pos-L.xyz".format(output_file_hierarchy_directory,file,perfect_project_name))
        else:
            xyz_file = str("{}{}{}{}-pos-1.xyz".format(output_file_hierarchy_directory,file,defect_project_name_stem,X)) #¡***¡
            new_xyz_file = str("{}{}{}{}-pos-L.xyz".format(output_file_hierarchy_directory,file,defect_project_name_stem,X))
            
        file = open(xyz_file, 'r') #¡***¡
        index = 0 #¡***¡
        for line in file: #¡***¡
            index += 1 #¡***¡
            if tot_atoms in line: #¡***¡
                j = index #¡***¡
                itr_start.append(j) #¡***¡ ## create array with line number of the each interation's starting line
        l = (itr_start[-1] - 1) #¡***¡ ## record line number of the last interation's starting line
        # print(X,state,total_atoms_in_calc[-1]) #¡***¡
        for n in total_atoms_in_calc: #¡***¡
            Lns = l + n #¡***¡
            L_itr_lines.append(Lns) #¡***¡     
        file.close #¡***¡
        file = open(xyz_file, 'r') #¡***¡
        output_file = open(new_xyz_file, 'w') #¡***¡
        for position, line in enumerate(file): #¡***¡ ## enumerate used to search file and return the text string found on line numbers specified by position
            if position in L_itr_lines: #¡***¡
                string = line #¡***¡
                output_file.write(string) #¡***¡ ## text string found on file's specified line numbers wrote to output_file

        exec(f'{system}_{state}_x, {system}_{state}_y, {system}_{state}_z = np.loadtxt(new_xyz_file, skiprows=2, usecols=(1,2,3), unpack=True)') #!*!

for X in '{...}': #!*! # list all defects in material being studied. If defects are related to addition of impurity atoms, list the symbols of impurities used in the different calculations rather than defect name
    system = str(X) #!*!
    for state in '{...}': #!*! # charge state of defects; for only one charge state comment out line and write replacement line state = {...} below. Remember to move code lines included within for state in {...} 
                        # back an indent which can be done by highlighting code and pressing shift tab. 
        # exec(f'{system}_{state}_defect_site_x = {...}') #!*! # if defect is vacancy, write x-coordinate of removed atom
        # exec(f'{system}_{state}_defect_site_y = {...}') #!*! # if defect is vacancy, write y-coordinate of removed atom
        # exec(f'{system}_{state}_defect_site_z = {...}') #!*! # if defect is vacancy, write z-coordinate of removed atom 
        for d in 'x', 'y', 'z': #!*!
            exec(f'{system}_{state}_diff_{d} = perf_{perfect_project_name}_{...}_{d} - {system}_{state}_{d}') #!*! # if state changed for perfect structure above to stop any iteration errors write what state was changed to in {...}
            exec(f'{system}_{state}_defect_site_{d} = {system}_{state}_{d}[{...}]') #!*! # if defect is a substitional or interstitial defect, write atom index of subsituted/inserted atom. If defect is a vacancy comment 
                                                                                    # out line and uncomment out the three lines prior to the for loop of d.
            exec(f'{system}_{state}_distance_{d} = {system}_{state}_defect_site_{d} - {system}_{state}_{d}') #!*!
        exec(f'{system}_{state}_tot_displacement = np.sqrt({system}_{state}_diff_x**2 + {system}_{state}_diff_y**2 + {system}_{state}_diff_z**2)') #!*! ## vector magnitude of displacement from perfect structure geometry
        exec(f'{system}_{state}_tot_distance = np.sqrt({system}_{state}_distance_x**2 + {system}_{state}_distance_y**2 + {system}_{state}_distance_z**2)') #!*! ## vector magnitude of distance from defect site
        exec(f'{system}_{state}_tot_distance_sorted = np.sort({system}_{state}_tot_distance)') #!*! ## distances from defect site sorted from smallest to largest 
        exec(f'{system}_{state}_tot_displacement_sorted = [x for _, x in sorted(zip({system}_{state}_tot_distance,{system}_{state}_tot_displacement))]') #!*! ## displacements from perfect structure geometry sorted based  
                                                                                                                                                         ## on distance from defect site from smallest to largest
        exec(f'{system}_{state}_tot_displacement_sorted2 = np.sort({system}_{state}_tot_displacement)') #!*! ## displacements from perfect structure geometry sorted from smallest to largest