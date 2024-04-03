import copy
import importlib
import os

from cplap_with_beef import Utility
importlib.reload(Utility)

from cplap_with_beef import Debug
importlib.reload(Debug)

# Manages all interaction with quantum espresso *.pwo file
class QE_Helper:

    def __init__(self,debug_obj=False):

        #use the debug_obj that has been passed in, or make a new one
        self.debug_obj = False
        if debug_obj != False and isinstance(debug_obj,object):
            self.debug_obj = debug_obj
        else:
            self.debug_obj = Debug.Debug(debug_bool = True,restart = True)

        self.util_obj = Utility.Utility(debug_obj=self.debug_obj)

        #Search through the QE output file

        #check that the QE info matches the expected species info
        #Look for "Cartesian axes"
        #Then started iterating through the elements and count them up "Crystallographic axes"

        #Grab 2000 energies
        #Look for "BEEFens 2000 ensemble energies"
        #Save all energies

        #Update energies
        #Divide by total number of atoms
        #save chemical potential as value per atom (lets say chp_per_atom)
        #this means for the co of a molecule would be (chp_per_atom*number of atoms in molecule)

    #Read through the qe output. Return a dictionary of raw qe info
    def read_qe_output(self,filename):

        result = {
            'job_done':False,
            'total_energy': False,
            'beef_ensemble':False,
            'atom_dict':False,
        }

        self.debug_obj.debug_log(f"Searching for {filename}")
        file_found = self.util_obj.checkIfFileExists(filename)
        self.debug_obj.debug_log(f"file_found: {file_found}")

        if file_found == True:

            self.debug_obj.debug_log(f"Current directory: {os.getcwd()}")
            self.debug_obj.debug_log(f"The file '{filename}' has been found! Huzzah!")

            job_complete = False
            total_energy_found = False
            total_energy = 0

            line_nums_atoms = [0,0]
            line_nums_beef_ensemble = [0,0]

            unit_type = ""
            ry_to_eV = 13.6057039763


            with open(filename) as file:

                line_count = 0

                line_info_cart_axes = 'Cartesian axes'
                line_info_crystal_axes = 'Crystallographic axes'
                line_info_beef_ens_header = 'BEEFens 2000 ensemble energies'
                line_info_min_energy = '!    total energy'
                line_info_job_done = 'JOB DONE.'

                #Loop through every line in the output file
                for line in file:

                    line_count += 1

                    #1. Search for minimum energy
                    if line_info_min_energy in line:
                        total_energy_found = True

                        unit = line.split(' ')[-1]
                        #self.debug_obj.debug_log(f"Energy unit: {unit}")

                        total_energy = float(line.split(" ")[-2])
                        #self.debug_obj.debug_log(f"Initial energy value: {total_energy}")

                        if "Ry" in unit:
                            #self.debug_obj.debug_log(f"Energy in Ry!")
                            unit_type = "Ry"
                            total_energy = total_energy*ry_to_eV

                    #2. Search for jobe done
                    if line_info_job_done in line:
                        job_complete = True
                        self.debug_obj.debug_log(f"Job complete found!")

                    #3. Search for 'Cartesian axes' header
                    if line_info_cart_axes in line:
                        line_nums_atoms[0] = line_count + 3

                    #4. Search for 'Crystallographic axes' header
                    if line_info_crystal_axes in line:
                        line_nums_atoms[1] = line_count - 2

                    #5. Search for 'BEEFens 2000 ensemble energies'
                    if line_info_beef_ens_header in line:
                        line_nums_beef_ensemble[0] = line_count + 1
                        line_nums_beef_ensemble[1] = line_count + 2000

                if total_energy_found and job_complete:
                    result['job_done'] = True
                    result['total_energy'] = total_energy
                    self.debug_obj.debug_log(f"Total energy = {total_energy} eV!")
                else:
                    return result

            #print(f"line_nums_atoms: {line_nums_atoms}")
            #print(f"line_nums_beef_ensemble: {line_nums_beef_ensemble}")

            #Now grab the atoms numbers and beef ensemble
            with open(filename) as file:

                line_count = 0

                atom_dict = {}
                beef_ens = []

                for line in file:
                    line_count += 1

                    if line_nums_atoms[0] <= line_count <= line_nums_atoms[1]:
                        atom_symbol = line.split()[1]
                        if atom_symbol in atom_dict:
                            atom_dict[atom_symbol] += 1
                        else:
                            atom_dict[atom_symbol] = 1

                    if line_nums_beef_ensemble[0] <= line_count <= line_nums_beef_ensemble[1]:
                        raw_energy = float(line.split()[0])
                        clean_energy = 0
                        if unit_type == "Ry":
                            clean_energy = raw_energy*ry_to_eV
                        else:
                            clean_energy = raw_energy
                        beef_ens.append(clean_energy)

                    if line_count > line_nums_beef_ensemble[1]:
                        break

                result['beef_ensemble'] = copy.deepcopy(beef_ens)
                result['atom_dict'] = copy.deepcopy(atom_dict)

        return result

#%%
