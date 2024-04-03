import copy
import importlib

import numpy as np
from numpy import genfromtxt

from cplap_with_beef import Utility
importlib.reload(Utility)

from cplap_with_beef import Debug
importlib.reload(Debug)

"""
The purpose of this file is to convert the change in chemical potential..
as calculated from Process_BEEF_Entry, and to create a complete chemical...
potential, \mu_final, where \mu_final = (\mu_ref + \delta\mu_FERE + \mu_BEEF)...
where \mu_ref is the per atom chemical potential of the reference element state...
(with BEEF there will be 2000 \mu_ref values), \delta\mu_FERE is the change in mu...
(this is a constant value for each element), and \mu_BEEF is the value generated via..
Process_BEEF_Entry.

Note: If the use of FERE values is requested, but corresponding element values are not found, the 
code will crash. Not elegant, but better to fail loudly.

"""
class Generate_Final_Values:

    def __init__(self,debug_obj=False):

        #use the debug_obj that has been passed in, or make a new one
        self.debug_obj = False
        if debug_obj != False and isinstance(debug_obj,object):
            self.debug_obj = debug_obj
        else:
            self.debug_obj = Debug.Debug(debug_bool = True,restart = True)

        self.util_obj = Utility.Utility(debug_obj=self.debug_obj)

        self.species_order = False
        self.BEEF_chemical_potential_reference_dict = {}

        self.FERE_values_by_element = {
            'Li': 0.21,
            'Be': 0.35,
            'N': -0.20,
            'O': 0.23,
            'Ga': 0.66,
            'In': 0.41,
            'Ba': 0.54,
            'Sn': 0.18,
            'I': 0,
            'Pb':0,
        }

    def get_FERE_value_by_element(self,element):
        element_found = False
        FERE_value = 0
        if element in self.FERE_values_by_element.keys():
            element_found = True
            FERE_value = self.FERE_values_by_element[element]
        return element_found,FERE_value

    def set_species_order(self,species_order):
        self.species_order = species_order

    def organize_final_chem_pot_calc(self):


        if self.species_order == False:
            self.debug_obj.debug_log("The species order must be set to generate final chemical potentials",type="ERROR")
            return False

        for specie in self.species_order:
            if self.check_if_BEEF_chemical_potential_reference_available(specie):
                continue
            else:
                self.debug_obj.debug_log(f"No raw energy data has been given for {specie}",type="Error")
                return False

        for beef_idx in range(0,2001):

            self.util_obj.changeDirectory("grid_exact")
            file_found,beef_data = self.read_in_beef_file(beef_idx)

            grid_shape = (0,len(self.species_order))
            grid_total_cp = np.zeros(grid_shape)

            if file_found:
                num_rows = len(beef_data)
                for row_index in range(0,num_rows):

                    temp_row = []

                    for col_index in range(0,len(self.species_order)):
                        specie = self.species_order[col_index]

                        #Collect all relevant chemical potenital info
                        delta_mu = beef_data[row_index][col_index]
                        mu_FERE_adjustment =  self.FERE_values_by_element[specie]
                        mu_reference = self.get_BEEF_chemical_potential_reference(specie,beef_idx)
                        """
                        print(f"row: {row_index},"
                              f"col: {col_index},"
                              f"delta_mu: {delta_mu},"
                              f"mu_FERE: {mu_FERE_adjustment}"
                              f"mu_ref: {mu_reference}")
                        """
                        total_cp = delta_mu + mu_reference + mu_FERE_adjustment
                        temp_row.append(total_cp)


                    grid_total_cp = np.append(grid_total_cp,[temp_row],axis=0)

                self.util_obj.change_directory_auto("grid_total")
                filename = f"beef_{beef_idx}_grid_total.csv"
                np.savetxt(filename, grid_total_cp, delimiter=",")

    # Read in a single chemical potential set for a given beef index
    # Parameter | beef_idx = the beef index
    def read_in_beef_file(self,beef_idx):

        generic_filename = "beef_#_exact.csv"

        precise_filename = generic_filename.replace("#",str(beef_idx))

        self.util_obj.change_directory_auto(name_of_desired_folder="grid_exact")
        file_found = self.util_obj.checkIfFileExists(precise_filename)
        #print(f"Beef {beef_idx} search result: {file_found}")

        if file_found:
            beef_data = genfromtxt(precise_filename, delimiter=',')
            return file_found,beef_data
        else:
            return False,False


    def set_BEEF_chemical_potential_reference(self,specie,list_of_raw_energies):
        self.BEEF_chemical_potential_reference_dict[specie] = list_of_raw_energies

    def get_BEEF_chemical_potential_reference(self,specie,beef_idx):
        raw_energies = self.BEEF_chemical_potential_reference_dict[specie]
        return raw_energies[beef_idx]

    def check_if_BEEF_chemical_potential_reference_available(self,specie):
        if specie in self.BEEF_chemical_potential_reference_dict.keys():
            return True
        else:
            return False



