
import copy
import importlib
import os
import numpy as np

from cplap_with_beef import Utility
importlib.reload(Utility)

from cplap_with_beef import Debug
importlib.reload(Debug)

from cplap_with_beef import QE_Helper
importlib.reload(QE_Helper)

from cplap_with_beef import CPLAP_Helper
importlib.reload(CPLAP_Helper)

from cplap_with_beef import Process_BEEF_Entry
importlib.reload(Process_BEEF_Entry)

from cplap_with_beef import Data_Analyzer
importlib.reload(Data_Analyzer)

from cplap_with_beef import Generate_Final_Values
importlib.reload(Generate_Final_Values)

class Organizer:

    # Primary object for user to interact with
    # Parameter| grid_spacing - how tight the grid will be for all results
    #            If too tight, code will be slow, especially when multiple species.
    #            If too loose, some beef indeces might be considered unstable as no point can 'fit' within envelope
    def __init__(self,debug_bool,grid_spacing = 0.01):

        self.phase_info = []
        self.beef_results = []

        #expected type
        self.expected_species_info = {}

        #use the debug_obj that has been passed in, or make a new one
        self.debug_obj = Debug.Debug(debug_bool = debug_bool,restart = True)
        self.util_obj = Utility.Utility(debug_obj=self.debug_obj)
        self.analyzer = Data_Analyzer.Data_Analyzer(debug_obj=self.debug_obj)
        self.final_cp_calculator = Generate_Final_Values.Generate_Final_Values(debug_obj=self.debug_obj)
        self.cplap_helper = None

        self.generate_folder_structure()

        self.grid_spacing = grid_spacing


    # Add different phases to the system. One phase per function call.
    # Parameter | type = 'primary','competing','element'
    # species_info | dictionary. Key = number. Value = { 'number': [number of atomic species],'symbol': [atomic type]}
    #                e.g. species_info = {'1': {'number': 1,'symbol': 'Ba'},
    #                                     '2': {'number': 1,'symbol': 'Sn'},
    #                                     '3': {'number': 3,'symbol': 'O'},
    #File path to quantum espresso *.pwo file of completed ensemble calculation
    def addPhase(self,type,species_info,dft_output_path):

        self.debug_obj.debug_log(f"Attempting to add phase from file: {dft_output_path}")

        temp_dict = {'type': type,
                     'species_info': species_info,
                     'dft_output_path': dft_output_path,
                     'raw_energies_per_atom': False,
                     'chemical_potentials': False}

        qeh = QE_Helper.QE_Helper(self.debug_obj)
        qe_results = qeh.read_qe_output(dft_output_path)

        acceptable_types = ['primary','competing','element']
        if type not in acceptable_types:
            self.debug_obj.debug_log(f"The provided phase type ({type}) is not allowed.",type="ERROR")
            return False

        # Check that the dft calculation ran as expected
        if qe_results['job_done'] == False:
            self.debug_obj.debug_log("The provided dft file indicates the job did not complete.",type="ERROR")
            return False

        # Check that the dft calc is consistent with the expected species ratio
        num_of_mocules_in_output = self.find_number_of_molecules_in_output(species_info,qe_results['atom_dict'])
        if num_of_mocules_in_output == False:
            self.debug_obj.debug_log("The provided dft file does not match the expected species.",type="ERROR")
            return False

        number_of_atoms_per_molecule = 0
        for atom_symbol in species_info:
            number_of_atoms_per_molecule += species_info[atom_symbol]

        #Add a check that if an element is being added that it only has one atom
        if type == 'element' and number_of_atoms_per_molecule != 1:
            self.debug_obj.debug_log("Element type phases must only have one atom type specified.",type="ERROR")
            return False

        #the beef ensemble represents the change w.r.t. to total energy. Must combine.
        raw_energies_per_atom = np.array(qe_results['beef_ensemble']) + qe_results['total_energy']

        #Add an EXTRA energy for the optimal result, which is not included in the BEEF ensemble
        raw_energies_per_atom = np.append(raw_energies_per_atom,qe_results['total_energy'])
        print(f"Length of raw_energies_per_atom: {len(raw_energies_per_atom)}")

        raw_energies_per_atom = raw_energies_per_atom/(number_of_atoms_per_molecule*num_of_mocules_in_output)
        temp_dict['raw_energies_per_atom'] = raw_energies_per_atom

        self.phase_info.append(temp_dict)

    def organize_CPLAP_interaction(self,dependent_species):

        self.cplap_helper = CPLAP_Helper.CPLAP_Helper(self.phase_info,dependent_species,debug_obj=self.debug_obj)
        self.cplap_helper.generate_cplap_input()
        self.cplap_helper.run_CPLAP_input()
        beef_list = self.cplap_helper.parse_CPLAP_output()
        self.debug_obj.debug_log(f"BEEF list length: {len(beef_list)}")

        #create a list of beef processor objects
        print(f"Lenth of BEEF list @ organize_CPLAP_interaction: {len(beef_list)}")

        for beef_entry in beef_list:

            beef_num = beef_entry['beef_num']
            temp_beef_processor = Process_BEEF_Entry.Process_BEEF_Entry(beef_num,self.grid_spacing,debug_obj=self.debug_obj)

            temp_beef_processor.set_list_of_inequalities(beef_entry['inequality_conditions'])
            temp_beef_processor.set_dependent_variable_calc(beef_entry['dependent_calc'])

            intersections_found = temp_beef_processor.set_list_of_intersections(beef_entry['intersections'])

            #only add to results if intersections exist (i.e. it's stable)
            if intersections_found == True:
                self.beef_results.append(temp_beef_processor)

    def process_all_BEEF_results(self):

        print(f"Processing {len(self.beef_results)} BEEF results\n")

        delete_list = []

        count = 1
        for beef_result in self.beef_results:
            print(f"count: {count}",end="\r")
            result = beef_result.organize_grid_densify()
            if result == False:
                failed_idx = beef_result.beef_idx
                self.debug_obj.debug_log(f"Beef failed: {failed_idx} ")
                delete_list.append(beef_result.beef_idx)
            count += 1

        for failed_idx in delete_list:
            for beef_result in self.beef_results:
                if failed_idx == beef_result.beef_idx:
                    self.beef_results.remove(beef_result)

    def process_BEEF_result(self,beef_idx):

        print(f"Processing BEEF {beef_idx}\n")

        for beef_result in self.beef_results:
            if beef_result.beef_idx == beef_idx:
                result = beef_result.organize_grid_densify()
                if result == False:
                    failed_idx = beef_result.beef_idx
                    self.debug_obj.debug_log(f"Beef failed: {failed_idx} ")
                    self.beef_results.remove(beef_result)

        """
        print(f"Beef #{beef_idx} "
              f"Intersections:")
        print(beef_result.get_list_of_intersections())
        print(f"Valid Intersection Neighbors: \n"
              f"{beef_result.get_valid_grid_points()}")
        """

    def generate_final_chemical_potentials(self):
        species_order = self.cplap_helper.get_species_order()
        self.final_cp_calculator.set_species_order(species_order)

        for specie in species_order:
            list_of_raw_energies = self.get_list_of_raw_energies(specie)
            self.final_cp_calculator.set_BEEF_chemical_potential_reference(specie,list_of_raw_energies)

        self.final_cp_calculator.organize_final_chem_pot_calc()


    # Search through the dictionary of phases for the matching atomic specie that is an element
    # Using the data from the quantum espresso calculation, return the raw energy per atom
    # Parameter | specie - the atomic specie of interest (e.g. 'Ba')
    def get_list_of_raw_energies(self,specie):
        for phase in self.phase_info:
            if phase['type'] == "element":

                number_of_species = len(phase['species_info'])
                first_specie = list(phase['species_info'].keys())[0]

                if number_of_species == 1 and first_specie == specie:
                    return phase['raw_energies_per_atom']

    def agregrate_all_beef_entries(self):
        self.analyzer.agregrate_beef_info()

    def plot_beef_entries(self):
        self.analyzer.generate_2d_heat_map()

    #Once all phases are added calculate the actual chemical potentials
    def generate_chemical_potentials(self):
        """
        temp_dict = {'type': type,
             'species_info': species_info,
             'dft_output_path': dft_output_path,
             'raw_energies_per_atom': False,
             'chemical_potentials': False}
        """
        element_dict = self._generate_chemical_potentials_of_elements()

        #phase_info is a list
        for phase in self.phase_info:

            if phase['type'] == 'primary' or phase['type'] == 'competing':

                number_of_atoms_per_molecule = 0
                for atom_symbol in phase['species_info']:
                    number_of_atoms_per_molecule +=  phase['species_info'][atom_symbol]

                phase['chemical_potentials'] = copy.deepcopy(phase['raw_energies_per_atom'])
                phase['chemical_potentials'] = phase['chemical_potentials']*number_of_atoms_per_molecule

                for atom_symbol in phase['species_info']:
                    atom_quantity = phase['species_info'][atom_symbol]
                    phase['chemical_potentials'] -= element_dict[atom_symbol]*atom_quantity

    def _generate_chemical_potentials_of_elements(self):
        #phase_info is a list
        element_dict = {}

        for phase in self.phase_info:

            if phase['type'] == 'element':

                #print(phase['species_info'])
                #print(phase['species_info'].keys())
                element_symbol = list(phase['species_info'].keys())[0]
                #element_symbol = phase['species_info'].keys()[0]
                raw_energies = phase['raw_energies_per_atom']
                element_dict[element_symbol] = raw_energies

        return element_dict



    def find_number_of_molecules_in_output(self,expected_dict,actual_dict):

        self.debug_obj.debug_log(f"expected dict: {expected_dict}")
        self.debug_obj.debug_log(f"actual_dict: {actual_dict}")

        atom_multiples = []
        for atom_symbol_expected in expected_dict:

            expected_num = expected_dict[atom_symbol_expected]
            actual_num = 0

            try:
                actual_num = actual_dict[atom_symbol_expected]
            except:
                self.debug_obj.debug_log(f"The atom type {atom_symbol_expected} was not found in the calculation.",type="WARNING")

            multiple = actual_num/expected_num
            atom_multiples.append(multiple)

        def check(list):
            return all(i == list[0] for i in list)

        if check(atom_multiples) == True:
            return atom_multiples[0]
        else:
            return False

    #This should return the chemical potential of the 'optimal' XC Functional
    #the 'optimal' result of the DFT calcs is not included in the BEEF ensemble
    #As a result, the 'optimal' energy is currently added to the end of the list of 2000
    #BEEF energies. This is hacky, improve later
    # Parameter |
    def get_optimal_chemical_potential(self,find_centroid = False,get_total_cp = False):
        species_order = self.cplap_helper.get_species_order()
        self.analyzer.set_species_order(species_order)

        return self.analyzer.get_packaged_grid_points(start_index = 2000,
                                                      stop_index = 2001,
                                                      find_centroid=find_centroid,
                                                      get_total_cp=get_total_cp)

    # This returns a dictionary of lists of dictionaries
    # e.g. {0: [{'Ga': 1.0, 'N': 2.0},{'Ga': 2.0, 'N': 3.0}],
    #       1: [{'Ga': 1.1, 'N': 2.0},{'Ga': 2.0, 'N': 3.0}}
    # The top level dictionary key indictaes the BEEF index
    # The top-level dictionary value (the list of dictionary) holds all valid grid points for that BEEF index
    def get_beef_ensemble_chemical_potentials(self,find_centroid = False,get_total_cp = False):
        #Get the species order and provide it to Process_BEEF_entry
        species_order = self.cplap_helper.get_species_order()
        self.analyzer.set_species_order(species_order)

        return self.analyzer.get_packaged_grid_points(start_index = 0,
                                                      stop_index = 2000,
                                                      find_centroid = find_centroid,
                                                      get_total_cp = get_total_cp)

    # Create the folders for storing all intermediate data
    def generate_folder_structure(self):

        main = "chemical_potential_info"
        input_on_deck = "input_on_deck"
        input_stored = "input_stored"
        results_stored = "results_stored"
        grid_exact = "grid_exact"
        grid_total = "grid_total"

        self.util_obj.createNewDirectory(main)
        self.util_obj.changeDirectory(main)

        path_main = os.getcwd()

        path_input_on_deck = f"{path_main}/{input_on_deck}"
        path_input_stored = f"{path_main}/{input_stored}"
        path_results_stored = f"{path_main}/{results_stored}"
        path_grid_exact = f"{path_main}/{grid_exact}"
        path_grid_total = f"{path_main}/{grid_total}"

        self.util_obj.createNewDirectory(path_input_on_deck)
        self.util_obj.createNewDirectory(path_input_stored)
        self.util_obj.createNewDirectory(path_results_stored)
        self.util_obj.createNewDirectory(path_grid_exact)
        self.util_obj.createNewDirectory(path_grid_total)

        self.util_obj.changeDirectory(path_main)

        filename_of_paths = ".file_paths.txt"

        # create a gaussian fit input file
        data = f"main = {path_main}\n" \
               f"input_on_deck = {path_input_on_deck}\n" \
               f"input_stored = {path_input_stored}\n" \
               f"results_stored = {path_results_stored}\n" \
               f"grid_exact = {path_grid_exact}\n" \
               f"grid_total = {path_grid_total}\n"

        with open(filename_of_paths, 'w') as f:
            f.write(data)
            f.close()








#%%