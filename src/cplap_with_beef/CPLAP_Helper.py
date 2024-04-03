import copy
import importlib
import os
import subprocess

import numpy as np

import re

from cplap_with_beef import Utility
importlib.reload(Utility)

from cplap_with_beef import Debug
importlib.reload(Debug)

# The CPLAP code by Buckeridge et al. is essential for CPLAP_With_BEEF.
# This class works directly with CPLAP.
class CPLAP_Helper:
    """
    param: debug_obj
    param: desired_species_order - the last item will become the dependent variable
    """
    def __init__(self,phase_info,dependent_species,debug_obj=False):


        #use the debug_obj that has been passed in, or make a new one
        self.debug_obj = False
        if debug_obj != False and isinstance(debug_obj,object):
            self.debug_obj = debug_obj
        else:
            self.debug_obj = Debug.Debug(debug_bool = True,restart = True)

        self.util_obj = Utility.Utility(debug_obj=self.debug_obj)

        self.input_on_deck = "input_on_deck"
        self.input_stored = "input_stored"
        self.results_stored = "results_stored"

        self.input_file_name_generic = "cplap_input_beef_#.dat"
        self.output_file_name_generic = "cplap_output_beef_#.dat"

        self.cplap_input_filename = "input.dat"
        self.cplap_results_filename = "results.dat"

        self.species_regex = 'mu_([A-Za-z]+)'
        self.intersection_value_regex = '(-?[0-9]+.[0-9]+)'

        self.phase_info = phase_info
        self.dependent_species = dependent_species
        self.independent_species = self.generate_species_order()

        self.debug_obj.debug_log(f"The dependent species: {self.dependent_species}")
        self.debug_obj.debug_log(f"The independent species: {self.independent_species}")

    # Create the desired order of the independent species
    # This is important to ensure that all future calculations are doing calculations in the
    # right order
    def generate_species_order(self):

        for phase in self.phase_info:

            if phase['type'] == "primary":

                primary_species = phase['species_info']
                species_order = list(primary_species.keys())

                if self.dependent_species in species_order:
                    species_order.remove(self.dependent_species)
                else:
                    self.debug_obj.debug_log("The provided dependent variable doesn't exist in the primary phase.",type="ERROR")
                    return False

                return species_order

    # Will return the correct order of the INDEPENDENT species, followed by the single DEPENDENT species
    def get_species_order(self):
        species_order = copy.deepcopy(self.independent_species)
        species_order.append(self.dependent_species)
        return species_order

    # For each BEEF-indexed DFT energy, create the input to use CPLAP
    # to determine the chemical potential stability envelope
    def generate_cplap_input(self):

        self.util_obj.change_directory_auto(self.input_on_deck)
        beef_ensemble_num = 2001

        for beef_num in range(0,beef_ensemble_num):

            #First find the primary info
            primary_species_num = 0
            primary_species_str = ""

            competing_total_num = 0

            competing_species_num = []
            competing_stoich_str = []

            #First go through and make the primary info
            for phase in self.phase_info:

                if phase['type'] == 'primary':
                    primary_species_num = len(list(phase['species_info']))

                    for specie in phase['species_info']:
                        number = phase['species_info'][specie]
                        primary_species_str = f"{primary_species_str}{number} {specie} "

                    primary_species_str = f"{primary_species_str}{phase['chemical_potentials'][beef_num]}"

                if phase['type'] == 'competing':
                    competing_total_num += 1

                    #Add the number of species to the list of numbers
                    temp_num = len(list(phase['species_info']))
                    competing_species_num.append(temp_num)

                    temp_str = ""

                    for specie in phase['species_info']:
                        number = phase['species_info'][specie]
                        temp_str = f"{temp_str}{number} {specie} "

                    temp_str = f"{temp_str}{phase['chemical_potentials'][beef_num]}"
                    competing_stoich_str.append(temp_str)


            #Print info
            self.debug_obj.debug_log(f"primary_species_num: {primary_species_num}")
            self.debug_obj.debug_log(f"primary_species_str: {primary_species_str}")
            self.debug_obj.debug_log(f"competing_total_num: {competing_total_num}")
            self.debug_obj.debug_log(f"competing_species_num: {competing_species_num}")
            self.debug_obj.debug_log(f"competing_stoich_str: {competing_stoich_str}")

            #Now create the full file
            data = f"#Number of species in primary phase\n" \
                   f"{primary_species_num}\n" \
                   f"#Primary phase relationship\n" \
                   f"{primary_species_str}\n" \
                   f"#Dependent variable or 'none\n" \
                   f"{self.dependent_species}\n" \
                   f"#Number of competing phases\n" \
                   f"{competing_total_num}\n" \
                   f"#number of species in competing phases followed by relationship\n"

            for i in range(0,competing_total_num):
                data = f"{data}{competing_species_num[i]}\n" \
                       f"{competing_stoich_str[i]}\n"

            cplap_input_filename = self.input_file_name_generic.replace("#",str(beef_num))

            with open(cplap_input_filename,"w") as f:
                f.write(data)
                f.close()

    # Work through the list of created CPLAP input files and run them
    def run_CPLAP_input(self,generate_grid = False,grid_increment = 0.01):

        self.input_file_pattern = '[a-z\_]+([0-9]+).dat'

        self.util_obj.change_directory_auto("input_on_deck")

        item_list = os.listdir()

        for item in item_list:

            try:
                a = re.search(self.input_file_pattern, item)
                beef_num = int(a.group(1))

                if isinstance(beef_num,int):
                    self.run_CPLAP_for_beef_num(beef_num,generate_grid,grid_increment)
            except:
                pass

    # CPLAP has a specific location where all files must be read from or created
    # This function copies files from input_on_deck into the CPLAP location
    # It then runs CPLAP and places the results in results_stored
    # The input is stored in input_stored
    def run_CPLAP_for_beef_num(self,beef_num,generate_grid = False,grid_increment = 0.01):

        # Check that there is a valid CPLAP environmental variable
        # Should be the location, not the executable

        loc_CPLAP_DIR = self.util_obj.checkEnvironmentalVars("CPLAP_DIR")

        address_dict = self.util_obj.fetch_directory_dict()

        if address_dict == False:
            self.debug_obj.debug_log("Address file not found. Cannot proceed.",type="WARNING")

        loc_main = ""
        loc_input_on_deck = ""
        loc_input_stored = ""
        loc_results_stored = ""
        summer = 0
        #fetch the all required locations for folder io
        #This logic shouldn't be here. Move to Utility
        for loc in address_dict:
            if loc == "main":
                loc_main = address_dict[loc]
                summer += 1
            if loc == "input_on_deck":
                loc_input_on_deck = address_dict[loc]
                summer += 1
            if loc == "input_stored":
                loc_input_stored = address_dict[loc]
                summer += 1
            if loc == "results_stored":
                loc_results_stored = address_dict[loc]
                summer += 1

        if summer != 4:
            self.debug_obj.debug_log("Some required folder locations could not be found",type="ERROR")

        try:
            #Change to CPLAP dir and delete old files
            self.util_obj.changeDirectory(loc_CPLAP_DIR)
            self.util_obj.deleteFiles("input.dat",precise = True)
            self.util_obj.deleteFiles("results.dat",precise = True)
            self.util_obj.deleteFiles("grid.dat",precise = True)

            #copy input_beef_x to CPLAP dir as input.dat
            input_on_deck_filename = self.input_file_name_generic.replace("#",str(beef_num))

            src_path = f"{loc_input_on_deck}/{input_on_deck_filename}"
            dst_path = f"{loc_CPLAP_DIR}/{self.cplap_input_filename}"
            self.util_obj.copyFiles(src_path,dst_path,delete_old = False) #Don't delete yet!

            # run CPLAP!
            # Will need to add more logic later to make grid with CPLAP
            if generate_grid == False:
                os.system("echo n | ./CPLAP > output")
                #subprocess.run(["echo","n","|","./CPLAP"])

            valid_results_string = "Results for system"
            results_are_valid = False

            #Check that the file isn't empty or gibberish. A rough check.
            with open(self.cplap_results_filename) as file:
               for line in file:
                   if valid_results_string in line:
                       results_are_valid = True

            #Check results are complete
            if results_are_valid == False:
                self.debug_obj.debug_log(f"CPLAP results file makes no sense for BEEF #{beef_num}",type="WARNING")
                #return False

            #Move results to output_stored. Rename to include beefnum. Delete old file.
            output_stored_filename = self.output_file_name_generic.replace("#",str(beef_num))
            src_path = f"{loc_CPLAP_DIR}/{self.cplap_results_filename}"
            dst_path = f"{loc_results_stored}/{output_stored_filename}"
            result_stored_filepath = dst_path
            self.util_obj.copyFiles(src_path,dst_path,delete_old = False)

            #If transfer of results is successful, delete source file
            result_is_stored = self.util_obj.checkIfFileExists(result_stored_filepath)
            if result_is_stored == True:
                self.util_obj.deleteFiles(src_path,precise=True)

            #Move input from input_on_deck to input_stored
            input_on_deck_filename = self.input_file_name_generic.replace("#",str(beef_num))
            src_path = f"{loc_input_on_deck}/{input_on_deck_filename}"
            dst_path = f"{loc_input_stored}/{input_on_deck_filename}"
            input_stored_filepath = dst_path
            self.util_obj.copyFiles(src_path,dst_path,delete_old = False)

            #If the transfer of the input is successful, delete source file
            input_is_stored =  self.util_obj.checkIfFileExists(input_stored_filepath)
            if input_is_stored == True:
                self.debug_obj.debug_log("Now attempting to delete input from on_deck")
                self.debug_obj.debug_log(f"src path: {src_path}")
                self.util_obj.deleteFiles(src_path,precise=True)

            #remove input from CPLAP
            self.util_obj.changeDirectory(loc_CPLAP_DIR)
            self.util_obj.deleteFiles("input.dat",precise = True)
            self.util_obj.deleteFiles("results.dat",precise = True)
            self.util_obj.deleteFiles("grid.dat",precise = True)

            #Go back to main folder
            self.util_obj.changeDirectory(loc_main)

            #If all good, return success
            if input_is_stored and result_is_stored:
                return True
            else:
                return False

        except:
            self.util_obj.changeDirectory(loc_main)
            self.debug_obj.debug_log(f"An unexpected error occurred.",type="WARNING")
            return False

    # Read the CPLAP results for ALL beef indeces
    def parse_CPLAP_output(self):

        self.debug_obj.debug_log("Made it to parse_CPLAP_output")

        self.input_file_pattern = '[a-z\_]+([0-9]+).dat'

        self.util_obj.change_directory_auto("results_stored")

        item_list = os.listdir()

        beef_list = []
        """
        num_to_run = 2000 #num_to_run
        
        for i in range(0,num_to_run):

            item = item_list[i]
        """
        for item in item_list:

            #try:
            print(f"item: {item}")
            a = re.search(self.input_file_pattern, item)
            beef_num = int(a.group(1))

            if isinstance(beef_num,int):

                self.debug_obj.debug_log(f"beef_num: {beef_num}")
                beef_entry = self.parse_CPLAP_output_for_beef_num(beef_num)

                beef_entry['dependent_calc'] = self.generate_dependent_variable_calc_func(beef_num)
                beef_list.append(beef_entry)
            #except:
                #self.debug_obj.debug_log("An exception was caught in parse_CPLAP_output")
                #pass

        return beef_list

    #Parse a SINGLE CPLAP result corresponding to a SINGLE BEEF index
    def parse_CPLAP_output_for_beef_num(self,beef_num):

        self.debug_obj.debug_log("Made it to parse_CPLAP_output for beef_num")

        beef_entry = {}
        beef_entry['beef_num'] = beef_num

        self.util_obj.change_directory_auto('results_stored')

        self.debug_obj.debug_log("Changed directory successfully")

        #find file that matches beef num
        output_filename = self.output_file_name_generic.replace('#',str(beef_num))
        self.debug_obj.debug_log(f"output_filename: {output_filename}")

        output_exists = self.util_obj.checkIfFileExists(output_filename)

        if output_exists == False:
            self.debug_obj.debug_log(f"file {output_filename} not found")
            return False

        system_stable = False
        system_binary = False

        line_loc_limiting_inequalities = -1
        line_loc_intersecting_points = -1
        line_loc_binary_points = -1
        line_loc_unstable = -1
        line_loc_competing_phase = -1

        marker_limiting_inequalities = "Limiting inequalities:"
        marker_intersecting_points = "Intersection points in chemical potential space:"
        marker_binary_points = "Solution for binary system"
        marker_unstable = "System is unstable - no solutions exist"
        marker_competing_phase = "Competing phase   1"

        line_locs_inequalities = [0,0]
        line_locs_intersections = [0,0]

        cplap_species_order = [] #list of strings
        beef_inequalities = [] #list of functions which will return true/false give certain coordinates
        intersections = [] #list of intersections, given in order as defined by desired_species_order

        self.debug_obj.debug_log(f"About to open file {output_filename}")

        with open(output_filename) as file:

            line_count = 0

            #iterate through to find line markers
            self.debug_obj.debug_log("Searching for line markers")
            for line in file:

                line_count += 1

                if marker_limiting_inequalities in line:
                    line_loc_limiting_inequalities = line_count

                if marker_intersecting_points in line:
                    line_loc_intersecting_points = line_count

                if marker_binary_points in line:
                    line_loc_binary_points = line_count
                    system_binary = True

                if marker_unstable in line:
                    line_loc_unstable = line_count

                if marker_competing_phase in line:
                    line_loc_competing_phase = line_count
                    system_stable = True

            if system_stable == True:
                line_locs_inequalities[0] = line_loc_limiting_inequalities + 2 #starts two lines down
                line_locs_inequalities[1] = line_loc_intersecting_points - 4 #starts four lines up
                if(system_binary == False):
                    line_locs_intersections[0] = line_loc_intersecting_points + 2
                    line_locs_intersections[1] = line_loc_competing_phase - 4
                    skip_to_inter_values = 2
                if(system_binary == True):
                    line_locs_intersections[0] = line_loc_binary_points + 5
                    line_locs_intersections[1] = line_loc_binary_points +7
                    skip_to_inter_values = 0

            if system_stable == False:
                line_locs_inequalities[0] = line_loc_limiting_inequalities + 2
                line_locs_inequalities[1] = line_loc_unstable - 4

        #Search for the inequalities
        with open(output_filename) as file:

            line_count = 0

            #print(f"line_locs_inequalities: {line_locs_inequalities}")
            #print(f"line_locs_intersections: {line_locs_intersections}")

            self.debug_obj.debug_log("Parsing inequality conditions")
            for line in file:

                line_count += 1

                if line_locs_inequalities[0] <= line_count <= line_locs_inequalities[1]:

                    self.debug_obj.debug_log("Found a inequality to check!")
                    line_info = line.split()
                    inequality_condition = self.parse_inequality_condition(line_info)
                    beef_inequalities.append(inequality_condition)

        #If stable, grab intersection points
        with open(output_filename) as file:

            if system_stable == True:

                self.debug_obj.debug_log("Parsing chemical potenital intersections")

                line_count = 0

                for line in file:

                    line_count += 1

                    #First line specifies the order of the provided info
                    if line_count == line_locs_intersections[0]:

                        line_info = line.split()

                        for item in line_info:
                            species_info = re.search(self.species_regex, item)
                            try:
                                species = species_info.group(1)
                                cplap_species_order.append(species)
                            except:
                                pass

                    #jump two lines to actual numbers
                    if line_locs_intersections[0] + skip_to_inter_values <= line_count <= line_locs_intersections[1]:

                        line_info = line.split()

                        new_coordinate = []

                        for item in line_info:
                            intersection_info = re.search(self.intersection_value_regex, item)
                            try:
                                intersection_val = float(intersection_info.group(1))
                                new_coordinate.append(intersection_val)
                            except:
                                pass

                        re_ordered_coordinate = self.reorder_CPLAP_list_of_values(cplap_species_order,new_coordinate)
                        intersections.append(re_ordered_coordinate)

        beef_entry['inequality_conditions'] = beef_inequalities
        beef_entry['intersections'] = intersections

        return beef_entry

    # Rather than directly use grid created by CPLAP, this code instead pulls the limiting inequalities
    # and stores a function of each limiting inequality. Later, it is possible to test points in chemical
    # potential space to see if they satisfy each function. If a point satisfies all functions of the
    # limiting inequalities, then it is within the stability envelope for that BEEF index.
    # Structured in this manner to allow the generation of aligned grid points across all BEEF-indeces
    # Parameter | beef_num : the beef index
    def generate_dependent_variable_calc_func(self,beef_num):

        self.debug_obj.debug_log("Entered generate_dependent_variable_calc_func")
        primary_idx = 0
        for i in range(0,len(self.phase_info)):
            item = self.phase_info[i]
            if item['type'] == 'primary':
                primary_idx = i

        beef_energy = self.phase_info[primary_idx]['chemical_potentials'][beef_num]

        independent_species = []
        independent_coefficients = []
        dependent_coefficient = 0

        for species in self.phase_info[primary_idx]['species_info']:

            independent_species.append(species)

            value = self.phase_info[primary_idx]['species_info'][species]
            independent_coefficients.append(value)

            if species == self.dependent_species:
                dependent_coefficient = value

        #This will strip out the dependent coefficient and leave the order as in the independent list
        independent_coefficients = self.reorder_CPLAP_list_of_values(independent_species,
                                                                     independent_coefficients)

        def temp(independent_coordinate_list):
            dependent_energy = beef_energy - np.dot(independent_coordinate_list,independent_coefficients)
            dependent_energy = dependent_energy/dependent_coefficient
            return dependent_energy

        return temp

    #the species_order is the desired order for the condition
    #This is because it is necessary to be consistent in which order any variables are provided
    #The CPLAP inequality list might not provide vairbales in the desired order
    def parse_inequality_condition(self,CPLAP_inequality_list):

        self.debug_obj.debug_log("Parsing a new inequality")
        self.debug_obj.debug_log(f"CPLAP_inequality_list: {CPLAP_inequality_list}")
        """ Single lines will be passed to this function. Take the lines create functions
           1.000 mu_Ba     1.000 mu_Sn  >    -9.7235
           0.666 mu_Ba    -0.333 mu_Sn  <    -1.2314
          -0.666 mu_Ba     0.333 mu_Sn  <     1.8372
           0.666 mu_Ba    -0.333 mu_Sn  <    -1.8863
          -0.333 mu_Ba     0.666 mu_Sn  <     1.0203
           1.000 mu_Ba                  <     0.0000
                           1.000 mu_Sn  <     0.0000
           1.000 mu_Ba                  >    -9.7235
                           1.000 mu_Sn  >    -9.7235
        """
        #test_A = ['1.000','mu_Ba','1.000','mu_Sn','>','-9.7235']

        LHS_max_idx = 0
        LHS_greater_than = False

        cplap_species_list = []
        cplap_coeffnt_list = []
        energy_value = 0

        for i in range(0,len(CPLAP_inequality_list)):
            item = CPLAP_inequality_list[i]

            if '>' in item or '<' in item:
                LHS_max_idx = i
                RHS_idx = i + 1
                energy_value = float(CPLAP_inequality_list[RHS_idx])

                if '>' in item:
                    LHS_greater_than = True

        #parse through the left hand side
        for i in range(0,LHS_max_idx):

            item_info = CPLAP_inequality_list[i]

            #check if species
            species_search_result = re.search(self.species_regex, item_info)

            try:
                species = species_search_result.group(1)
                cplap_species_list.append(species)
            except:
                pass

            #check if number
            number_search_result = re.search(self.intersection_value_regex, item_info)

            try:
                coefficient = float(number_search_result.group(1))
                cplap_coeffnt_list.append(coefficient)
            except:
                pass


        self.debug_obj.debug_log(cplap_species_list)
        self.debug_obj.debug_log(cplap_coeffnt_list)
        self.debug_obj.debug_log(LHS_greater_than)
        self.debug_obj.debug_log(energy_value)

        reordered_coefficient_list = self.reorder_CPLAP_list_of_values(cplap_species_list,cplap_coeffnt_list)
        inequality_function = self.generate_inequality_condition(reordered_coefficient_list,
                                                                 LHS_greater_than,
                                                                 energy_value)
        return inequality_function

    # coefficient_list - list of coefficients. This MUST be in the same order as the
    # LHS_greater_than - boolean. True if LHS (coord * coeff) is greater than RHS (energy)
    # energy_val - energy value of limit
    def generate_inequality_condition(self,coefficient_list,lhs_greater_than,energy_val):

        self.debug_obj.debug_log("Generating inequality condition")

        if lhs_greater_than == True:
            def temp(independent_coordinate_list):
                return np.dot(independent_coordinate_list,coefficient_list) > energy_val
            return temp
        else:
            def temp(independent_coordinate_list):
                return np.dot(independent_coordinate_list,coefficient_list) < energy_val
            return temp


    # desired_species_order - the desired order of species
    # cplap_species_order - the species order that CPLAP gives
    # cplap_coodinate - the CPLAP ccoridnate (not in the desired order)
    # Note! This will only give the coordinates of the independent variables!
    def reorder_CPLAP_list_of_values(self,cplap_species_order,cplap_coordinate):

        #first check that both species list have the same items

        cplap_species_length = len(cplap_species_order)
        cplap_coordinate_length = len(cplap_coordinate)

        if cplap_species_length != cplap_coordinate_length:
            self.debug_obj.debug_log(f"The CPLAP species and coordinate list are not the same length",type="WARNING")
            self.debug_obj.debug_log(f"cplap_species_order: {cplap_species_order}")
            self.debug_obj.debug_log(f"cplap_coordinate: {cplap_coordinate}")

        """
        if set(self.independent_species) != set(cplap_species_order):
            self.debug_obj.debug_log(f"independent_species: {self.independent_species}")
            self.debug_obj.debug_log(f"cplap_species_order: {cplap_species_order}")
        """

        length_of_lists = len(self.independent_species)

        reconfigured_coordinate = []

        for i in range(0,length_of_lists):

            try:
                desired_species = self.independent_species[i]

                index_of_desired_species_in_cplap_list = cplap_species_order.index(desired_species)

                value_of_desured_species_in_coordinate = cplap_coordinate[index_of_desired_species_in_cplap_list]

                reconfigured_coordinate.append(value_of_desured_species_in_coordinate)
            except:
                reconfigured_coordinate.append(0)

        return reconfigured_coordinate












#%%
