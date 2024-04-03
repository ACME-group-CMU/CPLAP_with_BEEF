import copy
import importlib

import numpy as np
from numpy import genfromtxt

from cplap_with_beef import Utility
importlib.reload(Utility)

from cplap_with_beef import Debug
importlib.reload(Debug)

class Data_Analyzer:
    """
    param: debug_obj
    param: desired_species_order - the last item will become the dependent variable
    """
    def __init__(self,debug_obj=False):

        #use the debug_obj that has been passed in, or make a new one
        self.debug_obj = False
        if debug_obj != False and isinstance(debug_obj,object):
            self.debug_obj = debug_obj
        else:
            self.debug_obj = Debug.Debug(debug_bool = True,restart = True)

        self.util_obj = Utility.Utility(debug_obj=self.debug_obj)
        self.beef_index_max = 2000

        self.species_order = False
        self.combined_beef_data = False

    def set_species_order(self,species_order):
        self.species_order = species_order

    def get_species_order(self):
        return self.species_order

    # This returns a nicely packaged list of viable grid points
    # returns list of dictionaries, where each dictionary is a single stable point
    # e.g. [{'Ga': 1, 'N': 2},{'Ga': 2, 'N': 3}]
    def get_packaged_grid_points(self,start_index,stop_index,find_centroid = False,get_total_cp = False):

        all_beefs = {}

        for beef_idx in range(start_index,stop_index):

            valid_grid_points = self._fetch_one_chemical_potenital_result(beef_idx=beef_idx,
                                                                          find_centroid = find_centroid,
                                                                          get_total_cp = get_total_cp)

            all_beefs[beef_idx] = valid_grid_points

        return all_beefs

    # This will select a single chemical potential and return it
    # This function is not meant to be called by the user, but is an internal support function
    # parameter: beef_idx = the index of the result desired
    #returns a list of dictionaries, where each dictionary is one possible set of possible chemical potentials
    #the list format is chosen because nothing meaningfuly about the index of the set
    # e.g. [{'Ga': 1, 'N': 2},{'Ga': 2, 'N': 3}]
    # Note! Sometime a certain beef_idx will not have any valid chemical potenitals points ...
    # .. in this case, the list will be empty
    def _fetch_one_chemical_potenital_result(self,beef_idx = None,find_centroid = False,get_total_cp = False):
        #Get the order of elements from CPLAP_Helper
        #CPLAP_Helper sets the order and is the source of truth

        all_points = []

        file_found,beef_data = self.read_in_beef_file(beef_idx = beef_idx,total_energy=get_total_cp)

        if file_found == False:
            print("File not found for beef {beef_idx}")
            return all_points

        #Get centroid estimate
        centroid_est = self.estimate_centroid_loc(beef_data)
        min_centroid_dist = 100000 #hacky, make better
        centroid_index = -1

        num_rows = len(beef_data)
        for row_index in range(0,num_rows):

            new_entry = {}
            for col_index in range(0,len(self.species_order)):
                new_entry[self.species_order[col_index]] = beef_data[row_index,col_index]

                if(find_centroid):
                    distance_to_centroid = np.linalg.norm(np.array(beef_data[row_index,:]) - np.array(centroid_est))
                    print()
                    if(distance_to_centroid < min_centroid_dist):
                        min_centroid_dist = distance_to_centroid
                        centroid_index = row_index
                    new_entry['centroid'] = False #adding to all locs

            all_points.append(new_entry)

        if(find_centroid):
            all_points[centroid_index]['centroid'] = True

        return all_points

    def estimate_centroid_loc(self,beef_data):
        centroid_estimate = []

        #Iterate through
        for col_index in range(0,len(self.species_order)):
            coord_mean = np.mean(beef_data[:,col_index])
            centroid_estimate.append(coord_mean)

        print(f"centroid_estimate: {centroid_estimate}")
        return centroid_estimate

    def agregrate_beef_info(self):

        self.util_obj.change_directory_auto("grid_exact")

        first_file = True

        self.beef_index_max = 2001

        for beef_idx in range(0,self.beef_index_max):

            print(f"beef_idx: {beef_idx}",end="\r")

            file_found,beef_data = self.read_in_beef_file(beef_idx = beef_idx)

            if file_found == False:
                continue

            if first_file == True:
                column_count = beef_data.shape[1]

                #col indexes for the files being read in
                self.single_file_ind_col_indexes = [0,column_count-1]
                self.single_file_dep_col_index = column_count-1
                self.number_of_independent_variables = column_count-1

                #col indexes for the aggregate file
                #creating two blank column
                self.agregate_file_ind_col_indexes = self.single_file_ind_col_indexes
                self.agregate_file_dep_col_index = [column_count,column_count+2000]
                self.agregate_file_count_col_index = column_count+2001
                self.agregate_file_columns = column_count + 2002

                self.combined_beef_data = np.empty(shape=[0,self.agregate_file_columns])

                first_file = False

            for row in beef_data:
                grid_point = row[0:self.single_file_ind_col_indexes[1]]
                #print(grid_point)
                row_index = self.find_row_where_grid_point_added(grid_point)
                #print(row_index)
                if row_index:
                    self.combined_beef_data[row_index,(self.agregate_file_dep_col_index[0] + beef_idx)] = row[self.single_file_dep_col_index]
                    self.combined_beef_data[row_index,self.agregate_file_count_col_index] += 1
                elif row_index == False:
                    temp = np.zeros(shape=[1,self.agregate_file_columns])
                    temp[0,0:self.agregate_file_ind_col_indexes[1]] = grid_point
                    temp[0,self.agregate_file_dep_col_index[0] + beef_idx] = row[self.single_file_dep_col_index]
                    temp[0,self.agregate_file_count_col_index] = 1
                    self.combined_beef_data = np.concatenate([self.combined_beef_data,temp])

        aggregated_filename = "combined_beef_data.csv"
        np.savetxt(aggregated_filename, self.combined_beef_data, delimiter=",")

    # Determine if a particular point has already been added to the agregate list
    # the dependent value can be provided or not (it will just be ignored)
    # This makes it possible to place points w/ the same ind. variables in the same row
    def find_row_where_grid_point_added(self,grid_point):

        num_ind = self.single_file_ind_col_indexes[1]
        matches = (self.combined_beef_data[:,0:num_ind] == grid_point[0:num_ind]).all(axis=1).nonzero()[0]
        if len(matches) > 0:
            return matches[0]
        else:
            return False

    #columns should be a list with two values. First is the x column, Second is the y column
    #column_names is list of strings
    # This will currently work with a system with three chemical potentials
    # If more than three are used, this will break because multiple rows will have identical 'x' and 'y' values
    # When you implement for more than three, maybe make a new data structure that combines these rows
    def generate_2d_heat_map(self,columns,column_names,limit):

        import numpy as np
        import matplotlib.pyplot as plt
        import scipy.interpolate as interp
        import matplotlib

        start = 3

        # select columns
        x_raw = self.combined_beef_data [:, columns[0]] # column 0
        y_raw = self.combined_beef_data [:, columns[1]] # column 1

        x = []
        y = []
        z = []

        for i in range(0,len(x_raw)):
            intensity = (self.combined_beef_data[i,start:start+limit] != 0).sum()
            if intensity > 0:
                x.append(x_raw[i])
                y.append(y_raw[i])
                z.append(intensity)

        # create a grid of coordinates
        x_grid, y_grid = np.meshgrid (np.linspace (min (x), max (x), 1000), np.linspace (max (y), min (y), 1000))

        # interpolate the z values at the grid points
        z_grid = interp.griddata ((x, y), z, (x_grid, y_grid), method='linear')

        # plot the image
        plt.imshow (z_grid, cmap='plasma', vmin=1, vmax=max(z),extent=[min(x),max(x),min(y),max(y)],interpolation='none')

        plt.xlabel(f"{column_names[0]}")
        plt.ylabel(f"{column_names[1]}")
        plt.title("Heat Map of BEEF-Derived \n Chemical Potential Envelopes")

        plt.colorbar (label="Number of Stable BEEF-Derived \n Chemical Potential Envelopes",shrink=0.7)


    def generate_standard_chempot_map(self,columns,column_names,beef_idx):

        import numpy as np
        import matplotlib.pyplot as plt
        import scipy.interpolate as interp
        import matplotlib

        start = 3

        # select columns
        x_raw = self.combined_beef_data[:, columns[0]] # column 0
        y_raw = self.combined_beef_data[:, columns[1]] # column 1
        z_raw = self.combined_beef_data[:, (start+beef_idx)] # column 2

        x = []
        y = []
        z = []

        print("printing all points")
        print(f"{len(x_raw)}")

        for i in range(0,len(x_raw)):

            if z_raw[i] != 0:
                x.append(x_raw[i])
                y.append(y_raw[i])
                z.append(z_raw[i])
                #print(f"{x_raw[i]},{y_raw[i]},{intensity}")

        # create a grid of coordinates
        x_grid, y_grid = np.meshgrid (np.linspace (min (x), max (x), 1000), np.linspace (max (y), min (y), 1000))

        # interpolate the z values at the grid points
        z_grid = interp.griddata ((x, y), z, (x_grid, y_grid), method='linear')

        # plot the image
        plt.imshow (z_grid, cmap='plasma', vmin=min(z), vmax=max(z),extent=[min(x),max(x),min(y),max(y)],interpolation='nearest')

        plt.xlabel(r"$\Delta \mu_{Pb}$ [eV]",fontsize=12)
        plt.ylabel(r"$\Delta \mu_{Br}$ [eV]",fontsize=12)
        plt.title("$\Delta \mu_{Cs}$ [eV] Envelope for BEEF #0",fontsize=15)


        #plt.colorbar (label=r"$\Delta \mu_{Cs}$ [eV]",shrink=0.7,weight="bold")
        plt.colorbar(shrink=0.7).set_label(label=r"$\Delta \mu_{Cs}$ [eV]",size=12)
        plt.show ()

    #by default, this will find the change in chemical potential (the CPLAP result)
    #if total_energy = True, the complete chemical potenial will be returned
    def read_in_beef_file(self,beef_idx,total_energy = False):

        generic_filename = ""

        if total_energy == False:
            generic_filename = "beef_#_exact.csv"
            self.util_obj.change_directory_auto(name_of_desired_folder="grid_exact")
        else:
            generic_filename = "beef_#_grid_total.csv"
            self.util_obj.change_directory_auto(name_of_desired_folder="grid_total")

        precise_filename = generic_filename.replace("#",str(beef_idx))
        file_found = self.util_obj.checkIfFileExists(precise_filename)

        if file_found:
            beef_data = genfromtxt(precise_filename, delimiter=',')
            return file_found,beef_data
        else:
            return False,False

    def read_in_agregrate_file(self,force_read_in = False):

        filename = "combined_beef_data.csv"

        if self.combined_beef_data is False or force_read_in == True:
            self.util_obj.change_directory_auto("grid_exact")

            file_exists = self.util_obj.checkIfFileExists(filename)

            if file_exists:
                self.combined_beef_data = genfromtxt(filename, delimiter=',')

        else:
            print("self.combined_beef_data already contains data.")

#%%
