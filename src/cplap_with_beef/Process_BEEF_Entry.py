import copy
import importlib
import os
import subprocess
import itertools
import numpy as np

import re

from cplap_with_beef import Utility
importlib.reload(Utility)

from cplap_with_beef import Debug
importlib.reload(Debug)

# From CPLAP_Helper we get the following:
#  - a list of intersection of the boundary conditions
#  - the boundary conditions (multiple inequlities)
# But we don't yet have a grid of points in chemical potential space
# This class generates this grid of points for the given beef index
class Process_BEEF_Entry:

    def __init__(self,beef_idx,grid_spacing,debug_obj=False):

        #use the debug_obj that has been passed in, or make a new one
        self.debug_obj = False
        if debug_obj != False and isinstance(debug_obj,object):
            self.debug_obj = debug_obj
        else:
            self.debug_obj = Debug.Debug(debug_bool = True,restart = True)

        self.util_obj = Utility.Utility(debug_obj=self.debug_obj)

        self.beef_idx = beef_idx
        self.grid_spacing = grid_spacing
        self.number_of_independent_digits = 0
        self.number_of_dependent_digits = 0
        self.extra_digits_to_store_wrt_independent_variables = 2
        self.determine_number_of_digits()

        self.list_of_intersections = False
        self.list_of_boundary_conditions = False
        self.dependent_variable_calc = False
        self.number_of_independent_coordinates = 0
        self.valid_grid_points = False

    # Define the points where to start the search for our grid of points
    # However these provided points will be (almost guarenteed) to not be on
    # our desired grid spacing
    def set_list_of_intersections(self,list_of_intersections):
        self.list_of_intersections = list_of_intersections
        try:
            self.number_of_independent_coordinates = len(self.list_of_intersections[0])
            return True
        except:
            return False

    # Provide the boundary conditions for the search for grid points
    def set_list_of_inequalities(self,list_of_boundary_conditions):
        self.list_of_boundary_conditions = list_of_boundary_conditions

    def set_dependent_variable_calc(self,dependent_variable_calc):
        self.dependent_variable_calc = dependent_variable_calc

    def get_list_of_intersections(self):
        return self.list_of_intersections

    def get_list_of_inequalities(self):
        return self.list_of_boundary_conditions

    def get_dependent_variable_calc(self):
        return self.dependent_variable_calc

    # This returns the grid points in a tabular format. No labels are set.
    # All rows but last are the independent species, with ordering set by CPLAP_Helper
    # The last row is the dependent species
    def get_valid_grid_points(self):
        return self.valid_grid_points

    # Determine the number of digits that must be stored for a given grid spacing
    # This will define both the number of independent and dependent digits to store
    def determine_number_of_digits(self):
        temp_value = self.grid_spacing
        count = 0
        while(temp_value < 1):
            temp_value *= 10
            count += 1

        if temp_value == 1:
            self.number_of_independent_digits = count
            self.number_of_dependent_digits = count + self.extra_digits_to_store_wrt_independent_variables
            return count
        if temp_value > 1:
            self.debug_obj.debug_log("The grid spacing must be an integer order of magnitude from 1 (e.g. 0.1,0.01)")
            return False

    # Given an exact value, round the value w/r/t the grid spacing
    # e.g. if grid spacing = 0.01, 15.451 becomes 15.45
    def floor_wrt_grid_spacing(self,exact_value):

        rounded_val = round(exact_value,self.number_of_independent_digits)

        if(rounded_val>exact_value):
            return (rounded_val-self.grid_spacing)
        else:
            return rounded_val

    # This will generate a list of intersection neighbors
    # This DOES NOT analyze if they are valid. It merely generates the list of points
    # If grid_buffers = 0, will only generate a single point (the closest one when rounded)
    # If grid buffers > 0, will create a number of grid points around the provided one
    def generate_intersection_neighbors(self,intersection_point,grid_buffers = 1):

        if grid_buffers == 0:

            closest_point = []

            for exact_value in intersection_point:
                floor = self.floor_wrt_grid_spacing(exact_value)
                ceiling = floor + self.grid_spacing

                if abs(exact_value-ceiling) < abs(exact_value-floor):
                    closest_point.append(ceiling)
                else:
                    closest_point.append(floor)

            return [closest_point]

        else:
            #generate a list of rounded intersection values w.r.t the grid spacing
            minimum_values = []
            for exact_value in intersection_point:
                min_val = self.floor_wrt_grid_spacing(exact_value)
                min_val -= (self.grid_spacing*grid_buffers)
                minimum_values.append(min_val)
            #print(f"minimum_values: {minimum_values}")

            num_dim = len(intersection_point) #the number of dimensions
            range_of_points = 2*grid_buffers+1
            ndim_grid = [range(range_of_points)]*num_dim

            addition_list = []
            for coord in itertools.product(*ndim_grid, repeat=1):
                new_addition = np.array(coord)*self.grid_spacing
                addition_list.append(new_addition)

            final_list = []
            for addition_item in addition_list:
                final_val = np.add(minimum_values,addition_item)
                final_list.append(final_val)

            return final_list

    # Determine if a particular point has already been added to the valid grid points
    # the dependent value can be provided or not (it will just be ignored)
    def check_if_grid_point_added(self,grid_point):
        num_ind = self.number_of_independent_coordinates
        matches = (self.valid_grid_points[:,0:num_ind] == grid_point[0:num_ind]).all(axis=1).nonzero()[0]
        if len(matches) > 0:
            return True
        else:
            return False

    # From the provided point, find the closest points that lies along the line as ...
    # specified by the index_of_interest
    # Will return to lists
    # list 1: booleans indicating if smaller and larger values have been found (e.g. [True,True]
    # list 2: floats indicating the closest smaller and larger value found along direction ...
    # for example if original_point = [3,4,5], index_of_interest = 2, then list_2 could = [3,7]
    # This would indicate that two valid neighboring points are [3,4,3] and [3,4,7]
    def find_closest_neighbors_along_direction(self,original_point,index_of_interest):

        original_point = original_point.round(decimals=self.number_of_independent_digits)

        # The starting list of indeces will include every row
        # The row index will be continually narrowed
        row_index_of_all_matches = range(0,np.shape(self.valid_grid_points)[0])
        number_of_dimensions = len(original_point)

        max_smaller_value = 0
        min_larger_value = 0

        smaller_value_found = True
        larger_value_found = True

        #find the points where all other columns match
        for i in range(0,number_of_dimensions):
            if i != index_of_interest:
                matching_rows = np.where(self.valid_grid_points[:,i].round(decimals=self.number_of_independent_digits) == original_point[i])[0]
                row_index_of_all_matches = set(row_index_of_all_matches).intersection(matching_rows)
                row_index_of_all_matches = list(row_index_of_all_matches)


        #print(f"row_index_of_all_matches: {row_index_of_all_matches}")

        row_indexes_with_larger_values = [i for i in row_index_of_all_matches
                                          if self.valid_grid_points[i,index_of_interest]
                                          > original_point[index_of_interest]]

        #print(f"row_indexes_with_larger_values: {row_indexes_with_larger_values}")

        try:
            min_larger_value = np.min(self.valid_grid_points[row_indexes_with_larger_values,index_of_interest])
        except:
            min_larger_value = original_point[index_of_interest]
            larger_value_found = False

        row_indexes_with_smaller_values = [i for i in row_index_of_all_matches
                                           if self.valid_grid_points[i,index_of_interest]
                                           < original_point[index_of_interest]]

        #print(f"row_indexes_with_smaller_values: {row_indexes_with_smaller_values}")

        try:
            max_smaller_value = np.max(self.valid_grid_points[row_indexes_with_smaller_values,index_of_interest])
        except:
            max_smaller_value = original_point[index_of_interest]
            smaller_value_found = False


        #print(f"original_point: {original_point}")
        #print(f"index_of_interest: {index_of_interest}")
        #print(f"[smaller_value_found,larger_value_found]: {[smaller_value_found,larger_value_found]}")
        #print(f"[max_smaller_value,min_larger_value]: {[max_smaller_value,min_larger_value]}")

        return [smaller_value_found,larger_value_found],[max_smaller_value,min_larger_value]

    # The dedicated function to add points to the validated list of points
    # Param: force_to_grid - if true, will find points near the provided points that are on the desired grid
    # Param: check_if_within_boundaries - check if within boundaries! False only if sure points are valid
    # Param: check_for_duplicates - remove duplicates if found
    # The params are all True by default, and have been added to make the function more efficient in case
    # ... the upstream logic makes the provided checks unnecessary
    def add_to_valid_grid_points(self,list_of_provided_points,
                                 force_to_grid = True,
                                 check_if_within_boundaries = True,
                                 check_for_duplicates = True,
                                 grid_buffers = 1):

        list_of_points_on_grid = []
        validated_points = []
        #print(f"Length of list_of_provided_points: {len(list_of_provided_points)}")

        initial_number_of_valid_points = self.count_valid_points()

        #make sure only independent variables are passed in
        provided_points_independent_variables = []
        for coordinate in list_of_provided_points:
            provided_points_independent_variables.append(coordinate[0:self.number_of_independent_coordinates])

        #Find points near the ungridded points that are on the desired grid
        if force_to_grid == True:
            for coordinate in provided_points_independent_variables:
                new_points_on_grid = self.generate_intersection_neighbors(coordinate,grid_buffers = grid_buffers)
                list_of_points_on_grid.extend(new_points_on_grid)
        else:
            list_of_points_on_grid = provided_points_independent_variables

        #Check if the points satisfy the known boundaries
        if check_if_within_boundaries == True:

            for coordinate in list_of_points_on_grid:

                all_pass = True
                for condition in self.list_of_boundary_conditions:

                    try:
                        if condition(coordinate) == False:
                            all_pass = False
                            break
                    except:
                        print(f"Exception! coordinate: {coordinate}")

                if all_pass == True:
                    validated_points.append(coordinate)
        else:
            validated_points = list_of_points_on_grid

        #Add the points to the list of valid grid points
        for coordinate in validated_points:

            dependent_value = self.dependent_variable_calc(coordinate)

            # The dependent value must be less than zero to be chemical stable
            # Must add this check because even though independent points are valid, the dependent
            # value could still be positive
            if dependent_value < 0:
                info_to_add = np.append(coordinate,dependent_value)
                self.valid_grid_points = np.append(self.valid_grid_points,[info_to_add],axis=0)

        #At end, remove any duplicates. This is very costly
        if check_for_duplicates == True:

            #round each row. Need to do this because of floating point numbers
            for i in range(0,len(self.valid_grid_points)):
                self.valid_grid_points[i] = self.valid_grid_points[i].round(decimals=self.number_of_dependent_digits)

            #find the unique rows
            unique_rows = np.unique(self.valid_grid_points[:,0:self.number_of_independent_coordinates],
                                    axis=0,
                                    return_index=True)
            unique_rows = unique_rows[1]

            self.valid_grid_points = [self.valid_grid_points[i,:] for i in range(0,len(self.valid_grid_points)) if i in unique_rows]
            self.valid_grid_points = np.array(self.valid_grid_points)

        final_number_of_valid_points = self.count_valid_points()
        number_of_points_added = final_number_of_valid_points - initial_number_of_valid_points
        #print(f"number_of_points_added: {number_of_points_added}")

        return number_of_points_added

        #print(f"valid_grid_points:\n{self.valid_grid_points}")

    # Return the number of valid grid points
    def count_valid_points(self):
        return (self.valid_grid_points.shape)[0]

    # From the list of intersections, find all valid grid points
    # Will continue to densify until no more points are found
    def organize_grid_densify(self):

        #print(f"beef index: {self.beef_idx}\n")
        grid_shape = (0,self.number_of_independent_coordinates + 1)
        self.valid_grid_points = np.zeros(grid_shape)

        # the intersections provided by CPLAP are not enough
        # Generate points BETWEEN the intersections to find all points ...
        # along the vertices
        self.interpolate_between_valid_points(self.list_of_intersections)
        number_of_starting_valid_points = self.count_valid_points()

        # In this case, the  stability window exists, but is super narrow
        # as a result, no points for the provided grid spacing sit on the grid
        if number_of_starting_valid_points == 0:
            return False

        keep_running = True

        while(keep_running == True):
            old_count = self.count_valid_points()
            #print(f"Current Number of points: {old_count}")

            self._densify_the_grid_one_run()

            new_count = self.count_valid_points()
            #print(f"Current Number of points: {new_count}")

            if(old_count == new_count):
                #print(f"All points found!")
                keep_running = False

        self.valid_grid_points = self.valid_grid_points[np.lexsort(np.transpose(self.valid_grid_points)[::-1])]

        self.util_obj.change_directory_auto("grid_exact")
        exact_result_filename = f"beef_{self.beef_idx}_exact.csv"
        np.savetxt(exact_result_filename, self.valid_grid_points, delimiter=",")


    # From the list of provided points, find points between every single point that are valid
    # Most useful for the starting list of intersections. This finds valid points on the edges of the
    # chemical potential space, which meaningfully speeds up the process
    def interpolate_between_valid_points(self,starting_list):

        #print("Entered interpolate_between_valid_points")

        new_points = []

        secondary_list = copy.deepcopy(starting_list)

        for primary_point in starting_list:

            # Make a new temporary list without any primary points
            # The secondary list will continue to get smaller and smaller
            secondary_list.remove(primary_point)

            #interpolate between the primary point and the secondary points
            for secondary_point in secondary_list:

                diff = np.subtract(secondary_point,primary_point)

                #determine the proper number of increments needed
                diff_length = np.linalg.norm(diff)
                number_of_steps = int(diff_length/self.grid_spacing) + 1

                #generate exact points between primary and secondary
                #these points are NOT on the grid
                for step in range(0,number_of_steps):
                    new_point = (step/number_of_steps)*diff + primary_point
                    new_points.append(new_point)

        self.add_to_valid_grid_points(new_points,force_to_grid = True,
                                                 check_if_within_boundaries = True,
                                                 check_for_duplicates = True,
                                                 grid_buffers=1)


    #For every valid point, search in every direction for more valid points
    #Called "one_run" since this function must be called multiple times
    #Otherwise some points might be missed
    def _densify_the_grid_one_run(self):

        self.debug_obj.debug_log(f"Densifying the grid!")
        old_coordinate_list = copy.deepcopy(self.valid_grid_points)

        #Adding these to help keep track of progress
        count_total = len(old_coordinate_list)
        count_one_hundredth = int(count_total/100)+1

        number_of_old_points = len(old_coordinate_list)
        number_of_dimensions = len(old_coordinate_list[0])

        count_nearest_neighbor_found = 0
        count_neighbor_found = 0
        count_pushing_boundary = 0

        for row_idx in range(0,number_of_old_points):

            if row_idx%count_one_hundredth==0:
                multiple = int(row_idx/count_one_hundredth)
                #print(f"{multiple}% done",end="\r")

            primary_point = old_coordinate_list[row_idx]

            for species_idx in range(0,number_of_dimensions):

                if species_idx == self.number_of_independent_coordinates-1:
                    break

                #check both directions
                for i in range(0,2):

                    multiplier = 1
                    if i == 1:
                        multiplier = -1

                    #If nearest neighbor already exists, break from loop and look in another direction
                    test_point = primary_point.copy()
                    test_point[species_idx] += self.grid_spacing*multiplier
                    if self.check_if_grid_point_added(test_point) == True:
                        count_nearest_neighbor_found += 1
                        break

                    search_result,closest_neighbors_along_direction = self.find_closest_neighbors_along_direction(primary_point,species_idx)
                    new_point = primary_point.copy()
                    new_points = []

                    #decide if it is safe to fill valid points without having to check boundaries)
                    fill_with_abandon = True
                    if multiplier == 1 and search_result[1] is False:
                        fill_with_abandon = False
                    if multiplier == -1 and search_result[0] is False:
                        fill_with_abandon = False

                    if fill_with_abandon == True:
                        count_neighbor_found += 1

                        while(closest_neighbors_along_direction[0] <= new_point[species_idx] <= closest_neighbors_along_direction[1]):
                            new_points.append(new_point.copy())
                            new_point[species_idx] += self.grid_spacing*multiplier

                        #Add without checking since guaranteed to be valid
                        #Already on grid, don't force to grid
                        #Don't worry about duplicates right now
                        self.add_to_valid_grid_points(new_points,
                                                      force_to_grid = False,
                                                      check_if_within_boundaries = False,
                                                      check_for_duplicates = False,
                                                      grid_buffers=0)

                        #break here because no need to check further (since a point was found)
                        break

                    # If logic gets here, then no neighbors found along direction
                    # Check that every new point is valid. Stop when no points added
                    count_pushing_boundary += 1
                    kick_start = True
                    points_added = 0
                    while(kick_start or points_added > 0):
                        kick_start = False
                        new_point[species_idx] += self.grid_spacing*multiplier
                        points_added = self.add_to_valid_grid_points([new_point],
                                                                     force_to_grid = True,
                                                                     check_if_within_boundaries = True,
                                                                     check_for_duplicates = False,
                                                                     grid_buffers=0)

        # At end of process check for duplicates
        # Currently quite time consuming
        self.add_to_valid_grid_points([],force_to_grid = False,
                                         check_if_within_boundaries = False,
                                         check_for_duplicates = True)

        #print(f"count_nearest_neighbor_found: {count_nearest_neighbor_found}")
        #print(f"count_neighbor_found: {count_neighbor_found}")
        #print(f"count_pushing_boundary: {count_pushing_boundary}")


#%%
