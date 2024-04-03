
import os
import importlib
import shutil

from cplap_with_beef import Debug
importlib.reload(Debug)

class Utility(object):

    def __init__(self,debug_obj=False):

        #use the debug_obj that has been passed in, or make a new one
        self.debug_obj = False
        if debug_obj != False and isinstance(debug_obj,object):
            self.debug_obj = debug_obj
        else:
            self.debug_obj = Debug.Debug(debug_bool = True,restart = True)

    def checkEnvironmentalVars(self,env_var):

        env_var_key = ""
        try:
            env_var_key = os.environ[env_var]
            self.debug_obj.debug_log(f"Environment variable {env_var}': {env_var_key}")
            return env_var_key
        except KeyError:
            message = f"The environment variable {env_var} not set. Set this before running this program"
            self.debug_obj.debug_log(message,type="Error")
            return False

    #filename = name of filepath/filename
    def checkIfFileExists(self,filename):
        file_found = False

        #searcht the local directory
        item_list = os.listdir()
        for item in item_list:
            if filename in item:
                file_found = True

        if file_found == False:
            #search the exact filename (might be a full path)
            try:
                file_found = os.path.isfile(filename)
            except:
                file_found == False

        return file_found

        #filename = name of filepath/filename
    def createNewDirectory(self,dir_path):
        result = False
        try:
            os.mkdir(dir_path)
            result = True
        except:
            self.debug_obj.debug_log(f"Current working directory: {os.getcwd()}")
            self.debug_obj.debug_log(f"Failed to create path: {dir_path}")
        return result

    #filename = name of filepath/filename
    def changeDirectory(self,dir_path):
        result = False
        try:
            os.chdir(dir_path)
            result = True
        except:
            self.debug_obj.debug_log(f"Current working directory: {os.getcwd()}")
            self.debug_obj.debug_log(f"Unable to change location to '{dir_path}'",type="WARNING")
        return result

    #If precise is False, if term_to_delete is in the file, the file is deleted
    #If precise is True, only the exact name will result in a deletion
    def deleteFiles(self,term_to_delete,precise = True):

        self.debug_obj.debug_log(f"Attempting to delete files.")

        delete_successful = False

        try:
            item_list = os.listdir()
            for item in item_list:

                if precise == True:
                    if term_to_delete == item:
                        os.remove(item)
                        delete_successful = True
                        self.debug_obj.debug_log(f"Deleted file: {item}")

                if precise == False:
                    if term_to_delete in item:
                        os.remove(item)
                        delete_successful = True
                        self.debug_obj.debug_log(f"Deleted file: {item}")

            if delete_successful == False:

                file_found = os.path.isfile(term_to_delete)
                if file_found == True:
                    os.remove(term_to_delete)
                    self.debug_obj.debug_log(f"Deleted file: {term_to_delete}")

        except:
            self.debug_obj.debug_log(f"Unable to delete any files.",type="WARNING")

    def copyFiles(self,src_path,dst_path,delete_old = False):

        shutil.copyfile(src_path, dst_path)

        copy_successful = self.checkIfFileExists(dst_path)

        if copy_successful:

            if delete_old == True:
                self.deleteFiles(src_path,precise=True)

                src_still_exists = self.checkIfFileExists(src_path)

                if src_still_exists == False:
                    delete_successful = True

        if copy_successful == True and delete_old == False:
            return True
        if copy_successful == True and delete_successful == True:
            return True
        else:
            if copy_successful == False:
                self.debug_obj.debug_log("File copy failed",type="WARNING")
            if delete_successful == False:
                self.debug_obj.debug_log("File delete failed",type="WARNING")
            return False

    def change_directory_auto(self,name_of_desired_folder,upstream_limit = 5):

        address_dict = self.fetch_directory_dict(upstream_limit)

        #Check that address_dict is a non-empty dict
        if address_dict:

            for loc in address_dict:

                if(name_of_desired_folder == loc):
                    self.debug_obj.debug_log(f"attempting to change folder to {address_dict[loc]}")
                    self.changeDirectory(address_dict[loc])
                    return True

        return False

    def fetch_directory_dict(self,upstream_limit = 5):

        filename_of_paths = ".file_paths.txt"
        file_found = False

        #Get current folder path
        original_path = os.getcwd()

        for i in range(0,upstream_limit):

            file_found = self.checkIfFileExists(filename_of_paths)

            if file_found:
                break

            os.chdir("..")

        #if found not found rteurn false
        if file_found == False:
            self.changeDirectory(original_path)
            return False

        address_dict = {}

        with open(filename_of_paths) as file:

            for line in file:

                line_info = line.split(' ')

                if line_info[1] == "=":
                    folder_name = line_info[0].strip()
                    folder_path = line_info[2].strip()
                    address_dict[folder_name] = folder_path


        #Go back to orginal folder
        self.changeDirectory(original_path)

        if address_dict:
            return address_dict
        else:
            return False