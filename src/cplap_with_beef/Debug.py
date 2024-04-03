import sys
from datetime import datetime
import os
import importlib

from cplap_with_beef import Utility
importlib.reload(Utility)


class Debug(object):
    
    def __init__(self,debug_bool,restart = False):

        #self.util_obj = Utility.Utility()

        self.debug = debug_bool
        self.errors_exist = False
        self.debug_data = ""

        self.debug_file = "debug.txt"
        #directory_dict = self.util_obj.fetch_directory_dict()
        #main_path = directory_dict['main']
        self.debug_file = f"{os.getcwd()}/chemical_potential_info/{self.debug_file}"

        now = datetime.now()

        if restart == True:

            #clear old file
            with open(self.debug_file, 'w') as f:
                f.write(f"##################################################################")
                f.write(f"\nLog of Defect Tracking Run")
                f.write(f"\nDate and time: {now}")
                f.write(f"\n##################################################################\n")
                f.close()

    def debug_log(self,message,type = None):

        whole_message = ""

        if type == None or type == "LOG":
            whole_message = f"\nLOG | {message}"

        if type == "WARNING" or type == "ERROR":
            whole_message = f"\n{type} | {message}\n"

        if type == "ERROR":
            self.errors_exist = True

        if self.debug == True:
            print(whole_message)

        #data = data + f"\n{message}"

        with open(self.debug_file, 'a') as f:
            f.write(whole_message)

    def kill_if_errors(self):

        if self.errors_exist:
            self.debug_log("Attempting to kill run.")
            sys.exit(1)
                
                
            