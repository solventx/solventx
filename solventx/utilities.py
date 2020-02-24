#Utility functions for Gym
import json

def read_config(file_name):
    """Load config json file and return dictionary."""
    
    with open(file_name, "r") as config_file:
        print(f'Reading configuration file:{config_file.name}')
        confDict = json.load(config_file)
        
    return confDict