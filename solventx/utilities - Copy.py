#Utility functions for Gym
import json
import os
import numpy as np
import math
from solventx import config
from solventx import templates

def read_config(file_name):
    """Load config json file and return dictionary."""
    
    print ('ola ola olaaaaaaa .......\n', file_name, 'ola ola olaaaaaaa ......\n')
    
    with open(file_name, "r") as config_file:
        print(f'Reading configuration file:{config_file.name}')
        confDict = json.load(config_file)
        
    return confDict



def get_process(obj):
    """Get products."""
    
    
    input_components = obj.confDict['modules']['input']
    strip_components = obj.confDict['modules']['output']['strip']


    n_components = len(input_components)
    config_key = ''
    
    print(f'Looping through following modules config:{list(config.valid_processes.keys())}')
    for key,config_dict in config.valid_processes.items():
        
       
        if set(input_components) == set(config_dict['input']):
            if set(strip_components) ==  set(config_dict['strip']):
                
                config_key = key
    
    if config_key:
        print(f'Found the following process config:{key}')            
    else:
        raise ValueError(f'No valid configuration found for input:{input_components},strip:{strip_components}!')            
   
    modules = config.valid_processes[config_key]['modules']
    x = []
   
    print(f'Process config {config_key}:Input:{input_components},Number of modules:{len(modules)}')
    print('Modules info:')
    for key,module in modules.items():
        x.extend(module['x'])
        print(f'Module {key}:{module["strip_group"]}')
    print(f'x0:{x}')

     
    return x, modules, n_components




def generate(confDict, N):
    
    home_dir = confDict['solventxHome']
    reeComps        = confDict['compositions']
    modulesData     = confDict['modules']      
    ree             = modulesData["input"]

    upper           = [reeComps[i]['upper'] for i in ree]
    lower           = [reeComps[i]['lower'] for i in ree]
    
    ree_mass_dict   = dict()
    
    for index in range(N): # N represents number of ree composition cases
        ree_mass_dict[str(index)] = {}
        ree_mass = [np.random.uniform(i,j) for i,j in zip(lower, upper)]  # select composition from bounds
        for item,jtem in zip(ree, ree_mass):
            ree_mass_dict[str(index)][item] = jtem
    
    with open(os.path.join(home_dir,'.\\solventx\\data\\json\\cases.json'),'w') as json_file:
        json.dump(ree_mass_dict, json_file)   
        
    return ree_mass_dict



def get_env_config_dict(config_file):
    """Read config file create confi dict."""
    
    assert 'json' in config_file, 'Config file must be a json file!'
    config_keys = templates.config_keys
    
    design_config = read_config(config_file)
    
    config_dict = {}
    for key in config_keys:
        if key in design_config.keys():
            config_dict.update({key:design_config[key]})
        else:
            raise ValueError(f'{key} not found in config JSON file!')
   
#    variable_config= config_dict['variable_config']
#    upper = [variable_config[jtem]['upper'] for jtem in variable_config.keys() ]
#    lower = [variable_config[jtem]['lower'] for jtem in variable_config.keys() ]
#    vtype = [variable_config[jtem]['type'] for jtem in variable_config.keys() ]
#
#    config_dict.update({'lower':lower})   
#    config_dict.update({'upper':upper})   
#    config_dict.update({'type':vtype})   
    
    return config_dict 

    