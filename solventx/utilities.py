#Utility functions for Gym
import json
import config


def read_config(file_name):
    """Load config json file and return dictionary."""
    
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









    