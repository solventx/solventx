# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 12:13:08 2020

@author: ciloeje
"""
import config

def get_process(obj):
    """Get products."""
    
    #print(self.process_config,self.sx_design.confDict)
    
#    input_components = input_list #self.sx_design.confDict['modules']['input']
#    strip_components = strip_list #self.sx_design.confDict['modules']['output']['strip']
#    #extraction_components = self.sx_design.confDict['modules']['output']['extraction']
    
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
        raise ValueError(f'No configuration found for input:{input_components},strip:{strip_components}!')            
   
    modules = config.valid_processes[config_key]['modules']
    x = []
   
    print(f'Process config {config_key}:Input:{input_components},Number of modules:{len(modules)}')
    print('Modules info:')
    for key,module in modules.items():
        x.extend(module['x'])
        print(f'Module {key}:{module["strip_group"]}')
    print(f'x0:{x}')
    
    return x,n_components