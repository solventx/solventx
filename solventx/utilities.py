#Utility functions for Gym
import json
import numpy as np
import math
import rbfopt
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
    
    with open('.\\solventx\\data\\json\\cases.json','w') as json_file:
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
   
    
    return config_dict 


def optimize(myb, iters=100, strategy='lhd_corr',target=1, restarts=3):
    
    # Specify optimization settings
    settings = rbfopt.RbfoptSettings(max_evaluations=iters, init_strategy=strategy, target_objval=target, max_noisy_restarts=restarts) #'lhd_corr')#, minlp_solver_path=minlp_solver, nlp_solver_path=nlp_solver) 

    myb.max_iter = settings.max_evaluations    
    alg = rbfopt.RbfoptAlgorithm(settings, myb) 
    rbfopt.RbfoptSettings()

    # call optimizer
    val, x, itercount,evalcount, fast_evalcount = alg.optimize()

    
    """ Get Results """   
    myb.evaluate(x) # Recycle all ree
    
    for i,j in zip(myb.var_space['mutable'], x):    
        print(i, '\t',j)

    
    """ store results in json   """
    result = dict()
    result["rees"] = myb.cases
    result["design x"] = [item for item in x]
    result["recovery"] = myb.recovery
    result["purity"] = myb.purity
    
    result["objective"] = {}
    for key,value in myb.recority.items():
        result["objective"][key] = [item for item in value]
    
    result["function value"] = val
    result["constraints"] = myb.constraint
    result["variable space"] = myb.var_space["mutable"]
    result["design"] = {item:jtem for (item,jtem) in zip(myb.var_space["mutable"].keys(),x)}
    
    
    return result
    