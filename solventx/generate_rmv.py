# -*- coding: utf-8 -*-
"""
Created on Sun May 24 12:16:35 2020

@author: ciloeje
"""
import json
import numpy as np

def generate(confDict, N):
    
    reeComps        = confDict['compositions']
    modulesData     = confDict['modules']      
    ree             = modulesData["input"]

    upper           = [reeComps[i]['upper'] for i in ree]
    lower           = [reeComps[i]['lower'] for i in ree]
    
    ree_mass_dict   = dict()
    
    for index in range(N):
        ree_mass_dict[str(N)] = {}
        ree_mass = [np.random.uniform(i,j) for i,j in zip(lower, upper)]  
        for item,jtem in zip(ree, ree_mass):
            ree_mass_dict[str(N)][item] = jtem
    
    with open('data/json/cases.json','w') as json_file:
        json.dump(ree_mass_dict, json_file)   
        
    return ree_mass_dict
        
