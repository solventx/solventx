# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 20:48:45 2020

@author: ciloeje
"""

# -*- coding: utf-8 -*-

import time
import json
import rbfopt


""" local imports """
from solventx import solventx as sx
from solventx import utilities as util


#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" 
Optimize solvent extraction configuration
"""    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
def optimize(confDict, confEnvDict, case):
    
    """ instantiate solvent extraction object """
    ree_mass = [item for item in case.values()]
    myb = sx.solventx(confDict, confEnvDict, ree_mass) 

    # Specify optimization settings
    settings = rbfopt.RbfoptSettings(max_evaluations=10, init_strategy='lhd_corr', target_objval=0, max_noisy_restarts=3) #'lhd_corr')#, minlp_solver_path=minlp_solver, nlp_solver_path=nlp_solver) 

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
    result["rees"] = case
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





#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" 
Section A: Test process simulation
"""    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

""" build configs """
config_file = ".\\design.json"
confDict = util.read_config(config_file)

config_env_file = ".\\env_design_config.json"
confEnvDict = util.get_env_config_dict(config_env_file)

""" generate cases """
cases = util.generate(confDict,1)

""" select case and create object """
ree_mass = [item for item in cases['0'].values()]
obj = sx.solventx(confDict, confEnvDict, ree_mass)

""" evaluate """
obj.evaluate(obj.design_variables) # 

""" Reward """
size = len(ree_mass)
reward = obj.objective(size)

""" Examine Specs """
purity, recovery = obj.constraintfun()
recority = obj.recority




#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" 
Section B: Call the Standard Optimizer
"""    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#t0 = time.time()
#
#config_file = ".\\design.json"
#confDict = util.read_config(config_file)
#
#config_env_file = ".\\env_design_config.json"
#confEnvDict = util.get_env_config_dict(config_env_file)

#""" generate cases """
#cases = util.generate(confDict,3)
#
#t1 = time.time()
##for now, pick only the first cases - later cases will use parallel submission
#resjsn = optimize(confDict, confEnvDict, cases['0'])
##
#
#with open('solventx/data/json/results.json', 'w') as json_file:
#    json.dump(resjsn, json_file)
#
#t2 = time.time()
#print('\nPrep Time', t1-t0)
#print('Optim Time', t2-t1)
#print('Total Time', t2-t0)
#
#
#print('\nRecovery: ', resjsn["recovery"]['Strip-0'])
#print('Purity: ', resjsn["purity"]['Strip-0'])


