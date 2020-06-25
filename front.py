# -*- coding: utf-8 -*-

import time
import json
#import rbfopt


""" local imports """
from solventx import solventx as sx
from solventx import utilities as util



""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" 
Optimize solvent extraction configuration
"""    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""



""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" 
Main Section of Code: Call the Optimizer
"""    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
for item in ['a','b','c','d','e','f']:

    t0 = time.time()
    
    config_file = ".\\design_"+item+".json"
    confDict = util.read_config(config_file)
    
    config_env_file = ".\\env_design_config.json"
    confEnvDict = util.get_env_config_dict(config_env_file)
    
    """ generate cases """
    cases = util.generate(confDict,3)
    
    t1 = time.time()
    
    #for now, pick only the first cases - later cases will use parallel submission
    ree_mass = [item for item in cases['0'].values()]
    
    """ instantiate solvent extraction object """
    myb = sx.solventx(confDict, confEnvDict, ree_mass) 
    myb.cases = cases['0']
    iters = 15 # Number of iteration, at least 100 
    
    resjsn = util.optimize(myb,iters)
    #
    
    with open('solventx/data/json/results_'+item+'.json', 'w') as json_file:
        json.dump(resjsn, json_file)
    
t2 = time.time()
print('\nPrep Time', t1-t0)
print('Optim Time', t2-t1)
print('Total Time', t2-t0)


#print('\nRecovery: ', resjsn["recovery"]['Strip-0'])
#print('Purity: ', resjsn["purity"]['Strip-0'])
#
#print('\nRecovery 2: ', resjsn["recovery"]['Strip-0'])
#print('Purity 2: ', resjsn["purity"]['Strip-2'])

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" 
Secondary Section of Code: When not optimizing
"""    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#
##for item in ['a','b','c','d','e','f']:
#for item in ['e']:
#    config_file = ".\\design_"+item+".json"
#    confDict = util.read_config(config_file)
#    
#    config_env_file = ".\\env_design_config.json"
#    confEnvDict = util.get_env_config_dict(config_env_file)
#    
#    """ generate cases """
#    cases = util.generate(confDict,3)
#    
#    """ select case and create object """
#    ree_mass = [item for item in cases['0'].values()]
#    obj = sx.solventx(confDict, confEnvDict, ree_mass)
#    
#    """ evaluate """
#    obj.evaluate(obj.design_variables) # 
#    fun = obj.objective()
#    pur, rec = obj.constraintfun()
#    
#    print ('Design case ', item)
#    print('Recovery',obj.recovery)
##    print ('Purity\n', obj.purity)
#    print ('Reward',obj.objective())
#    print('$$$$$$$$$$$$$$$$$$$$\n')









""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" Action methods """
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


""" Option 1 : ~ 1 sec; use for initial training 
    You improve purity by directly adding pure strip target elements to the scrub column
    you set the recycle ratio limit to a maximum of 1, even though it can be higher than 1
    you require high recovery to ensure your column is actually doing work"""

#t = time.time()
#obj.evaluate_open(obj.design_variables) # 
#print('Option 1 time', time.time() - t)














#
##for key, value in obj.var_space['mutable'].items():
##    print (key, ': ', value)
##print()
##
##for key, value in obj.var_space['immutable'].items():
##    print (key, ': ', value)
##print()
##
##for key, value in obj.mod_space.items():
##    print (key, ': ', value)
##print()
##
##for key, value in obj.x_space.items():
##    print (key, ': ', value)
##print()
#
#
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#""" Action methods """
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#
#
#""" Option 1 : ~ 1 sec; use for initial training 
#    You improve purity by directly adding pure strip target elements to the scrub column
#    you set the recycle ratio limit to a maximum of 1, even though it can be higher than 1
#    you require high recovery to ensure your column is actually doing work"""
#
#t = time.time()
#obj.evaluate_open(obj.design_variables) # 
#print('Option 1 time', time.time() - t)
#
#
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#""" Processing outputs """
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#
#print()
#""" Recovery - dictionary of rare earth output stream fractional recovery values by column
#      recovery values = amount of ree in raffinate stream/amount of ree in original feed stream """
#obj.recovery # 
#
#
#""" Purity - dictionary of rare earth output stream fractional compositions by column
#     # purity values = amount of ree in raffinate/total ree in raffinate stream """
#obj.purity # 
#
#
#""" Target - dictionary of target rare earth IDs by final column
#"""
#obj.target_rees # 
#
#
#""" Streams - dictionary of all streams leaving the three columns in each module
#"""
#obj.streams # 
#
#
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#""" visualize outcomes """
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#### plot results - only for testing config.valid_processes == f
#### test to see how recovery and purity vary with Extractant concentration
##
##import numpy as np
##import pandas as pd
##recov = []
##pur = []
##for ha in np.arange(0.3,0.8, 0.05):
##    obj.design_variables[obj.var_space['mutable']['(HA)2(org)-0']] = ha
##    obj.evaluate_open(obj.design_variables)
##    
##    recov.append(obj.recovery['Strip-0'])
##    pur.append(obj.purity['Strip-0'])
##
##df = pd.DataFrame(recov, columns=['Nd Recovery','Pr Recovery'])   
##df['Purity'] = pur
##df['Index'] = np.arange(0.3,0.8, 0.05)
##
##
##import matplotlib.pyplot as plt
##import pandas as pd
##
### gca stands for 'get current axis'
##ax = plt.gca()
##
##df.plot(kind='line',x='Index',y='Nd Recovery',ax=ax)
##df.plot(kind='line',x='Index',y='Pr Recovery',ax=ax)
##df.plot(kind='line',x='Index',y='Purity',ax=ax)
##
##plt.show()

