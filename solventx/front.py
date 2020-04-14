# -*- coding: utf-8 -*-


""" imports """
import solventx_B as sx
import utilities as util



""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" 
here the main portion of the code begins. There are three evaluate options:
    evaluate_open(): once thru, rapid computation: approximate results
    evaluate_loop(lim=val): partial recycle implemented; slower computation
    evaluate_loop(lim=0): full recycle implemented - most accurate result;
    slowest computation
    
    The definition of recycle differs between option 1 and the others: in option 
    1, recycle is the ratio of the scrub feed to the primary extraction feed
    this value can be significantly greater than 1; in the other options, 
    recycle is the actual fraction of the strip exit stream recycled to the 
    scrub column, and takes values between 0 and 1
    
    This is work in progress. model still needs to be updated for accurate
    results
"""    
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
   
    
""" Provide variable list   """

#config_file = "solventx\\design.json"
config_file = "D:\\github\\solventx\\solventx\\design_.json"
""" instantiate solvent extraction object """

obj = sx.solventx(config_file) 

obj.get_process()
print ('\nnumber of components', obj.num_input ,'\n')


""" create variable space parameters """
obj.create_var_space(input_feeds=1)

#for key, value in obj.var_space['mutable'].items():
#    print (key, ': ', value)
#print()
#
#for key, value in obj.var_space['immutable'].items():
#    print (key, ': ', value)
#print()
#
#for key, value in obj.mod_space.items():
#    print (key, ': ', value)
#print()
#
#for key, value in obj.x_space.items():
#    print (key, ': ', value)
#print()


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" Action methods """
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


""" Option 1 : ~ 1 sec; use for initial training 
    You improve purity by directly adding pure strip target elements to the scrub column
    you set the recycle ratio limit to a maximum of 1, even though it can be higher than 1
    you require high recovery to ensure your column is actually doing work"""

#t = time.time()
obj.evaluate_open(obj.design_variables) # 
#print('Option 1 time', time.time() - t)


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" Processing outputs """
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

print()
""" Recovery - dictionary of rare earth output stream fractional recovery values by column
      recovery values = amount of ree in raffinate stream/amount of ree in original feed stream """
obj.recovery # 


""" Purity - dictionary of rare earth output stream fractional compositions by column
     # purity values = amount of ree in raffinate/total ree in raffinate stream """
obj.purity # 


""" Target - dictionary of target rare earth IDs by final column
"""
obj.target_rees # 


""" Streams - dictionary of all streams leaving the three columns in each module
"""
obj.streams # 


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" visualize outcomes """
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
### plot results - only for testing config.valid_processes == f
### test to see how recovery and purity vary with Extractant concentration
#
#import numpy as np
#import pandas as pd
#recov = []
#pur = []
#for ha in np.arange(0.3,0.8, 0.05):
#    obj.design_variables[obj.var_space['mutable']['(HA)2(org)-0']] = ha
#    obj.evaluate_open(obj.design_variables)
#    
#    recov.append(obj.recovery['Strip-0'])
#    pur.append(obj.purity['Strip-0'])
#
#df = pd.DataFrame(recov, columns=['Nd Recovery','Pr Recovery'])   
#df['Purity'] = pur
#df['Index'] = np.arange(0.3,0.8, 0.05)
#
#
#import matplotlib.pyplot as plt
#import pandas as pd
#
## gca stands for 'get current axis'
#ax = plt.gca()
#
#df.plot(kind='line',x='Index',y='Nd Recovery',ax=ax)
#df.plot(kind='line',x='Index',y='Pr Recovery',ax=ax)
#df.plot(kind='line',x='Index',y='Purity',ax=ax)
#
#plt.show()

