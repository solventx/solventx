# -*- coding: utf-8 -*-
"""
front page

Created on Wed Aug 14 08:07:36 2019

@author: ciloeje
"""


""" imports """
import solventx as sx
import config
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

config_file = "D:\\github\\solventx\\v4\\solventx-refactoring_siby_nw\\design.json"
""" instantiate solvent extraction object """

obj = sx.solventx(config_file) 
#print('object:',obj)

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

#design_variables = obj.x


""" Option 1 : ~ 1 sec; use for initial training 
    You improve purity by directly adding pure strip target elements to the scrub column
    you set the recycle ratio limit to a maximum of 1, even though it can be higher than 1
    you require high recovery to ensure your column is actually doing work"""

#t = time.time()
obj.evaluate_open(obj.design_variables) # 
#print('Option 1 time', time.time() - t)
##


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" Processing outputs """
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

print()
""" Recovery - dictionary of rare earth output stream fractional recovery values by column
      recovery values = amount of ree in raffinate stream/amount of ree in original feed stream """
obj.recovery # open loop recycle calculation - less accurate


""" Purity - dictionary of rare earth output stream fractional compositions by column
     # purity values = amount of ree in raffinate/total ree in raffinate stream """
obj.purity # open loop recycle calculation - less accurate



#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#""" visualize what you can't see """
#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
### plot results - will implemen later 


