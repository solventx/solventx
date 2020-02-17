# -*- coding: utf-8 -*-
"""
front page

Created on Wed Aug 14 08:07:36 2019

@author: ciloeje
"""


""" imports """
import solventx as sx


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

# for option 1 (eval_open)
x1 = [0.65,1e-5,1e-5,1.2,5.9,1.2,1.6,5,6,6,4]
x2 = [1e-5,1e-5,1.2,6,1,4,0.51,1,8,2]
x3 = [1e-5,1e-5,1.2,6,1,1,0.01,4,9,2]

## for options 2 and 3 (eval_loop)
#x1 = [0.65,1e-5,1e-5,1.2,5.9,1.2,1.6,0.89,6,6,4]
#x2 = [1e-5,1e-5,1.2,6,1,4,0.89,1,8,2]
#x3 = [1e-5,1e-5,1.2,6,1,1,0.89,4,9,2]

x = x1+x2+x3 # input provided by AI environ


""" instantiate solvent extraction object """
obj = sx.solventx() 

""" create variable space parameters """
obj.create_var_space(x,components=len(obj.ree_mass), input_feeds=1)
# myblack.var_names


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" Action methods """
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


""" Option 1 : ~ 1 sec; use for initial training"""
#t = time.time()
obj.evaluate_open(x) # 
#print('Option 1 time', time.time() - t)


""" Option 2 - partially closed recycle loop """
#t = time.time()
#myblack.evaluate_loop(x, lim=len(myblack.ree_strip)) # Recycle strip target ree
#print('Option 2 time', time.time() - t)


""" Option 3 - closed recycle loop: ~15 secs; use for main training """
#obj.evaluate_loop(x,lim=0) # Recycle all ree


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" Processing outputs """
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


""" Recovery - dictionary of rare earth output stream fractional recovery values by column
      recovery values = amount of ree in raffinate stream/amount of ree in original feed stream """
obj.recovery_ # closed loop recycle calculation - more accurate
obj.recovery # open loop recycle calculation - less accurate


""" Purity - dictionary of rare earth output stream fractional compositions by column
     # purity values = amount of ree in raffinate/total ree in raffinate stream """
obj.purity_ # closed loop recycle calculation - more accurate
obj.purity # open loop recycle calculation - less accurate


""" Raffinates - dictionary of output stream stream molar amounts/flows by column """
obj.raffinates_
obj.raffinates


""" max_recov - dictionary of the rare earth recovery value that corresponds to the highest purity  by column """
obj.max_recov_
obj.max_recov


""" max_pur - dictionary of maximum purity values by column """
obj.max_pur_
obj.max_pur

""" ree_max - dictionary of rare earths that have the maximum composition by column """
obj.ree_max_
obj.ree_max # both return the same REE acronyms
     



""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""" visualize what you can't see """
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
## plot results - will implemen later 


