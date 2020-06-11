# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 12:43:20 2020

@author: splathottam
"""

DEFAULT_module0_x = [0.65,1e-4,1e-4,2,2,4,6,.8,8,10,4] #
DEFAULT_module1_x =      [1e-5,1e-5,1.2,6,1,4,0.51,1,8,2]
DEFAULT_module2_x =      [1e-5,1e-5,1.2,6,1,1,0.01,4,9,2]

mutable_var_names   = ['H+ Extraction', 'H+ Scrub', 'H+ Strip', 'OA Extraction', 'OA Scrub', 'OA Strip', 'Recycle', 'Extraction', 'Scrub', 'Strip']


valid_processes = {'a':{'input':['Nd','Pr','Ce','La'],
                        'strip':['Nd','Ce'],
                        'modules':{"0":{"strip_group":["Nd","Pr"], "ext_group":["Ce","La"],"x":DEFAULT_module0_x, "mvn":mutable_var_names},
                                  "1":{"strip_group":["Nd"],"ext_group":["Pr"],"x":DEFAULT_module1_x, "mvn":mutable_var_names},
                                  "2":{"strip_group":["Ce"],"ext_group":["La"],"x":DEFAULT_module2_x, "mvn":mutable_var_names}}
                        },
                   'b':{'input':['Nd','Pr'],
                        'strip':['Nd'],
                        'modules':{"1":{"strip_group":["Nd"],"ext_group":["Pr"],"x":DEFAULT_module0_x[:1]+DEFAULT_module1_x, "mvn":mutable_var_names}}
                        },
                   'c':{'input':['Ce','La'],
                        'strip':['Ce'],
                        'modules':{"2":{"strip_group":["Ce"],"ext_group":["La"],"x":DEFAULT_module0_x[:1]+DEFAULT_module2_x, "mvn":mutable_var_names}}
                        },
                   'd':{'input':['Nd','Pr','Ce','La'],
                        'strip':['Nd'],
                        'modules':{"0":{"strip_group":["Nd","Pr"],"ext_group":["Ce","La"],"x":DEFAULT_module0_x, "mvn":mutable_var_names},
                                  "1":{"strip_group":["Nd"],"ext_group":["Pr"],"x":DEFAULT_module1_x, "mvn":mutable_var_names}}
                        },                   
                   'e':{'input':['Nd','Pr','Ce','La'],
                        'strip':['Ce'],
                        'modules':{"0":{"strip_group":["Nd","Pr"],"ext_group":["Ce","La"],"x":DEFAULT_module0_x, "mvn":mutable_var_names},
                                  "2":{"strip_group":["Ce"],"ext_group":["La"],"x":DEFAULT_module2_x, "mvn":mutable_var_names}}
                        },                   
                   'f':{'input':['Nd','Pr','Ce','La'],
                        'strip':[''],
                        'modules':{"0":{"strip_group":["Nd","Pr"],"ext_group":["Ce","La"],"x":DEFAULT_module0_x, "mvn":mutable_var_names}}
                        },                   
                     }
