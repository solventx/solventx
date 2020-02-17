# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:55:46 2019

@author: ciloeje
"""

class result_struct():

    def __init__(self, x, status, msg, fun):
        self.x = x
        self.status = status
        self.message = msg
        self.fun     = fun