'''
Created on Oct 31, 2021

@author: mzwier
'''

from enum import IntEnum

class CalcType(IntEnum):
    UNKNOWN = 0
    OTHER   = 1
    
    ENERGY = 10
    
    OPTIMIZATION = 20
    
    HESSIAN = 30
    
    EXCITATION = 40
    
    TRANSITION_STATE = 50
    
    REACTION_PATH = 60

class Calculation:
    def __init__(self):
        self.calc_type = CalcType.UNKNOWN
        
        self.geometries = None
        