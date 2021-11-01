'''
Created on Oct 31, 2021

@author: mzwier
'''

from enum import IntEnum
from .method import QMMethod

class CalcType(IntEnum):
    UNKNOWN = 0
    
    ENERGY = 10
    
    OPTIMIZATION = 20
    
    HESSIAN = 30
    
    EXCITATION = 40
    
    TRANSITION_STATE = 50
    
    REACTION_PATH = 60

class Calculation:
    def __init__(self, method=None):
        self.calc_type = CalcType.UNKNOWN
        
        self.prev = None
        self.next = None
        
        self.provenance = {}
        self.method = method
        self.atom_data = None
        self.geometries = None
        