'''
Created on Oct 31, 2021

@author: mzwier
'''

class CalcStep:
    def __init__(self, *, provenance=None, method=None, geometry=None, properties=None):
        
        # Any of these may be references to another CalcStep's corresponding object
        self.provenance = provenance or {}
        self.method = method or {}
        self.geometry = geometry
        self.properties = properties or {}

class CalcTree:
    def __init__(self, calcs, hints=None):
        # Root calculations
        self.calcs = calcs
        
        # Hints about how to process this tree
        # Examples might include 'optimization', 'hessian', 'scan'
        self.hints = set(hints) if hints else set()


        
    