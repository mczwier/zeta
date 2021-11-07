'''
Created on Nov 5, 2021

@author: mzwier
'''

class PropertySet(dict):
    def __init__(self, iterable=None, method=None, provenance=None):
        super().__init__(iterable or {})
        
        self.method = method
        self.provenance = provenance
