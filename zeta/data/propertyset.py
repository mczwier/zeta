'''
Created on Nov 5, 2021

@author: mzwier
'''

class PropertySet(dict):
    '''A PropertySet contains the results of calculations at a given geometry.
    Results are stored dict-like, and references to method and provenance
    objects are also provided as attributes.'''
    def __init__(self, iterable=None, method=None, provenance=None):
        super().__init__(iterable or {})
        
        self.method = method
        self.provenance = provenance
