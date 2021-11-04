'''
Created on Oct 31, 2021

@author: mzwier
'''


class Method:
    def __init__(self):
        self.name = None
        self.details = {}

class QMMethod(Method):    
    def __init__(self):
        super().__init__()
        self.charge = None
        self.multiplicity = None
        self.basis_name = None
        self.aux_basis_names = {}
        self.additional_keywords = []
                
        
        