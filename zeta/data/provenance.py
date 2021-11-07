'''
Created on Nov 5, 2021

@author: mzwier
'''

class Provenance(dict):
    def __init__(self, iterable=None):
        super().__init__(iterable or {})
