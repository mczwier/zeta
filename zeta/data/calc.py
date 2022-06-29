'''
Created on Oct 31, 2021

@author: mzwier
'''

class CalcNode:
    def __init__(self):
        self.parents = []
        self.provenance = {}
        self.geometry = None
        self.method = {}
        self.results = {}
        
class Provenance(dict):
    pass

