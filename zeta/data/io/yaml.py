'''
Created on Nov 25, 2021

@author: mzwier
'''

import yaml

def write_yaml(calcs, outfile):
    '''Write the given calculation graph to `outfile` in yaml format.'''
    yaml.dump(calcs, outfile, Dumper=yaml.CDumper)
    
def read_yaml(infile):
    '''Read a calculation graph from `infile` and return.'''
    return yaml.load(infile)
