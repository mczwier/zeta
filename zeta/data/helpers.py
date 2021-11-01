'''
Created on Oct 31, 2021

@author: mzwier
'''

_multmap = {'singlet': 1,
            'doublet': 2,
            'triplet': 3,
            'quadruplet': 4,
            'quintet': 5,
            'quintuplet': 5,
            'sextet': 6,
            'sextuplet': 6,
            'septet': 7,
            'septuplet': 7,
            'octet': 8,
            'octuplet': 8,
            'nonet': 9,
            'nonuplet': 9}

def normalize_multiplicity(mult):
    try:
        mult = mult.lower()
        return _multmap[mult]
    except (AttributeError,AttributeError):
        try:
            return int(mult)
        except ValueError:
            raise ValueError('unknown multiplicity {!r}'.format(mult))    


 