'''
Created on Jan 13, 2023

@author: mzwier
'''

# Quantum numbers

_l_to_subshell = {0: 's',
                 1: 'p',
                 2: 'd',
                 3: 'f',
                 4: 'g',
                 5: 'h',
                 6: 'i',
                 7: 'k',
                 8: 'l',
                 9: 'm',
                 10: 'n',
                 11: 'o',
                 12: 'q',
                 13: 'r',
                 14: 't',
                 15: 'u',
                 16: 'v',
                 17: 'w',
                 18: 'x',
                 19: 'y',
                 20: 'z'}

_subshell_to_l = dict((v,k) for (k,v) in _l_to_subshell.items())

def l_to_subshell(l):
    return _l_to_subshell[l]

def subshell_to_l(subshell):
    return _subshell_to_l[subshell.lower()]