import cython
from libc.stdlib cimport atof

@cython.boundscheck(False)
@cython.wraparound(False)
def _assign_positions_cy(lines, Py_ssize_t natoms, double[:,:] coords):
    for iatom in range(natoms):
        coords[iatom,0] = float(lines[iatom][0:18])
        coords[iatom,1] = float(lines[iatom][18:36])
        coords[iatom,2] = float(lines[iatom][36:54])