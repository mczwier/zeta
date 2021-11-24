import cython

@cython.boundscheck(False)
@cython.wraparound(False)
def _assign_positions_cy(lines, Py_ssize_t natoms, Py_ssize_t fieldwidth, double[:,:] coords):
    for iatom in range(natoms):
        coords[iatom,0] = float(lines[iatom][0:fieldwidth])
        coords[iatom,1] = float(lines[iatom][fieldwidth:2*fieldwidth])
        coords[iatom,2] = float(lines[iatom][2*fieldwidth:3*fieldwidth])
