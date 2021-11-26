'''Data model.'''


from . import calc, geom, helpers, method, elements, io

from .calc import Calculation, CalcStep
from .geom import Geometry, normalize_atoms
from .helpers import normalize_multiplicity

import numpy as np

# The version of this data model. Major schema updates can break stored
# data in such a way that an explicit conversion is required. Minor schema
# updates, on the other hand, might be handled by special cases in the code,
# or have no functional impact at all.
schema_major_version = 1
schema_minor_version = 1

atomicnumber_dtype = np.uint8

del np