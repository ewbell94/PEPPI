#!/usr/bin/env python
import os
import sys
import numpy as np
from pdbchains import PDBchains
#seq=[]
native=PDBchains()
native.Read('12asA.pdb', ca_only=True,allow_alt_loc =False, allow_insert = False)
header=native.template_headers
seq=native.sequence
print(header)
