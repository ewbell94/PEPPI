#!/usr/bin/env python

import os
from alignment import TMalign,TMscore
transfromMatrix, tmout = TMalign('rank1_chain1.pdb', '12asA.pdb',modelIsFile=True,nativeIsFile=True)
tmscorematrix, tmscore = TMscore('rank1_chain1.pdb', '12asA.pdb',modelIsFile=True,nativeIsFile=True)
print(transfromMatrix)
print(tmout)
print(tmscorematrix)
print(tmscore)

