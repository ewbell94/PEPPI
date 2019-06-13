#!/usr/bin/env python

import sys
import numpy as np
sys.path.append("/nfs/amino-home/zhng/zhanglab_programs/TACOS/source/PDB/TMalign")
import tmalign

threetoone={"GLY":"G","ALA":"A","LEU":"L","MET":"M","PHE":"F",
            "TRP":"W","LYS":"K","GLN":"Q","GLU":"E","SER":"S",
            "PRO":"P","VAL":"V","ILE":"I","CYS":"C","TYR":"Y",
            "HIS":"H","ARG":"R","ASN":"N","ASP":"D","THR":"T"}

f=open(sys.argv[1])
seqA=""
coordA=[]
for line in f:
    if line[0:4]=="ATOM" and line[12:16].strip()=="CA":
        seqA+=threetoone[line[17:20]]
        coordA.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
f.close()
coordA=np.asarray(coordA)

f=open(sys.argv[2])
seqB=""
coordB=[]
for line in f:
    if line[0:4]=="ATOM" and line[12:16].strip()=="CA":
        seqB+=threetoone[line[17:20]]
        coordB.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
f.close()
coordB=np.asarray(coordB)

result=tmalign._tmalign(seqB,coordB,seqA,coordA)
print("%.5f"%result[8])

if len(argv) > 3 and argv[3]=="-o":
    
