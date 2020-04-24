#!/usr/bin/env python

import sys
import numpy as np
sys.path.append("/nfs/amino-home/zhng/zhanglab_programs/TACOS/source/PDB/TMalign")
import tmalign

threetoone={"GLY":"G","ALA":"A","LEU":"L","MET":"M","PHE":"F",
            "TRP":"W","LYS":"K","GLN":"Q","GLU":"E","SER":"S",
            "PRO":"P","VAL":"V","ILE":"I","CYS":"C","TYR":"Y",
            "HIS":"H","ARG":"R","ASN":"N","ASP":"D","THR":"T",
            "SEC":"U","PYL":"O"}

onetothree={value:key for key,value in threetoone.items()}

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
print("%.5f"%result[6])
#print(result)
if len(sys.argv) > 3 and sys.argv[3]=="-o":
    rotmat=np.asmatrix(result[0][:3,:3])
    transvec=np.asmatrix(result[0][:3,3])
    coordA=(np.asmatrix(coordA)-transvec)*rotmat
    
    alignment=result[2]
    queryalign=result[7]
    tempalign=result[5]
    queryi=0
    tempi=0
    alignedquery=[]
    alignedtemp=[]
    for i in range(len(alignment)):
        if alignment[i]!=" ":
            alignedquery.append(queryi)
            alignedtemp.append(tempi)
        if queryalign[i]!="-":
            queryi+=1
        if tempalign[i]!="-":
            tempi+=1

    #print(queryi,tempi)
    #print(len(alignedquery),len(alignedtemp))
    o=open(sys.argv[4],"w")
    a=open("%s_all"%sys.argv[4],"w")
    for i in range(coordA.shape[0]):
        pdbline="ATOM  %5d  CA  %s A %3d    %8.3f%8.3f%8.3f\n"%(i+1,onetothree[seqA[i]],i+1,coordA[i,0],coordA[i,1],coordA[i,2])
        a.write(pdbline)
        if i in alignedquery:
            o.write(pdbline)
        
    a.write("TER\n")
    o.write("TER\n")
    for i in range(coordB.shape[0]):
        pdbline="ATOM  %5d  CA  %s B %3d    %8.3f%8.3f%8.3f\n"%(i+5001,onetothree[seqB[i]],i+1,coordB[i,0],coordB[i,1],coordB[i,2])
        a.write(pdbline)
        if i in alignedtemp:
            o.write(pdbline)
    a.write("TER\n")
    o.write("TER\n")
    a.close()
    o.close()

        
