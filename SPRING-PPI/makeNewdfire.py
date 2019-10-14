import matplotlib.pyplot as p
import numpy as np
from copy import deepcopy

def polyTail(row):
    localmax=0
    closeint=-1
    while row[localmax]-row[localmax+1]<0.001:
        if closeint < 0 and row[localmax]>0.01:
            closeint=localmax-1
        localmax+=1
    if closeint < 0:
        closeint=localmax-1
    farint=2*localmax-closeint
    X=np.asmatrix([[localmax**2,localmax,1.],[closeint**2,closeint,1.],[farint**2,farint,1.]])
    b=np.asmatrix([[row[localmax]],[0],[0]])
    vals=np.linalg.pinv(X)*b
    a=vals[0,0]
    b=vals[1,0]
    c=vals[2,0]
    out=deepcopy(row)
    for i in range(localmax):
        out[i]=vals[0,0]*i*i+vals[1,0]*i+vals[2,0]
    return out

f=open("bin/dfire.txt")
flatd=[float(line) for line in f]
f.close()

dfire=[]
for i in range(21):
    box=[]
    for j in range(21):
        row=flatd[20*21*i+20*j:20*21*i+20*j+20]
        box.append(row)
    dfire.append(box)

for i in range(20):
    for j in range(20):
        #p.plot([x for x in range(len(dfire[i][j]))],[x for x in dfire[i][j]])
        dfire[i][j]=polyTail(dfire[i][j])
        #p.plot([x for x in range(len(dfire[i][j]))],[x for x in dfire[i][j]])
        #p.ylim(-0.5,0.5)
        #p.xlabel("Inter-residue Distance")
        #p.ylabel("DFIRE Score")
        #p.show()


flatnewd=[]
for i in range(21):
    for j in range(21):
        flatnewd+=dfire[i][j]

g=open("bin/newdfire.txt","w")
for x in flatnewd:
    g.write(str(x)+"\n")
g.close()
