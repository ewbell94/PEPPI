from sys import argv
import numpy as np

try:
    pdb1=open(argv[1])
except:
    print("First PDB file not found!")
    exit(1)

try:
    pdb2=open(argv[2])
except:
    print("Second PDB file not found!")
    exit(1)

correspondence={}
coords1=[]
for line in pdb1:
    if line[:4]=="ATOM" and line[12:16].strip()=="CA":
        
        coords1.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        label=line[21]+line[22:26].strip()
        correspondence[label]=[len(coords1)-1,-1]
pdb1.close()        

chain1len=-1
coords2=[]
for line in pdb2:
    if chain1len < 0 and line[:3]=="TER":
        chain1len=len(coords2)
    elif line[:4]=="ATOM" and line[12:16].strip()=="CA":
        coords2.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
        label=line[21]+line[22:26].strip()
        if label in correspondence:
            correspondence[label][1]=len(coords2)-1
pdb2.close()

totallen=len(coords2)
chain2len=totallen-chain1len
        
matrix1=[]
matrix2=[]
for key in correspondence.keys():
    if correspondence[key][1] >= 0:
        matrix1.append(coords1[correspondence[key][0]])
        matrix2.append(coords2[correspondence[key][1]])

matrix1=np.asmatrix(matrix1)
matrix2=np.asmatrix(matrix2)

transvec1=np.sum(matrix1,axis=0)/matrix1.shape[0]
transvec2=np.sum(matrix2,axis=0)/matrix2.shape[0]
#print(transvec1,transvec2)
matrix1-=transvec1
matrix2-=transvec2
matrix1=matrix1.T
matrix2=matrix2.T
H=matrix2*matrix1.T
U,S,Vt=np.linalg.svd(H)
d=np.linalg.det(U*Vt)
R=U*np.diag(np.asarray([1,1,d]))*Vt
rmsd=1/np.sqrt(matrix1.shape[1])*np.linalg.norm(R*matrix1-matrix2)
print(rmsd)

