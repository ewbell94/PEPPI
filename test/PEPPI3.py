#!/usr/bin/env python
import numpy as np
from math import floor

outdir="/nfs/amino-home/ewbell/PEPPI/test"
peppidir="/nfs/amino-home/ewbell/PEPPI"
bindir=peppidir+"/bin"

f=open(outdir+"/protcodeA.csv")
protsA=[line.split(",")[0] for line in f]
f.close()

f=open(outdir+"/protcodeB.csv")
protsB=[line.split(",")[0] for line in f]
f.close()

supported=["SPRING","STRING","SEQ","CT","SPRINGNEG"]
scorelens=[]
scoremat=[]
for prog in supported:
    print(prog)
    singleprog=[]
    f=open(outdir+"/PPI/"+prog+"res.txt")
    for line in f:
        if line.strip()=="" or line[:4]!="prot":
            print("Weird line: "+line.strip())
            continue
        parts=line.strip().split(",")
        chains=parts[0].split("-")
        if chains[0] not in protsA or chains[1] not in protsB:
            print("Irregular pair: "+parts[0])
            continue
        index=protsA.index(chains[0])*len(protsB)+protsB.index(chains[1])
        #print(chains,index)
        scores=[float(n) if n!="?" else np.NaN for n in parts[1:]]
        if len(singleprog)==0:
            singleprog=np.empty([len(protsA)*len(protsB),len(scores)])
            singleprog[:]=np.NaN
            scorelens.append(len(scores))
        singleprog[index]=np.array(scores)
    if len(singleprog)==0:
        singleprog=np.empty([len(protsA)*len(protsB),1])
        singleprog[:]=np.NaN
        scorelens.append(1)
    if len(scoremat)==0:
        scoremat=singleprog
    else:
        scoremat=np.hstack((scoremat,singleprog))
    f.close()
#print(scoremat)
from pickle import load
from sklearn.neighbors import KernelDensity
from math import sqrt
from math import log
from scipy.stats import norm
from scipy.stats import beta

supported=[["SPRING","kde"],["STRING","beta"],["SEQ","kde"],["CT","kde"],["SPRINGNEG","kde"]]
def extractData(resdir,supported):
    points={}
    for i in range(len(supported)):
        dataset=[]
        f=open("%s/%sres.txt"%(resdir,supported[i][0]))
        for line in f:
            parts=line.strip().split(",")
            dataset.append([parts[0].split("/")[1]]+parts[1:])
        f.close()
        for datum in dataset:
            try:
                points[datum[0]].append(datum[1:])
                if i==0:
                    print(datum[0])
                    print("Duplicate key exists")
                    exit(1)
            except:
                points[datum[0]]=[[] for n in range(i)]
                points[datum[0]].append(datum[1:])
        for key in points.keys():
            if len(points[key]) < i+1:
                points[key].append([])
    return points

def calcLR(model,point,supported=supported):
    lr=0.
    for i in range(len(supported)):
        for j in range(len(point[i])):
            if point[i][j]=="nan":
                continue
            if supported[i][1]=="normal":
                lr+=norm.logpdf(float(point[i][j]),model[i][j][1][0],sqrt(model[i][j][1][1]))-norm.logpdf(float(point[i][j]),model[i][j][0][0],sqrt(model[i][j][0][1]))
            elif supported[i][1]=="beta":
                if point[i][j]=="0" or point[i][j]=="0.0":
                    lr+=log((1.-model[i][j][1][0])/(1.-model[i][j][0][0]))
                else:
                    lr+=log(model[i][j][1][0]/model[i][j][0][0])
                    lr+=beta.logpdf(float(point[i][j]),model[i][j][1][1],model[i][j][1][2])-beta.logpdf(float(point[i][j]),model[i][j][0][1],model[i][j][0][2])
            elif supported[i][1]=="kde":
                lr+=model[i][j][1].score([[float(point[i][j])]])-model[i][j][0].score([[float(point[i][j])]])
    return lr

model=load(open("/nfs/amino-home/ewbell/PEPPI/bin/model_multiD"))
#print(model)
scoredPoints=[] 

f=open(outdir+"/allres.txt","w")
for i in range(len(scoremat)):
    if i % 100000 == 0:
        print(str(i)+"/"+str(len(scoremat)))
    if sum(np.isnan(scoremat[i])) == len(scoremat[i]):
        continue
    aind=int(floor(i/len(protsB)))
    bind=i-(aind*len(protsB))
    f.write(protsA[aind]+"-"+protsB[bind])
    j=0
    point=[]
    for l in scorelens:
        progscore=[str(n) for n in scoremat[i][j:j+l]]
        point.append(progscore)
        f.write(";"+",".join(progscore))
        j+=l
    f.write("\n")
    if len(point) != len(supported):
        print("Skipped "+protsA[aind]+"-"+protsB[bind])
        continue
    scoredPoints.append((protsA[aind]+"-"+protsB[bind],calcLR(model,point)))
f.close()

scoredPoints.sort(key=lambda x:x[1])

outfile=open(outdir+"/LR.csv","w")
for point in scoredPoints[::-1]:
    outfile.write("%s,%.3f\n"%point)

outfile.close()
