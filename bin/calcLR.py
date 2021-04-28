from pickle import load
from sklearn.neighbors import KernelDensity
from math import sqrt
from math import log
from sys import argv
from scipy.stats import norm
from scipy.stats import beta

#supported=[["SPRING","kde"],["STRING","beta"],["SEQ","kde"],["CT","kde"],["TMSEARCH","kde"]]
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
            if point[i][j]=="?":
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
resfile=open(argv[1])
scoredPoints=[]
for line in resfile:
    parts=line.strip().split(";")
    point=[x.split(",") for x in parts[1:]]
    if len(point) != len(supported):
        print("Skipped "+parts[0])
        continue
    print(parts[0])
    scoredPoints.append((parts[0],calcLR(model,point)))

scoredPoints.sort(key=lambda x:x[1])

outfile=open(argv[2],"w")
for point in scoredPoints[::-1]:
    outfile.write("%s,%.3f\n"%point)



