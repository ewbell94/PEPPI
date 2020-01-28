import os.path
from sys import argv
from math import log
from math import sqrt
from scipy.stats import norm
from scipy.stats import beta

peppidir="/nfs/amino-home/ewbell/PEPPI"
supported=[["SPRING","normal"],["STRING","beta"]]

def loadModel():
    f=open(peppidir+"/bin/model.csv")
    model=[]
    for line in f:
        module=[]
        for sub in line.strip().split(";"):
            submod=[]
            for pn in sub.split("/"):
                submod.append([float(x) for x in pn.split(",")])
            module.append(submod)
        model.append(module)
    f.close()
    return model

def evaluatePoint(model,point,supported=supported):
    lr=1.
    for i in range(len(supported)):
        for j in range(len(point[i])):
            if point[i][j]=="?":
                continue
            if supported[i][1]=="normal":
                lr*=norm.pdf(point[i][j],model[i][j][1][0],sqrt(model[i][j][1][1]))/norm.pdf(point[i][j],model[i][j][0][0],sqrt(model[i][j][0][1]))
            elif supported[i][1]=="beta":
                if point[i][j]==0.0:
                    lr*=(1.-model[i][j][1][0])/(1.-model[i][j][0][0])
                else:
                    lr*=model[i][j][1][0]/model[i][j][0][0]
                    lr*=beta.pdf(point[i][j]/1000.,model[i][j][1][1],model[i][j][1][2])/beta.pdf(point[i][j]/1000.,model[i][j][0][1],model[i][j][0][2])
    return log(lr)

model=loadModel()
#print(model)
for inputdir in argv[1:]:
    datapoint=[]
    for prog in supported:
        res=[]
        if os.path.exists(inputdir+"/"+prog[0]+"/res.txt"):
            f=open(inputdir+"/"+prog[0]+"/res.txt")
            for line in f:
                if line.strip()=="?":
                    res.append("?")
                else:
                    res.append(float(line))
            f.close()
        datapoint.append(res)
    #print(datapoint)
    try:
        print(evaluatePoint(model,datapoint))
    except:
        print(-1.0)
