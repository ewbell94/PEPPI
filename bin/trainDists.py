import numpy as np
import os
from sklearn.neighbors import KernelDensity
from math import sqrt
from math import log
from random import sample
from random import choice
from random import seed
from pickle import dump
from scipy import std
from scipy.stats import norm
from scipy.stats import beta

supported=[["SPRING","kde"],["STRING","beta"],["SEQ","kde"],["CT","kde"],["SPRINGNEG","kde"]]
truepath=os.path.dirname(os.path.abspath(__file__))

def extractData(supported,truepath):
    data=[]
    labels=[]
    f=open(truepath+"/../lib/trainNB.txt")
    for line in f:
        parts=line.split("\t")
        labels.append(float(parts[1]))
        point=[x.split(",") for x in parts[0].split(";")]
        data.append(point)
    return (data,labels)

def trainDist(data,dist):
    if dist=="beta":
        params=beta.fit(data,floc=0,fscale=1)
        return [params[0],params[1]]
    elif dist=="kde":
        kde=KernelDensity(bandwidth=(max(data)-min(data))/25,kernel="gaussian")
        kde.fit([[i] for i in data])
        return kde
 
def trainModel(trainData,trainLabels,param):
    supported=param
    distParams=[]
    negatives=[trainData[i] for i in range(len(trainLabels)) if trainLabels[i]==0]
    positives=[trainData[i] for i in range(len(trainLabels)) if trainLabels[i]==1]
    for i in range(len(supported)):
        ndata=[x[i] for x in negatives]
        pdata=[x[i] for x in positives]
        featureParams=[]
        for j in range(len(trainData[0][i])):
            print("%s%d"%(supported[i][0],j))
            neg=[float(x[j]) for x in ndata if len(x) > 0 and x[j]!="?" and x[j]!=""]
            pos=[float(x[j]) for x in pdata if len(x) > 0 and x[j]!="?" and x[j]!=""]
            nparams=[]
            pparams=[]
            if supported[i][1]=="beta":
                nlen=float(len(neg))
                plen=float(len(pos))
                neg=[x for x in neg if x!=0.]
                pos=[x for x in pos if x!=0.]
                nparams=[len(neg)/nlen]
                pparams=[len(pos)/plen]
            if supported[i][1]=="kde":
                nparams=trainDist(neg,supported[i][1])
                pparams=trainDist(pos,supported[i][1])
            else:
                nparams+=trainDist(neg,supported[i][1])
                pparams+=trainDist(pos,supported[i][1])
            featureParams.append([nparams,pparams])
            
        if len(featureParams) > 0:
            distParams.append(featureParams)
    return distParams

trainData,trainLabels=extractData(supported,truepath)
trainedModel=trainModel(trainData,trainLabels,supported)
    
dump(trainedModel,open(truepath+"/model_multiD","w"))
