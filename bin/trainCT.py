from sklearn.neural_network import MLPClassifier
from pickle import dump
import os

data=[]

def trainModel(trainData,trainLabels,param):
    model=MLPClassifier(hidden_layer_sizes=param)
    model.fit(trainData,trainLabels)
    return model

labels=[]
thispath=os.path.dirname(os.path.abspath(__file__))
f=open(thispath+"/../lib/CTtrainvec.txt")
for line in f:
    parts=line.split("\t")
    data.append([float(x) for x in parts[0].split(",")])
    labels.append(float(parts[1]))
f.close()

model=trainModel(data,labels,(1000))
dump(model,open(thispath+"/CTNN","w"))
