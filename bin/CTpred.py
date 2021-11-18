import os
from pickle import load
from sys import argv
#This script simply predicts interaction probability from the precalculated CTNN

<<<<<<< HEAD
model=load(open("/home/ewbell/PEPPI/bin/CTNN"))
=======
thispath=os.path.dirname(os.path.abspath(__file__))
model=load(open(thispath+"/CTNN"))
>>>>>>> master
f=open(argv[1])
point=[float(i) for i in f]
f.close()
posp=model.predict_proba([point])[0][1]
print(posp)
