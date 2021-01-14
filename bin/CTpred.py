from pickle import load
from sys import argv

model=load(open("/home/ewbell/PEPPI/bin/CTNN"))
f=open(argv[1])
point=[float(i) for i in f]
f.close()
posp=model.predict_proba([point])[0][1]
print(posp)
