#!/home/ewbell/miniconda2/bin/python

import os;
import sys;
import commands;
import numpy


##### usage 
# ./TMsearch.py REF_model1.pdb out.fasta

####### set variables, some of you may change #######
librarydir="/oasis/projects/nsf/mia181/zhanglab/library"
pdbdir="/home/ewbell/SPRINGDB/monomers/";
tmalign=librarydir+"/contact/NEW/DMPfold/TMalign";
tmaligncpp=librarydir+"/contact/NEW/DMPfold/TMaligncpp";
#list_path=pdbdir+'/list';
#tmalign="/home/zhengwei/bin/TMalign"
#list_path='/home/zhengwei/amino-zhengwei/LOMETS/DMPmod/list2'

query=sys.argv[1]
out_path=sys.argv[2]

listtag='all'
list_path=pdbdir+'/list';

if len(sys.argv)>=4:
    listtag=sys.argv[3]

'''
if listtag=="domain":
    list_path="/nfs/amino-home/zhng/local_library/HHD_summary/domain.list"
elif listtag=="domain60":
    list_path="/nfs/amino-home/zhng/local_library/HHD_summary/domain60.list"
elif listtag=="domain50":
    list_path="/nfs/amino-home/zhng/local_library/HHD_summary/domain50.list"
elif listtag=="domain40":
    list_path="/nfs/amino-home/zhng/local_library/HHD_summary/domain40.list"
else:
    list_path=pdbdir+'/list';
'''
list_path="/home/ewbell/SPRINGDB/70CDHITstruct.mono"
print 'list-tag:',listtag
print 'Now searching list is',list_path

########
#    >template_name TMscore_nz_by_query TMscore_nz_by_template
#    XXXXX-XX----XXXX   # query
#    XXX-XXXXXXXX--X-  # template
#    sorted by TMscore_nz_by_query
#    using TMaligncpp do fast search (worse accuracy), then use TMalign do full align with Top 300 template
########
topN=5000    ### do full TMalign with Top300 templates
#get list file and remove duplicate items

os.system(tmaligncpp+" -dir1 "+pdbdir+" "+list_path+" -suffix .pdb "+query+" -outfmt 2 -fast >fast_search.tab")

tmpoutfile=file("tmpids.list",'w')
tmpoutfile.write(str(topN)+"\n")
tems=[]
tms=[]
tmpinfile=file("fast_search.tab",'r')
tmpinfile.readline()
lines=tmpinfile.readlines()
tmpinfile.close()
for line in lines:
#    print line
    if not line.startswith("Total"):
        tmpstrs=line.strip('\n').split('\t')
        tems.append(tmpstrs[0].split(".")[0])
        tms.append(float(tmpstrs[3]))

tmqs=numpy.array(tms)
indexlist0 = numpy.argsort(-tmqs);

cc=0
for idd in indexlist0:
    cc+=1
    tmpoutfile.write(tems[idd]+"\n");
    if cc>=topN:
        break

tmpoutfile.close()

idlist=[]

idfile=file("tmpids.list",'r')
idlines=idfile.readlines()[1:]
idlist_tmp=[]
for idline in idlines:
    idlist_tmp.append(idline.strip('\n'));

idlist=list(dict.fromkeys(idlist_tmp))

print 'total enerties:',len(idlist)

tmscoreq_list=[]
tmscoret_list=[]
alignments=[]

tmscoreq=-1
tmscoret=-1
alignment='-'
tmpfile=file(out_path+'.unsort','w')

for id in idlist:
    #### do TM-align with input structure
    (s,o)=commands.getstatusoutput(tmalign+' '+pdbdir+'/'+id+'.pdb '+query)
    #print o
    if o.__contains__("if normalized by length of Chain_2"):
        strs=o.split('\n');
        for stri in range(0,len(strs)):
            strtmp=strs[stri];
            if strtmp.__contains__("if normalized by length of Chain_1"):
                tmscoret=float(strtmp.split(' ')[1])
            if strtmp.__contains__("if normalized by length of Chain_2"):
                tmscoreq=float(strtmp.split(' ')[1])
            if strtmp.__contains__(" denotes aligned residue pairs of d < 5.0 A"):
                break;
        alignment=strs[stri+1]+'\n';
        alignment+=strs[stri+3]+'\n';

    if tmscoreq!=-1 and tmscoret!=-1 and alignment!='X':
        tmscoreq_list.append(tmscoreq)
        tmscoret_list.append(tmscoret)
        alignments.append(alignment)
        tmpfile.write('>'+id+' '+str(tmscoreq)+' '+str(tmscoret)+'\n')
        tmpfile.write(alignment);

    ###### reset #####
    tmscoreq=-1
    tmscoret=-1
    alignment='-'

tmpfile.close()

#print tmscoreq_list
#print tmscoret_list
#print alignments

tmqlist=numpy.array(tmscoreq_list)
indexlist = numpy.argsort(-tmqlist);

outfile=file(out_path,'w')

for index in indexlist:
    outfile.write(idlist[index]+' '+str(tmscoreq_list[index])+' '+str(tmscoret_list[index])+'\n')
    #outfile.write(alignments[index])

outfile.close()
