#!/usr/bin/env python
# -*- coding: utf-8 -*-
doc='Calculate TMscore No suposition'

import os
import sys
import pdbchains
import numpy as np
import re
import getpass
import subprocess
from subprocess import Popen, PIPE
from pdbchains import PDBchains

bindir=os.path.dirname(os.path.realpath(__file__))
def TMscore_rTmscore(modelfile,nativefile,outputdirectory):
    """
    input file: pdb file
    """
    #read input file
    native = PDBchains()
    model = PDBchains()
    native.Read(nativefile,ca_only = True)
    model.Read(modelfile,ca_only = True)
    # TMscore,RMSD,rTMscore,TMscor1,TMscore2
    modelA,modelB,nativeA,nativeB,TMscore,RMSD=Calculate_Complex_TMscore_superpose(modelfile,nativefile,outputdirectory)
    TMscoreA=CalculateTMscore(modelA,nativeA)

    TMscoreB=CalculateTMscore(modelB,nativeB)
    rTMscore=2.0/(1.0/TMscoreA +1.0/TMscoreB)

    ##Calulate LRMSD
    LRMSD = Calulate_LRMSD(model,native)
    
    ##Calulate IRMSD,fnat,acc
    I_RMSD,fnat,acc,f1score=Get_Interface_RMSD(modelfile,nativefile)
    return TMscore,RMSD,TMscoreA,TMscoreB,rTMscore,LRMSD,I_RMSD,fnat,acc,f1score

def CalculateTMscore(model,native):
    
    """
    Calculates the TMscore between the modelCoords and nativeCoords.
    The TMscore was obtained in the paper: Scoring Function for 
    Automated Assessment of Protein Structure Template Quality.
    Yang Zhang and Jeffrey Skolnick
    TMscore = (1/Length) *aligned_sum(1/1+(di/d0)^2)
    d0 = 1.24*(nativeLength - 15)^(1/3) - 1.8 
    """
    abrevTocode = { "GLY": "G", "ALA": "A", "VAL": "V","LEU": "L",
                    "ILE": "I", "MET": "M", "PHE": "F","TRP": "W",
                    "PRO": "P", "SER": "S", "THR": "T","CYS": "C",
                    "TYR": "Y", "ASN": "N", "GLN": "Q","ASP": "D",
                    "GLU": "E", "LYS": "K", "ARG": "R","HIS": "H",
                    "UNK": "X", "ASX": "B", "GLX": "Z","MSE": "M",
                    "SEC":"U","PYL":"K"}
    x = []
    y = []
    z = []
    # model coordinate
    x1 = []
    y1 = []
    z1 = []
    
    n1=[]
    modelLength=len(model.sequence[0])
    nativeLength=len(native.sequence[0])

    for index in xrange(0,modelLength):
        for index1 in xrange(0,nativeLength):
            resNum = model.atom_info[0][index].res_num                    # residue number
            resNum1 = native.atom_info[0][index1].res_num
            if resNum == resNum1:
               if model.atom_info[0][index].res_name == native.atom_info[0][index1].res_name:
                  model_coord = model.atom_info[0][index].coord
                  native_coord = native.atom_info[0][index1].coord
                  x.append(native_coord[0])
                  y.append(native_coord[1])
                  z.append(native_coord[2])
                  x1.append(model_coord[0])
                  y1.append(model_coord[1])
                  z1.append(model_coord[2])
                  n1.append(resNum1)

    
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    x1 = np.array(x1)
    y1 = np.array(y1)
    z1 = np.array(z1)
    if nativeLength > 15:
       d0=1.24*((nativeLength - 15)**(1.0/3)) - 1.8
       if d0 < 0.5:
          d0=0.5
    else:
       d0=0.5

    di=np.sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)#distance
    tmpScore = np.sum(1/(1+(di/d0)**2))
    TMscore=tmpScore/nativeLength
    return TMscore

def Calculate_Complex_TMscore_superpose(modelfile,nativefile,outputdirectory):
    """
    Calculate complex TMscore

    input: modelfile,nativefile
    outputfile: monomer modelA/modelB; monomer nativeA/nativeB; Superpose file
    """
    native = PDBchains()
    model = PDBchains()
    native.Read(nativefile,ca_only = True)
    model.Read(modelfile,ca_only = True)
    #*************Calcute complex TMscore*******
    if model.num_chains ==2 and native.num_chains==2:
       cmd=[bindir+'/TMscore','-c',modelfile,nativefile]
       p=subprocess.Popen(cmd,stdin=PIPE, stdout=PIPE)
       stdout, stderr=p.communicate()
    else:
        
       cmd=['TMscore',modelfile,nativefile]
       p=subprocess.Popen(cmd,stdin=PIPE, stdout=PIPE)
       stdout, stderr=p.communicate()

    RMSD_pattern=re.compile("RMSD of  the common residues=")
    TMscore_pattern=re.compile("TM-score    =")
    delespace=re.compile('\s+')
    start=0
    """
    transformationMatrix (4x4 Numpy Array): 4x4Transformation matrix
    ----------------------------------------------------------------
    u    - u(i,j) is   rotation  matrix 
    t    - t(i)   is translation vector
    
    for example:
    -------- rotation matrix to rotate Chain-1 to Chain-2 ------
     i          t(i)         u(i,1)         u(i,2)         u(i,3)
     1     28.6668792618   0.2314398185   0.9407221819  -0.2479463387
     2     24.3935993589  -0.4747761081   0.3316732967   0.8152180514
     3     -7.5158130253   0.8491308837  -0.0709549201   0.5233948239
     
     rotating Chain_1 from (x,y,z) to (X,Y,Z)
    
     X(i)=t(1)+u(1,1)*x(i)+u(1,2)*y(i)+u(1,3)*z(i)
     Y(i)=t(2)+u(2,1)*x(i)+u(2,2)*y(i)+u(2,3)*z(i)
     Z(i)=t(3)+u(3,1)*x(i)+u(3,2)*y(i)+u(3,3)*z(i)

    """
    transformationMatrix=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    for line in stdout.splitlines():
        #print line
        if RMSD_pattern.findall(line):
           RMSD=float(delespace.split(line)[-1])
        if TMscore_pattern.findall(line):
           TMscore=float(delespace.split(line)[2])
        if re.findall("u\(i\,1\)",line):
            start=1
            i=0
        elif line.startswith('Superposition'):
             start=0
        if len(line)>1 and start==1:
           matrix=delespace.split(line)
           transformationMatrix[i][0]=matrix[2]
           transformationMatrix[i][1]=matrix[3]
           transformationMatrix[i][2]=matrix[4]
           transformationMatrix[i][3]=matrix[5]
           i +=1
    rotationMatrix=np.array(transformationMatrix)[1:4,1:4]       # u
    rotationMatrix=rotationMatrix.astype('float64')
    transformationVector=np.array(transformationMatrix)[1:4,0]   # t
    transformationVector=transformationVector.astype('float64')

    if model.num_chains==2:
       model.coord[0]=np.dot(model.coord[0],rotationMatrix.T)+transformationVector
       model.coord[1]=np.dot(model.coord[1],rotationMatrix.T)+transformationVector
    
       # generate template Superpose file
       modelA=PDBchains()
       modelA.Append(model,chains = [0] )
       modelC=modelA
       modelB=PDBchains()
       modelB.Append(model,chains = [1] )
       modelC.Append(model,chains=[1])
       modelC.Write(outputdirectory+'/'+'Superpose.pdb')
    
    elif model.num_chains==1:
       modelA=model
       modelB=model
       modelC=modelA
       modelC.Append(model,chains=[0])
       modelC.Write(outputdirectory+'/'+'Superpose.pdb')

    
       # generate native monomer A and B
    if native.num_chains==2:
       nativeA=PDBchains()
       nativeB=PDBchains()
       nativeA.Append(native,chains=[0])
       nativeB.Append(native,chains=[1])

    elif native.num_chains==1:
         nativeA=model
         nativeB=model


    return modelA,modelB,nativeA,nativeB,TMscore,RMSD


def Calulate_LRMSD(model,native):
    """
    L-RMSD us the ligand RMSD calculated after optimal 
    superposition of the receptor structure to the native.
    receptor: the larger protein chain
    ligand: the small protein chain
    """
    if os.getenv('SLURM_JOBID'):
       workdir='/scratch/'+getpass.getuser()+'/'+os.getenv('SLURM_JOBID')
    else:
        workdir='/tmp/'+getpass.getuser()
    if not os.path.exists(workdir):
       os.mkdir(workdir)
   # native = PDBchains()
    #model = PDBchains()
    #native.Read(nativefile,ca_only = True)
    #model.Read(modelfile,ca_only = True)
    nativeALength=native.sequence[0]
    nativeBLength=native.sequence[1]
    modelALength=model.sequence[0]
    modelBLength=model.sequence[1]
    
    ReceptorNative = PDBchains()
    LigandNative = PDBchains()
    ReceptorModel = PDBchains()
    LigandModel = PDBchains()

    ##  temp file
    ReceptorModelfile=workdir+'/'+'ReceptorModel.pdb'
    LigandModelfile=workdir+'/'+'LigandModel.pdb'

    ReceptorNativefile=workdir+'/'+'ReceptorNative.pdb'
    LigandNativefile=workdir+'/'+'LigandNative.pdb'
    ######################################################
    if nativeALength <= nativeBLength:
       ReceptorModel.Append(model,chains=[0])
       LigandModel.Append(model,chains=[1])
       ReceptorNative.Append(native,chains=[0])
       LigandNative.Append(native,chains=[1])
       ReceptorModel.Write(ReceptorModelfile)
       LigandModel.Write(LigandModelfile)
       ReceptorNative.Write(ReceptorNativefile)
       LigandNative.Write(LigandNativefile)

    else:
       ReceptorNative.Append(native,chains=[1])
       LigandNative.Append(native,chains=[0])

       ReceptorModel.Append(model,chains=[1])
       LigandModel.Append(model,chains=[0])

       ReceptorModel.Write(ReceptorModelfile)
       LigandModel.Write(LigandModelfile)
       ReceptorNative.Write(ReceptorNativefile)
       LigandNative.Write(LigandNativefile)

    
    cmd = [bindir+'/RMSD',ReceptorModelfile,ReceptorNativefile]
    p=subprocess.Popen(cmd,stdin=PIPE, stdout=PIPE)
    stdout, stderr=p.communicate()
    ReceptorMatrix=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    delespace=re.compile('\s+')
    start=0
    for line in stdout.splitlines():
        #print line
        if re.findall("u\(i\,1\)",line):
           start=1
           i=0
        elif line.startswith('Superposition'):
           start=0
        if len(line)>1 and start==1:
           matrix=delespace.split(line)
           ReceptorMatrix[i][0]=matrix[2]
           ReceptorMatrix[i][1]=matrix[3]
           ReceptorMatrix[i][2]=matrix[4]
           ReceptorMatrix[i][3]=matrix[5]
           i +=1
    ReceptorRotationMatrix=np.array(ReceptorMatrix)[1:4,1:4]       # u
    ReceptorRotationMatrix=ReceptorRotationMatrix.astype('float64')
    ReceptorVector=np.array(ReceptorMatrix)[1:4,0]         # t
    ReceptorVector=ReceptorVector.astype('float64')

    LigandModel.coord[0]=np.dot(LigandModel.coord[0],ReceptorRotationMatrix.T)+ReceptorVector
    
    Ligand_model_Supose=PDBchains()
    Ligand_model_Supose.Append(LigandModel,chains=[0])
    Ligand_model_Supose.Write(workdir+'/'+'Ligand_model_Supose.pdb')
    LRMSD=Calculate_RMSD(Ligand_model_Supose,LigandNative)
    cmd='rm -rf '+workdir
    os.system(cmd)
    return LRMSD




def Calculate_RMSD(model,native):
    """
    Calculate RMSD no superpose

    RMSD= sqrt(sum(di**2)/L)

    di: distance
    L: 
    """
   # native = PDBchains()
    #model = PDBchains()
    #native.Read(nativefile,ca_only = True)
    #model.Read(modelfile,ca_only = True)


    x = []
    y = []
    z = []
    # model coordinate
    x1 = []
    y1 = []
    z1 = []
    n1=[]

    modelLength=len(model.sequence[0])
    nativeLength=len(native.sequence[0])
    for index in xrange(0,modelLength):
        for index1 in xrange(0,nativeLength):
            resNum = model.atom_info[0][index].res_num                    # residue number
            resNum1 = native.atom_info[0][index1].res_num
            if resNum == resNum1:
               model_coord = model.atom_info[0][index].coord
               native_coord = native.atom_info[0][index1].coord
               x.append(native_coord[0])
               y.append(native_coord[1])
               z.append(native_coord[2])
               x1.append(model_coord[0])
               y1.append(model_coord[1])
               z1.append(model_coord[2])
               n1.append(resNum1)
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    x1 = np.array(x1)
    y1 = np.array(y1)
    z1 = np.array(z1)
    #print len(x),x[0],y[0],z[0]
    #print len(x1),x1[0],y1[0],z1[0]
    di=np.sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)#distance
    RMSD=np.sqrt(np.sum(di**2)/modelLength)
    return RMSD

    
    
def Get_Interface_RMSD(modelfile,nativefile,discutoff=10): 
    """
    calulate the interface of the chainA and chainB 
    input:  modelfile
            native file

    Zi Liu
    07/04/2018
    """
    #################################################################################
    # Interface RMSD use program RMSD.f
    # RMSD.f download from https://zhanglab.ccmb.med.umich.edu/TM-score/
    # Ref. Y. Zhang, J. Skolnick, Scoring function for automated assessment of protein
    # structure template quality, Proteins, 57: 702-710 (2004)
    ##################################################################################

    if os.getenv('SLURM_JOBID'):
       workdir='/scratch/'+getpass.getuser()+'/'+os.getenv('SLURM_JOBID')
    else:
       workdir='/tmp/'+getpass.getuser()
    if not os.path.exists(workdir):
       os.mkdir(workdir)

    ####################################
    # step 1. read native PDB file
    #         Calculate native (chainA and chainB interface)
    ####################################
    tmppdb=PDBchains()
    tmppdb.Read(nativefile,ca_only = True)
    pdbchainB=PDBchains()
    pdbchainA=PDBchains()
    pdbchainA.Append(tmppdb,chains=[0])
    pdbchainA.Renumber(start_res_num = 1,chains=[0])
    pdbchainB.Append(tmppdb,chains=[1])
    pdbchainB.Renumber(start_res_num = len(pdbchainA.sequence[0])+1,chains=[0])

    # for chainA coordinate
    ca_posA   = pdbchainA.ca_pos[0]
    chainAlen = len(ca_posA)
    CAcoordA = pdbchainA.coord[0][ca_posA,:]

    # for chainB coordinate
    ca_posB   = pdbchainB.ca_pos[0]
    chainBlen = len(ca_posB)
    CAcoordB = pdbchainB.coord[0][ca_posB,:]
    numberCa = CAcoordA.shape[0]  # row
    #######################################################################
    # Calculate native interface/contacet parts 

    ######################################################################
    distance_list=[]  # for row is chainA position,colum is chainB position
    contact_list=[]
    for i in xrange(numberCa):
        distance=np.sqrt(np.sum( (CAcoordB - CAcoordA[i])**2, axis = 1))
        isContact = (distance<=discutoff)
        distance_list.append(distance)
        contact_list.append(isContact)
    distanceMatrix = np.array(distance_list)
    contactMatrix = np.array(contact_list)
    # for chainA interface/contact residues number
    ChainAInterfaceCa=[i for i in xrange(0,chainAlen) if True in contactMatrix[i,:] ]
    # for chainB interface/contact residues number
    ChainBInterfaceCa=[i for i in xrange(0,chainBlen) if True in contactMatrix[:,i] ]
    numContacts = np.sum(contactMatrix) # the number of interface/contact
    #print numContacts
   # row=contactMatrix.shape[0]
    #colum=contactMatrix.shape[1]
    #if not ChainAInterfaceCa or not ChainBInterfaceCa:
       #print("No interface between the two chains,please check the distance cutoff")
       #exit()
    ## ChainA/chainB interface/contact residues,coordinate
    #chainAresnum = [] # chainA interface/contact residues number
    #chainBresnum = [] # chainB interface/contact residues number
    
    #for i in ChainAInterfaceCa:
        #ca_pos = pdbchainA.ca_pos[0][i]
        #resnum = pdbchainA.atom_info[0][ca_pos].res_num
        #chainAresnum.append(resnum)
    
    #for i in ChainBInterfaceCa:
        #ca_pos = pdbchainB.ca_pos[0][i]
        #resnum = pdbchainB.atom_info[0][ca_pos].res_num
        #chainBresnum.append(resnum)

    #NativeInterface=pdbchainA.SliceResNum(0, chainAresnum)
    #NativeInterface.Append(pdbchainB.SliceResNum(0,chainBresnum))
    ## generate native interface file
    #nativeInfile=workdir+'/NativeInterface.pdb'
    #NativeInterface.Write(nativeInfile,ter_split=False)

    ######################################################################
    ## for model interface 
    ######################################################################
    tmpmodel=PDBchains()
    tmpmodel.Read(modelfile,ca_only = True)
    modelchainA=PDBchains()
    modelchainB=PDBchains()
    modelchainA.Append(tmpmodel,chains=[0])
    modelchainB.Append(tmpmodel,chains=[1])
    modelchainArenum=[] # model chainA interface residue number
    modelchainBrenum=[] # model chainB interface residue number
    # for model chainA
    #for i in modelchainA.ca_pos[0]:
        #res_num = modelchainA.atom_info[0][i].res_num
        #for j in chainAresnum:
            #if int(j)==int(res_num):
               #modelchainArenum.append(res_num)
    ## for model chainB
    #for i in modelchainB.ca_pos[0]:
        #res_num = modelchainB.atom_info[0][i].res_num
        #for j in chainBresnum:
            #if int(j)==int(res_num):
               #modelchainBrenum.append(res_num) 
    #ModelInterface=modelchainA.SliceResNum(0, modelchainArenum)
    #ModelInterface.Append(modelchainB.SliceResNum(0,modelchainBrenum))
    
    ## generate model interface file
    #modelInfile=workdir+'/ModelInterface.pdb'
    #ModelInterface.Write(modelInfile,ter_split=False)
    ########################################################################
    ## Calculate RMSD
    #######################################################################
    #IRMSD_pattern=re.compile("RMSD of the common residues=")
    #cmd=['RMSD',modelInfile,nativeInfile]
    #p=subprocess.Popen(cmd,stdin=PIPE,stdout=PIPE)
    #stdout, stderr=p.communicate()
    #for line in stdout.splitlines():
        #if IRMSD_pattern.match(line):
           #IRMSD=line.split()[-1]
           #print 'IRMSD: ',IRMSD
    
    ######################################################################
    # Calculate fnat and ACC
    # fraction of native contacts (fnat)
    # fnat is the number of antive (correct) residue-residue contacts 
    # in the predicted complex divided by the number of contacts in the target complex
    # Accuracy of the predicted interface contacts
    # ACC is the fraction of correct predictions divided by the total number of predictions.
    ######################################################################
    
    #########################################################################
    # Calculate model interface parts
    #########################################################################
    
    # for model chainA coordinate
    m_ca_posA   = modelchainA.ca_pos[0]
    m_chainAlen = len(m_ca_posA)
    m_CAcoordA = modelchainA.coord[0][m_ca_posA,:]
    # for model chainB coordinate
    m_ca_posB   = modelchainB.ca_pos[0]
    m_chainBlen = len(m_ca_posB)
    m_CAcoordB = modelchainB.coord[0][m_ca_posB,:]
    
    m_numberCaA = m_CAcoordA.shape[0] # row
    dis_list=[]  # for row is model chainA position,colum is model chainB position
    cont_list=[] # 
    for i in xrange(m_numberCaA):
        distance1=np.sqrt(np.sum( (m_CAcoordB - m_CAcoordA[i])**2, axis = 1))
        isContact1 = (distance1<=discutoff)
        dis_list.append(distance1)
        cont_list.append(isContact1)
    distMatrix = np.array(dis_list)
    contMatrix = np.array(cont_list)
    numConts_model = np.sum(contMatrix) # the number of model interface/contact
    ############################################
    # calculate residues-residues interface pairs
    # residues-residues number
    #############################################
    if numContacts >=3:
       row=contactMatrix.shape[0]
       colum=contactMatrix.shape[1]
       # ChainA/chainB interface/contact residues,coordinate
       chainAresnum = [] # chainA interface/contact residues number
       chainBresnum = [] # chainB interface/contact residues number
    
       for i in ChainAInterfaceCa:
           ca_pos = pdbchainA.ca_pos[0][i]
           resnum = pdbchainA.atom_info[0][ca_pos].res_num
           chainAresnum.append(resnum)
    
       for i in ChainBInterfaceCa:
           ca_pos = pdbchainB.ca_pos[0][i]
           resnum = pdbchainB.atom_info[0][ca_pos].res_num
           chainBresnum.append(resnum)

       NativeInterface=pdbchainA.SliceResNum(0, chainAresnum)
       NativeInterface.Append(pdbchainB.SliceResNum(0,chainBresnum))
       # generate native interface file
       nativeInfile=workdir+'/NativeInterface.pdb'
       NativeInterface.Write(nativeInfile,ter_split=False)

       ######################################################################
       ## for model interface 
       ######################################################################
       #tmpmodel=PDBchains()
       #tmpmodel.Read(modelfile,ca_only = True)
       #modelchainA=PDBchains()
       #modelchainB=PDBchains()
       #modelchainA.Append(tmpmodel,chains=[0])
       #modelchainB.Append(tmpmodel,chains=[1])
       #modelchainArenum=[] # model chainA interface residue number
       #modelchainBrenum=[] # model chainB interface residue number
       # for model chainA
       for i in modelchainA.ca_pos[0]:
           res_num = modelchainA.atom_info[0][i].res_num
           for j in chainAresnum:
               if int(j)==int(res_num):
                  modelchainArenum.append(res_num)
       # for model chainB
       for i in modelchainB.ca_pos[0]:
           res_num = modelchainB.atom_info[0][i].res_num
           for j in chainBresnum:
               if int(j)==int(res_num):
                  modelchainBrenum.append(res_num) 
       ModelInterface=modelchainA.SliceResNum(0, modelchainArenum)
       ModelInterface.Append(modelchainB.SliceResNum(0,modelchainBrenum))
       # generate model interface file
       modelInfile=workdir+'/ModelInterface.pdb'
       ModelInterface.Write(modelInfile,ter_split=False)
       #######################################################################
       # Calculate RMSD
       ######################################################################
       if len(modelchainArenum)+len(modelchainBrenum)!=0:
          IRMSD_pattern=re.compile("RMSD of the common residues=")
          cmd=[bindir+'/RMSD',modelInfile,nativeInfile]
          p=subprocess.Popen(cmd,stdin=PIPE,stdout=PIPE)
          stdout, stderr=p.communicate()
          for line in stdout.splitlines():
              #print line
              if IRMSD_pattern.match(line):
                 IRMSD=line.split()[-1]
                 IRMSD=float(IRMSD)
                 #print 'IRMSD: ',IRMSD
       else:
          IRMSD=-1
          #print IRMSD
    else:
      
        IRMSD=-1  # native no interface 
        #print IRMSD
    ############################################# 
    native_interface_list=[] 
    for i in xrange(0,chainAlen):
        for j in xrange(0,chainBlen):
            if contactMatrix[i,j]:
               resnumA= pdbchainA.atom_info[0][i].res_num
               resnumB = pdbchainB.atom_info[0][j].res_num
               native_interface_list.append([resnumA,resnumB])
    model_interface_list=[]
    for i in xrange(0,m_chainAlen):
        for j in xrange(0,m_chainBlen):
            if contMatrix[i,j]:
               modelresnumA=modelchainA.atom_info[0][i].res_num
               modelresnumB=modelchainB.atom_info[0][j].res_num
               model_interface_list.append([modelresnumA,modelresnumB])
    correct_count=0
    for i in native_interface_list:
        for j in model_interface_list:
            if i==j:
               #print i,j
               correct_count+=1
    if numContacts >0 and numConts_model >0:   # if the number of natve/model(prediction) interface is below 3 
       fnat=float(correct_count)/numContacts  # recall
       acc=float(correct_count)/numConts_model # precision
       if fnat ==0 and acc==0:
          f1score=0
       else:
          f1score=2*((fnat*acc)/(fnat+acc)) 
    else:
       fnat=0
       acc=0
       f1score=0
    #print IRMSD,fnat,acc,f1score
    cmd='rm -rf '+workdir
    os.system(cmd)

    
    


    return IRMSD,fnat,acc,f1score





if __name__=="__main__":

   modelfile=sys.argv[1]
   nativefile=sys.argv[2]
   outputdirectory=sys.argv[3]
   #Get_Interface_RMSD(modelfile,nativefile)
   #RMSD=Calculate_RMSD(modelfile,nativefile)   
   #print RMSD
   TMscore,RMSD,TMscoreA,TMscoreB,rTMscore,LRMSD,IRMSD,fnat,acc,f1score=TMscore_rTmscore(modelfile,nativefile,outputdirectory)
   #LRMSD=Calulate_LRMSD(modelfile,nativefile)
   print 'TM-score: ',TMscore, '\nRMSD: ',RMSD,'\nTM-scoreA: ',TMscoreA,'\nTM-scoreB: ',TMscoreB,'\nrTM-score: ',rTMscore,'\nLRMSD: ',LRMSD, \
   '\nIRMSD',IRMSD,'\nfnat',fnat,'\nacc',acc,'\nf1score',f1score
   #print TMscore,RMSD,TMscoreA,TMscoreB,rTMscore,LRMSD

    


    



