#!/usr/bin/env python
"""
The alignment.py module provides a set of functions for sequence and 
structure alignment.

NWalign: Aligns to protein sequences using blosum scoring matrix
TMalign: Aligns to protein structures.
TMscore: Aligns a model to native structure given an alignment.
KabschFunction: Uses protein sequence alignment for structure superposition and
    accesing structural similarity.
CalculateTMscore: Evaluates the TMscore function score.
KabschTransformation: Performs the kabsch rotation with translation
KabschRotation: Performs kabsch rotation algorithm.
CalculateRMSD: Calculates the Root Mean Square Deviation between 
    two sets of points

Many of the functions returns several items.  They are stored and
returned as a Pandas Series which is a structured data array.  
A quick google search will inform of how to use.
"""
#python 2.x and 3.x cross compatibility
from __future__ import print_function, division
import sys
if sys.version_info[0] >= 3:
    xrange = range

import os
ROOT_PATH = os.path.dirname(os.path.realpath(__file__) ) + '/' 
sys.path.append(ROOT_PATH+"/NWalign/")
sys.path.append(ROOT_PATH+"/TMalign/")
sys.path.append(ROOT_PATH+"/TMscore/")
import nwalign                  #fortran wrapper located in PDB/NWalign
import tmalign                  #function wrapper located in PDB/TMalign 
import TMscoreMod               #function wrapper located in PDB/TMscore
import numpy as np
import pandas as pd
from pdbchains import PDBchains


def NWalign(sequenceA, sequenceB,sequenceAIsFile=False,sequenceBIsFile=False):
    """
    NWalign functions aligns two proteins sequences (seqA and seqB).
    
    Arguments:
        sequenceA (str): Protein sequence of first chain or fasta
            file if sequenceAIsFile = True.
        sequenceB (str): Protein sequence of second chain or fasta
            file if sequenceBIsFile = True.
        sequenceAIsFile (optional,bool): If true. sequenceA is expected
            to be the name of a fasta file.
        sequenceBIsFile (optional,bool): If true. sequenceB is expected
            to be the name of a fasta file.

    Returns (Pandas series consisting of data below):
        seqidA (float): Sequence identity for first chain normalized by length
        seqidB (float): Sequence identity for second chain normalized by length
        seqid  (float): Sequence identity of the alignment
        numAlign (int): Number of aligned residues
        seqA_align (str): Alignment for first sequence
        align (str): string of ' ' and ':' representing unalinged
            and aligned residues
        seqB_align (str): Alignment of second sequence
    Example:
        out = NWalign('GAGG','AAGG')
        out = NWalign('1abcA.fasta', 'AAGG', sequenceAIsFile = True)
        out['seqidA'] == out[0]
        out['numAlign'] == out[2]
    """
    if sequenceAIsFile:
        pdbTmp = PDBchains()
        (header,sequences) = pdbTmp.ReadFasta(sequenceA)
        sequenceA = sequences[0]

    if sequenceBIsFile:
        pdbTmp = PDBchains()
        (header, sequences) = pdbTmp.ReadFasta(sequenceB)
        sequenceB = sequences[0]

    out = nwalign.compares(sequenceA, sequenceB)
    lengthMax = max(len(sequenceA),len(sequenceB) )*1.0
    seqidA = out[0]
    seqidB = out[1]
    seqid = max(seqidA,seqidB)
    numAlign = out[2]
    seqA_align = ''.join(out[3])
    align = ''.join(out[4])
    seqB_align = ''.join(out[5])
    output = pd.Series([seqidA,seqidB,seqid,numAlign,seqA_align,align,
                        seqB_align],
        index=["seqidA", "seqidB","seqid", "numAlign", "seqA_align","align",
                "seqB_align"])
    return output

def TMalign(model, native, modelIsFile = False, nativeIsFile = False, 
            modelchain = 0, nativechain = 0,applyTransform = False): 
    """
    TMalign does structural superposition for two protein chains.

    Arguments:
        model (file/PDBchains):  If modelIsFile = True.  The model input should 
            be the full name and path of a pdb file.  Else model is expected to 
            use the PDBchains data structure.
        native (file/PDBchains):  If nativeIsFile = True.  The native input 
            should be the full name and path of a pdb file.  Else model is 
            expected to use the PDBchains data structure.
        modelIsFile (bool,optional):  Default is False.  If true model should 
            be the full path of a pdb file.  Else model should be stored using 
            PDBchains structure
        nativeIsFile (bool,optional): Default is False.  If true native should 
            be the full path of a pdb file.  Else native should be stored using 
            PDBchains structure.
        modelchain (int, optional):  The chain in the pdb file or PDBchains 
            that should be used. Default is the first chain. 
            Second chain = 1, Third = 2 ., etc.
        nativechain (int, optional): The chain in the pdb file or PDBchains 
            that should be used. Default is the first chain.
        applyTransform (bool, optional): Apply the rotation and transformation 
            operations on the model PDB.
    Returns (transformationMatrix, pandasSeries):
        transformationMatrix (4x4 Numpy Array): 4x4 Transformation matrix.   
        seqid (float): Sequence Identity of alignment. Normalized by number 
                        aligned.
        numAlign (int): number of aligned residues
        rmsd (double): Root mean squared deviation score
        TMscore_normModel (double): TMscore normalized by length of model
        TMscore_normNative (double): TMscore normalized by length of native
        seqModel_align (string): alignment of model sequence 
        alignment (str): string of ' ' and ':' representing unaligned and 
                        aligned regions.
        seqNative_align (string): alignemnt of native seqeunce.

        Aligment example
        seqNative_alignment: G-ARP-VSGA
        alignment:           : : : : ::
        seqModel_alignment:  GGA-PKL-PA
    Example:
        x = PDBchains()
        x.Read('pdb1.pdb')
        (transformationMatrix,TMout) = TMalign(x,'pdb2.pdb',nativeIsFile = True)
        TMout['seqid'] == TMout[0]
        TMout['rmsd'] == TMout[2]
        TMout['seqNative_alignment'] == TMout[7]
    """
    seq = []
    CAcoord = []
    data_tmp = [ [model, modelIsFile, modelchain] , 
                [native, nativeIsFile, nativechain] ]
    for tmp in data_tmp:
        isfile = tmp[1]
        if isfile:
            fileName = tmp[0]
            chainNum = tmp[2]       
            pdb = PDBchains()
            pdb.Read( fileName,ca_only = True, allow_alt_loc = False, 
                        allow_insert = False )
            seq.append( pdb.sequence[chainNum ] )
            CAcoord.append( pdb.coord[chainNum] )
        else:
            chainNum = tmp[2]
            pdb = tmp[0]
            seq.append( pdb.sequence[chainNum] )
            CAcoord.append( pdb.coord[chainNum][pdb.ca_pos[chainNum],: ] )  

    #out = tmalign._tmalign(seqB, CAcoordB, seqA, CAcoordA)
    #chainA is considered nativie
    out  = tmalign._tmalign(seq[0],CAcoord[0], seq[1],CAcoord[1])
    transformationMatrix = out[0]
    seqid = out[1]
    alignment = out[2]
    numAlign = out[3]
    rmsd = out[4]
    seqModel_align = out[5]
    TMscore_normModel = out[6]
    seqNative_align = out[7]
    TMscore_normNative = out[8]

    if applyTransform:
        model.Transformation(transformationMatrix, chains = [modelchain])

    output = pd.Series([seqid,numAlign,rmsd,TMscore_normModel,
                TMscore_normNative,seqModel_align,alignment,seqNative_align],
        index=["seqid", "numAlign", "rmsd", "TMscore_normModel",
        "TMscore_normNative","seqModel_align","alignment","seqNative_align"])
    return transformationMatrix, output

def TMscore(model, native, modelchain = 0, nativechain = 0, modelIsFile = False,
    nativeIsFile = False, useNWalign = False, applyTransform=False):
    """
    The TMscore uses the given alignment inorder to determine coordinates
    to use in the Kabsch rotation function along with returning the TMscore.  
    If either the modelAlignment or nativeAlignment are missing, they are 
    calculated in this function using NWalign.  The rotation translation 
    matrices returned superimpose the model to the native. 

    Requires:
        model (PDBchains object/file): PDB file read in using the pdbchains 
            module. Or the file name of a PDB file when modelIsFile is set 
            to True
        modelchain (opt,int): The default chain used is the first one read in.
            by PDBchains.Read.  However any chain that is stored can be used
            for the calculation by changing this value.
        native (PDBchains object/file): PDB file read in using the pdbchains 
            module.Or the file name of a PDB file when nativeIsFile is 
            set to True
                nativechain (optional,int): The default chain used is the first 
                        one read in. by PDBchains.Read.  However any chain that 
                        is stored can be used for the calculation by changing 
                        this value.
        modelIsFile (optional, bool): If model is string containing file name.  
            This should be set to true.  Default expects model to be PDBchains 
            object.
        nativeIsFile (optional,bool): If native is string containing file name. 
            This should be set to true.  Default expects native to be PDBchains 
            object.
        useNWalign (optional, bool): Generate initial alignment using NWalign 
            before running TMscore.  Default is false and common residue numbers
            are used for the alignment.  The native residue numbers should go 
            from 1 to L
        applyTransform (bool, optional): Apply the rotation and transformation 
            operations on the model PDB.
    Returns (transformationMatrix, panda series):
        transformationMatrix (4x4 Numpy Array): 4x4 Transformation matrix.
        TMscoreModel  (float): TMscore normalized by length of model
        TMscoreNative (float): TMscore normalized by length of native
        rmsd (float): Root mean squared deviation between model and native.
        numAlign (int): Number of aligned residues from optimal TMscore 
                superposition.
    Example:
        x = PDBchains()
        x.Read('pdb1.pdb')
        (transformationMatrix, TMout) = TMalign('pdb2.pdb',x ,modelIsFile=True)
        TMout['TMscoreModel'] == TMout[0]
        TMout['rmsd'] == TMout[2]

    """
    if modelIsFile:
        fileName = model
        model = PDBchains()
        model.Read(fileName)
    if nativeIsFile:
        fileName = native
        native = PDBchains()
        native.Read(fileName)

    modelLen = len(model.sequence[modelchain])
    nativeLen = len(native.sequence[nativechain])
    if useNWalign:
        modelSequence  = model.sequence[modelchain]
        nativeSequence = native.sequence[nativechain]
        (seqIDmodel,seqIDnative,seqid,numAlign,modelAlignment,align,nativeAlignment) = NWalign(modelSequence,nativeSequence)
        tmpnative = PDBchains()
        tmpnative.Append(native, [nativechain] )
        tmpnative.Renumber()
        native = tmpnative
        tmpmodel = PDBchains()
        tmpmodel.Append(model, [modelchain] )
        tmpmodel.MakeModel(0, nativeAlignment, modelAlignment)
        model = tmpmodel    
    
    x1 = []
    y1 = []
    z1 = []
    n1 = []
    for i in model.ca_pos[modelchain]:
        coordinates = model.atom_info[modelchain][i].coord
        res_num = model.atom_info[modelchain][i].res_num
        x1.append(coordinates[0])
        y1.append(coordinates[1])
        z1.append(coordinates[2])
        n1.append(res_num)
    if len(n1) != len(set(n1)):
        print("TMscore error.  Repeated residue number in model")
        repeatsDic = {}
        for i in set(n1):
            repeatsDic[i] = 0
        for i in n1:
            repeatsDic[i] += 1
        repeats = []
        for key in repeatsDic:
            if repeatsDic[key] > 1:
                repeats.append(key)
        print("repeated residue numbers: ",repeats)
        raise ValueError
    x1 = np.array(x1)
    y1 = np.array(y1)
    z1 = np.array(z1)
    n1 = np.array(n1)

    x2 = []
    y2 = []
    z2 = []
    n2 = []
    for i in native.ca_pos[nativechain]:
        coordinates = native.atom_info[nativechain][i].coord
        res_num = native.atom_info[nativechain][i].res_num
        x2.append(coordinates[0])
        y2.append(coordinates[1])
        z2.append(coordinates[2])
        n2.append(res_num)

    if len(n2) != len(set(n2)):
        print("TMscore error.  Repeated residue number in native")
        repeatsDic = {}
        for i in set(n2):
            repeatsDic[i] = 0
        for i in n2:
            repeatsDic[i] += 1
        repeats = []
        for key in repeatsDic:
            if repeatsDic[key] > 1:
                repeats.append(key)
        print("repeated residue numbers: ",repeats)
        raise ValueError    
    
    x2 = np.array(x2)
    y2 = np.array(y2)
    z2 = np.array(z2)
    n2 = np.array(n2)

    TM,rmsd,numAlign,u,t = TMscoreMod.tmscore(x1,y1,z1,n1,x2,y2,z2,n2)  
    transformationMatrix = np.identity(4)
    for i in xrange(0,3):
        transformationMatrix[i][3] = t[i]
        for j in xrange(0,3):
            transformationMatrix[i][j] = u[i][j]
    TMscoreNative = TM
    TMscoreModel = (TM*nativeLen)/(1.0*modelLen)
    output = pd.Series([TMscoreModel, TMscoreNative, rmsd,numAlign],
                index=["TMscoreModel","TMscoreNative","rmsd","numAlign"])

    if applyTransform:
        model.Transformation(transformationMatrix, chains = [modelchain])

    
    return transformationMatrix, output
    
    
def KabschFunction(model, native, modelchain = 0, nativechain = 0, 
    modelIsFile = False, nativeIsFile = False, modelAlignment = '',
    nativeAlignment = '',applyTransform = False):
    """
    The Kabsch Function uses the given alignment inorder to determine 
    coordinates to use in the Kabsch rotation function along with returning the 
    TMscore.  If either the modelAlignment or nativeAlignment are missing, they 
    are calculated in this function using NWalign.  The rotation translation 
    matrices returned superimpose the model to the native. 

    Requires:
        model (PDBchains object/file): PDB file read in using the pdbchains 
            module. Or the file name of a PDB file when modelIsFile is 
            set to True
        modelchain (opt,int): The default chain used is the first one read in.
            by PDBchains.Read.  However any chain that is stored can be used
            for the calculation by changing this value.
        native (PDBchains object/file): PDB file read in using the pdbchains 
            module.Or the file name of a PDB file when nativeIsFile is 
            set to True
                nativechain (optional,int): The default chain used is the first 
                        one read in. by PDBchains.Read.  However any chain that 
                        is stored can be used for the calculation by changing 
                        this value.
        modelAlignment  (optional,str): The model portion of the sequence 
                        alignment.  
        nativeAlignment (optional,str): The native portion of the sequence 
                        alignment.
            Aligment example
                modelAlignment : G-ARP-VSGA
                alignment      : : : : : ::
                nativeAlignment: GGA-PKL-PA
        modelIsFile (optional, bool): If model is string containing file name.  
            This should be set to true.  Default expects model to be PDBchains 
            object.
        nativeIsFile (optional,bool): If native is string containing file name. 
            This should be set to true.  Default expects native to be PDBchains 
            object.
        applyTransform (optional,bool): Apply the rotation and translation 
            matricies to the model chain.  Default is False. 
    Returns (transformationMatrix, panda series):
        transformationMatrix (4x4 Numpy Array): 4x4 Transformation matrix.
        TMscoreModel  (float): TMscore normalized by length of model
        TMscoreNative (float): TMscore normalized by length of native
        rmsd (float): Root mean squared deviation between model and native.
    Example:
        x = PDBchains()
        x.Read('pdb1.pdb')
        (transformationMatrix, TMout) = TMalign('pdb2.pdb',x ,modelIsFile=True)
        TMout['TMscoreModel'] == TMout[0]
        TMout['rmsd'] == TMout[2]

    """
    if modelIsFile:
        fileName = model
        model = PDBchains()
        model.Read(fileName)
    if nativeIsFile:
        fileName = native
        native = PDBchains()
        native.Read(fileName)
    
    #if either has empty string obtain alignment using NWalign
    if not modelAlignment or not nativeAlignment:
        modelSequence  = model.sequence[modelchain]
        nativeSequence = native.sequence[nativechain]
        (seqIDmodel,seqIDnative,seqid,numAlign,modelAlignment,align,nativeAlignment) = NWalign(modelSequence,nativeSequence)
    modelAlignPos  = [] #residue positions in model where there is alignment to native
    nativeAlignPos = [] #residue positions in native where there is alignment to model

    mcount = 0 #model ca place holder
    ncount = 0 #native ca place holder
    for i in xrange(0,len(modelAlignment)):
        if modelAlignment[i] != '-' and nativeAlignment[i] != '-':  
            modelAlignPos.append(mcount)
            nativeAlignPos.append(ncount)
            mcount += 1
            ncount += 1

        elif modelAlignment[i] != '-':
            mcount += 1
        else:
            ncount += 1

    #obtain CA coordinates of aligment residues for model and native pdbs
    modelCApos = model.ca_pos[ modelchain ]
    CaAligned = [modelCApos[i] for i in modelAlignPos]
    modelCaAligned  = model.coord[ modelchain ][CaAligned,:] 

    nativeCApos = native.ca_pos[ nativechain ]
    CaAligned = [nativeCApos[i] for i in nativeAlignPos]
    nativeCaAligned = native.coord[nativechain][CaAligned,:]

    (transformationMatrix, rmsd) = KabschTransformation(modelCaAligned,
                                            nativeCaAligned)

    #apply rotation translation to modelCaAligned for calculating TMscore
    rotation = transformationMatrix[:3,:3]
    translation = transformationMatrix[:3,3].reshape((1,3))
    modelCaAligned = np.dot(modelCaAligned, np.transpose(rotation)) +translation

    (TMscoreModel,TMscoreNative) = CalculateTMscore( modelCaAligned, 
                len(model.sequence[modelchain]),
                nativeCaAligned, len(native.sequence[nativechain])   )      
        
    if applyTransform:
        model.Transformation(transformationMatrix, chains = [modelchain] )
    output = pd.Series([TMscoreModel, TMscoreNative, rmsd],
                    index=["TMscoreModel","TMscoreNative","rmsd"])
    return (transformationMatrix, output)


def CalculateAlignmentScores(model, native, modelchain = 0, nativechain = 0, 
        modelIsFile = False, nativeIsFile = False, modelAlignment = '',
        nativeAlignment = ''):
    """
    Calculates the TMscore and RMSD of model and native without any 
    transformations. Uses the current position of model and native for the 
    alignment score calculations 

    Requires:
        model (PDBchains object/file): PDB file read in using the pdbchains 
            module. Or the file name of a PDB file when modelIsFile is 
            set to True
        modelchain (opt,int): The default chain used is the first one read in.
            by PDBchains.Read.  However any chain that is stored can be used
            for the calculation by changing this value.
        native (PDBchains object/file): PDB file read in using the pdbchains 
            module. Or the file name of a PDB file when nativeIsFile is 
            set to True
                nativechain (optional,int): The default chain used is the first 
                        one read in. by PDBchains.Read.  However any chain that 
                        is stored can be used for the calculation by changing 
                        this value.
        modelAlignment  (optional,str): The model portion of the sequence 
                        alignment.  
        nativeAlignment (optional,str): The native portion of the sequence 
                        alignment.
            Aligment example
                modelAlignment : G-ARP-VSGA
                alignment      : : : : : ::
                nativeAlignment: GGA-PKL-PA
        modelIsFile (optional, bool): If model is string containing file name.  
            This should be set to true.  Default expects model to be 
            PDBchains object.
        nativeIsFile (optional,bool): If native is string containing file name. 
            This should be set to true.  Default expects native to be PDBchains 
            object.
        applyTransform (optional,bool): Apply the rotation and translation 
            matricies to the model chain.  Default is False. 
    Returns panda series:
        TMscoreModel  (float): TMscore normalized by length of model
        TMscoreNative (float): TMscore normalized by length of native
        rmsd (float): Root mean squared deviation between model and native.
    Example:
        x = PDBchains()
        x.Read('pdb1.pdb')
        (transformationMatrix, TMout) = TMalign('pdb2.pdb',x ,modelIsFile=True)
        TMout['TMscoreModel'] == TMout[0]
        TMout['rmsd'] == TMout[2]

    """
    if modelIsFile:
        fileName = model
        model = PDBchains()
        model.Read(fileName)
    if nativeIsFile:
        fileName = native
        native = PDBchains()
        native.Read(fileName)
    
    #if either has empty string obtain alignment using NWalign
    if not modelAlignment or not nativeAlignment:
        modelSequence  = model.sequence[modelchain]
        nativeSequence = native.sequence[nativechain]
        (seqIDmodel,seqIDnative,seqid,numAlign,modelAlignment,align,nativeAlignment) = NWalign(modelSequence,nativeSequence)
    modelAlignPos  = [] #residue positions in model where there is alignment to native
    nativeAlignPos = [] #residue positions in native where there is alignment to model

    mcount = 0 #model ca place holder
    ncount = 0 #native ca place holder
    for i in xrange(0,len(modelAlignment)):
        if modelAlignment[i] != '-' and nativeAlignment[i] != '-':  
            modelAlignPos.append(mcount)
            nativeAlignPos.append(ncount)
            mcount += 1
            ncount += 1

        elif modelAlignment[i] != '-':
            mcount += 1
        else:
            ncount += 1

    #obtain CA coordinates of aligment residues for model and native pdbs
    modelCApos = model.ca_pos[ modelchain ]
    CaAligned = [modelCApos[i] for i in modelAlignPos]
    modelCaAligned  = model.coord[ modelchain ][CaAligned,:] 

    nativeCApos = native.ca_pos[ nativechain ]
    CaAligned = [nativeCApos[i] for i in nativeAlignPos]
    nativeCaAligned = native.coord[nativechain][CaAligned,:]

    rmsd = CalculateRMSD(modelCaAligned,nativeCaAligned)
    (TMscoreModel,TMscoreNative) = CalculateTMscore( modelCaAligned, 
                                    len(model.sequence[modelchain]),
        nativeCaAligned, len(native.sequence[nativechain])   )      

    output = pd.Series([TMscoreModel, TMscoreNative, rmsd],
                        index=["TMscoreModel","TMscoreNative","rmsd"])
    return output



def CalculateTMscore(modelCoords,modelLength, nativeCoords, nativeLength, 
        distanceCutoff = -1,d0=-1):
    """
    Calculates the TMscore between the modelCoords and nativeCoords.
    The TMscore was obtained in the paper: Scoring Function for 
    Automated Assessment of Protein Structure Template Quality.
    Yang Zhang and Jeffrey Skolnick

    TMscore = (1/Length) *aligned_sum(1/1+(di/d0)^2)
    d0 = 1.24*(nativeLength - 15)^(1/3) - 1.8 

    Requires:
        modelCoords and nativeCoords (model and native coordinates)
        are N*3 numpy arrays with same shape.  Additionaly the 
        modelCoords should be superimposed onto the nativeCoords.
        The points in model and native coords should represent the
        one to one alignment.
    Arguments:
        modelCoords  (N*3 numpy array): The calpha coordinates of a
            PDB file that are aligned to the native.
        modelLength (int): Length of the Model chain.
        nativeCoords (N*3 numpy array): The calpha coordinates of a
            PDB file that are aligned to the model.
        nativeLength (int): Length of the Native chain.
        distanceCutoff (opt,float): Distance required for TMscore alignment
            calculation between residues.  If set to negative, a function
            dependant on native length is used to determine the 
            distanceCutoff else the user defined value is used.
        d0 (opt,float): d0 parameter to calculate TMscore
    Returns (TMscoreModel, TMscoreNative):
        TMscoreModel  (float): TMscore normalized by the model length.
        TMscoreNative (float): TMscore normalized by the native length.
    """
    if len( modelCoords ) == 0 or len(nativeCoords) == 0:
        return 0

    if d0 <= 0:
        if nativeLength > 15:
            d0 = 1.24*((nativeLength - 15)**(1.0/3)) - 1.8
            if d0 < 0.5:
                d0 = 0.5
        else:
            d0 = 0.5
 
    #calculate distance between each aligned residue    
    squaredDistance = np.sum( (modelCoords-nativeCoords)**2, axis = 1)

    #remove residues distances greater than residue cutoff
    if distanceCutoff < 0:
        distanceCutoff = 1.5*(nativeLength**0.3) + 3.5  
    squaredCutoff = distanceCutoff**2
    distance = np.sqrt(squaredDistance[ squaredDistance <= squaredCutoff ]  )

    tmpScore=  np.sum(1/(1+(distance/d0)**2))
    TMscoreModel = tmpScore/modelLength
    TMscoreNative = tmpScore/nativeLength
    return (TMscoreModel, TMscoreNative)

def CalculateRMSD(X,Y,distanceCutoff = -1):
    """
    Calculates Root Mean Squared Deviation between point sets X and Y.
    
    Requires:
        X and Y are both numpy arrays with the same shape.
    Arguments:
        X ( N*D numpy array ): N is the number of rows.  D is the
            number of dimensions.
        Y ( N*D numpy array ): N is the number of rows.  D is the
            number of dimenstions.
        distanceCutoff (opt, float):  Only conisder points whose
            distance is below the cutoff.  
        Each row is a point in euclidean space.
    Returns (rmsd).
        rmsd (float): Root mean squared deviation.
    """
    squaredDistance = np.sum( (X-Y)**2, axis = 1 )
    #remove residue with distances greater than distanceCutoff
    if distanceCutoff > 0:
        squaredDistance = squaredDistance[ squaredDistance < distanceCutoff**2 ]
    
    rmsd = np.sqrt( np.sum(squaredDistance)/len(squaredDistance) ) 
    return rmsd 

def KabschTransformation(X,Y):
    """
    The Kabsch algorithm: for optimal rotation translation transformation matrix
    between a set of points.  The algorithm starts with two sets of points X and
    Y centers them on the origin and returns the optimal rotation and 
    translation matrix for superposition of X onto Y.  

    Requires:
        X and Y are expected to be N*D Numpy arrays where N is the 
        number of rows and D is the dimension.  The X and Y matricies
        should be of the same size.   

    Arguments:
        X ( N*D numpy array )
        Y ( N*D numpy array )
    Returns (transformationMatrix, rmsd):
        transformationMatrix ( D+1 * D+1 numpy array ): Rotation 
            Translation Transformation Matrix.  For 3*N 
            coordinates returns 4X4 rotation translation 
            matrix.
        rmsd ( float ): Root Mean Square Deviation

    Algorithm Explanation:

    The kabsch algorithm finds the optimal rotation between two sets of
    coordinates with at least 3 data points each.  Below I present the
    derivation of the kabsch alogirthm used to find the optimal rotation
    matrix between two sets of points inorder to minimize the rmsd
    (root mean squared deviation).

    The kabsch algorithm seeks to find the rotation matrix U that minimzes 
    the funciton:

    E = sum( |xn -U*yn|^2 ) where xn and yn are a set of 3D coordinates.
    The function E can be expanded to:

    E = sum(|xn|^2 + |yn|^2) - 2*sum( xn*U*yn).  Which can be rewritten
    as E = E0 - 2*sum( xn*U*yn).  Notice E0 does not depend on U.  In
    order to minimize E we need to maximize 2*sum( xn*U*yn).

    The trick is in noticing sum( xn*U*yn) = Trace(X*U*Y^t) where X and Y 
    are Nx3 matricies where the columns are the x,y and z coordinates and
    each row is a seperate data point, and Y^t is the transpose of Y.
    Using matrix notation are goal is to maximize L where
    L = Trace(X*U*Y^t)

    Using the cyclic property of the trace notice that
    L =Trace(X*U*Y^t)   'trace of NxN matrix'
      = Trace(U*Y^t*X)  'trace of 3x3 matrix'
      = Trace(U*R).

    U is still unknown but R is known and its the correlation matrix
    calculated as R=Y^t*X.  Using Singular Value Decomposition R can
    be decomposed into three matrices.  R = V*S*W^t. V and W are
    orthogonal matrices and S is a diaganol containing the singular
    values oi.

    Using SVD we rewrite the Trace(U*R).
    Trace(U*R) = Trace(U*V*S*W^t) = Trace(S*W^t*U*V).  We can group
    the matrices W^t*U*V into the matrix T.  Now we have
    Trace(S*T) = o1*T11 + o2*T22 + o3*T33.  Knowing that T is the
    product of orthoganal matrices each value of Tii <= 1.  When
    Tii = 1 and T = I (identity matrix) the Trace(S*T) is maximized
    Thus Tmax = W^t*U*V = I

    U = W*V^t.

    The kabsch algorithm can return two possible answers a rotation
    matrix (right handed coordinate system) and a reflectin of
    a rotation matrix (left handed coordinate system).  If the
    determinant of U is 1 its right handed and the problem is solved.
    If the determinant of U is -1 the solution is left handed.  In
    order to get a right handed system the third column of U needs
    to be multiplied by -1.
    """
    
    #Place x,y on origin by subtracting center of mass (com).
    xcom = np.mean(X, axis = 0)
    ycom = np.mean(Y, axis = 0)
    x = X - xcom
    y = Y - ycom    
 
    #Calculate Optimal rotation matrix
    U = KabschRotation(x,y)
    #coordinate dimension + 1
    Dimen = x.shape[1] + 1

    #create D+1xD+1 transformation matrix.
    #translation matrix for X to origin
    XtoOrigin = np.identity(Dimen)
    #translation matrix for X on origin to Y center of mass
    OriginToY = np.identity(Dimen)
    Rmatrix = np.identity(Dimen) #rotation matrix around origin
    for i in xrange(0,Dimen-1):
        XtoOrigin[i][Dimen-1] = -1*xcom[i]
        OriginToY[i][Dimen-1] = ycom[i]
        for j in xrange(0,Dimen-1):
            Rmatrix[i][j] = U[i][j]

    transformationMatrix = np.dot(OriginToY, np.dot(Rmatrix,XtoOrigin) )    
    #(U*x.T).T = (x*U.T) where .T == Transpose.
    x = np.dot(x,np.transpose(U))
    rmsd = CalculateRMSD(x,y)
    return (transformationMatrix, rmsd)

def KabschRotation(x,y):
    """
    Returns optimal rotation matrix U.

    Requires x and y are same shape and centered on origin.

    Arguments:
        x and y ( N*D numpy array): Both have to be same shape and
            centered on origin.

    Returns:
        U (D*D numpy array): Optimal rotation matrix for x superposition
            onto y.
    """
    C = np.dot(np.transpose(x), y)
    V, S, Wt = np.linalg.svd(C)
    #U is optimal rotation about origin
    U = np.dot(np.transpose(Wt),np.transpose(V))
    #if U left handed make right handed
    if np.linalg.det(U) < 0: #U left handed make right handed
        Wt[2,:] *= -1
        U = np.dot(np.transpose(Wt),np.transpose(V))
    return U
