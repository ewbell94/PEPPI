#!/usr/bin/env python
"""
The pdbchains.py module provides data structures and methods for
handling protein structure and sequence data in PDB or fasta format.

Warning in PDBchains. coord and atom_info.x,atom_info.y and atom_info.z.
share the same reference.  If you set pdb.coord[0] = random_numpy_array.
atom_info values will be different then those in coord.  Which will make
the writing to files incorect.  Only use the provided methods to alter
member data, or you can make a copy as to not alter the class structure.

Inorder to properly change the values in atom info or coord and keep the
references intact. The changes need to happen by element. ie:

pdb.coord[0][0][0] = 37.546
is okay.  Now pdb.atom_info[0][0].coord[0] = 37.546
"""
#python 2.x and 3.x cross compatibility
from __future__ import print_function, division
import sys
if sys.version_info[0] >= 3:
    xrange = range 

import numpy as np
import math
import pandas as pd
import re
import copy
import gzip

def ReadFasta(fin):
    """
    Reads and stores sequence and header information from
    a file in fasta format.

    Arguments:
        fin (str): fasta file
    Returns:
        header_list ([str]): list of headers for each fasta sequence
        sequence_list ([str]): list of sequences.
    """
    header_list = []
    sequence_list = []
    header = ''
    sequence = ''
    f = open(fin)
    for line in f:
        if line[0] == '>':
            header = line.strip()
            break

    for line in f:
        if line[0] == '>':
            header_list.append(header)
            sequence_list.append(sequence)
            sequence = ''
            header = line.strip()
        else:
            sequence += line.strip()
    header_list.append(header)
    sequence_list.append(sequence)
    return (header_list, sequence_list)

def WriteFasta(fout, sequences, headers,resPerLine = 60, append = False):
    """
    Writes a set of sequences and headers out to file in fasta format.

    Requires the list of sequences and headers to be of the same size.

    Arguments:
        fout (str): Name of output file.
        sequences ([str]): List of sequences.
        headers ([str]): List of sequence headers.
        resPerLine (optional, int): Maximum portion of sequence writter per 
            line.  Default write 60 per line.
        append (optional, bool): If true appends data to an already existing
            file.
    """
    write_type = 'w'
    if append:
        write_type = 'a'
    outfile = open(fout, write_type)
    for i in xrange(0, len(headers) ):
        headers[i] = headers[i].strip()
        if headers[i][0] != '>':
            outfile.write('> ') 
        outfile.write(headers[i])
        outfile.write('\n')
        first = 0
        last = resPerLine
        while len(sequences[i]) > first:
            outfile.write( sequences[i][first:last] )
            outfile.write( '\n' )
            first += resPerLine
            last  += resPerLine
    outfile.close()

class PDBstruct:
    """
    This class is a simple structure to contain information in a pdb 
    line or template pdb line.It is used for data storage in the PDBchains
    class.  The data stored can either be a nucleotide or amino acid.

    Attributes:
        atom_num    (int):       Atom number
        atom_name   (string):    Atom name. ie Ca for C-alpha or P 
                        for Phospahte in Nucleotide
        alt_loc     (character): Alternate location idicator
        res_name    (string):    Residue Name is ALA
        chain_id    (character): PDB Chain Identifier
        res_num     (integer):   Residue Number
        res_insert  (character): Insertion of Residue
        coord       ([float,float,float]) Holdes x,y,z coordinates
        coord[0]    (float):    X coordinate
        coord[1]    (float):    Y coordinate
        coord[2]    (float):    Z coordinate
        occup       (double):    Occupancy
        temperature (double):    Temperature Factor
        seg         (character): Segment Identifier
        ele_symb    (string):    Element Symbol
        temp_num    (integer):   Template Residue Number
        temp_res    (string):    Template Residue Name
    """
    def __init__(self,atom_num = None, atom_name=None, alt_loc=None,
                res_name=None, chain_id=None, res_num=None, res_insert=None,
                coord = None, occup=None, temperature=None, 
                seg=None, ele_symb = None, temp_num = '', temp_res = ''):
        self.atom_num    = atom_num
        self.atom_name   = atom_name
        self.alt_loc     = alt_loc
        self.res_name    = res_name
        self.chain_id    = chain_id
        self.res_num     = res_num
        self.res_insert  = res_insert
        self.coord       = coord
        self.occup       = occup
        self.temperature = temperature
        self.seg         = seg 
        self.ele_symb    = ele_symb
        self.temp_num    = temp_num
        self.temp_res    = temp_res


    def __repr__(self):
        return "PDBstruct()"
    def __str__(self):
        return "member of PDBstruct class"

class PDBchains:
    """ 
    This class contains attributes and data methods for handling PDB 
    structures.
    
    Attributes:
        num_chains (int): number of chains currently stored in the object
        fin ([str]): The full file path for the input file.
        fname ([str]): The pdb file name read in for each chain.
        coord (list of N*3 Numpy Array): Each item in the list is the xyz
            coordinates of a pdb chain.  For example coord[0] is a
            N*3 numpy array of the first chain.  coord[0][0] would 
            resturn the 1*3 x y z coordinates of the first atom in
            the first chain.  
        ca_pos (list of N*1 Numpy Array): For protein chains or DNA/RNA.  
            For a protein chain each chain has a list of the indexes 
            corresponding to C alpha atoms in the coord list.
            For DNA/RNA the ca_pos corresponds to the C3'
        atom_info ([[PDBstruct]]): This contains the information for each
            atom in each protein chain.  For example 
            atom_info[chain_num][atom_num] is a given atom in the given
            pdb chain.
        isProtein ([bool]): True if chain protein, False if DNA/RNA.
        codeToabrev (dictionary): Dictionary converts the one amino acid 
            letter code to the three letter abbreviation. 
        abrevTocode  (dictionary): Dictionary converts the three letter 
            amino acid abbreviation to the one letter code. 
            
            AMINO ACIDS Data Structure Standard 20
                Amino Acid          Abbreviation    Code 
                --------------------------------------------
                Alanine         ALA     A
                Arginine        ARG     R
                Asparagine      ASN     N
                Aspartic Acid       ASP     D
                Cysteine        CYS     C
                Glutamine       GLN     Q
                Glutamic Acid       GLU     E
                Glycine         GLY     G
                Histidine       HIS     H
                Isoleucine      ILE     I
                Leucine         LEU     L
                Lysine          LYS     K
                Methionine      MET     M
                Phenylalanine       PHE     F
                Proline         PRO     P
                Serine          SER     S
                Threonine       THR     T
                Tryptophan      TRP     W
                Tyrosine        TYR     Y
                Valine          VAL     V

            Ambiguous Residues
                -------------------------------------------
                                ASP/ASN ambiguous       ASX             B
                                GLU/GLN ambiguous       GLX             Z
                Selenomethione      MSE     M
                Selenocystine       SEC     C 
                Unknown         UNK     X
                

            Deocyribonucleotides/Ribonueclotide Data Structure
            Nucleotide          Abbrevaition    Code
            ----------------------------------------------------
            Adenosine           DA/A        A
            Cytosine            DC/C        C
            Guanine             DG/G        G
            Uracil              U       U
            Thymine             DT      T
            Inosine             DI/I        I
        
        sequence ([str]): Contains the sequence for each protein chain
        init_header (str): Stores header for init data.  Used only for writing
            init file.
        pdbHeader (str): Contains header information for pdb to file.
        template_header ([str]): Stores template header for each template in
            init file.  Used only for writing init file. 

            example: 1h4w pdb has all data fields. 
        ATOM   1312  CA  GLY A 188      59.195  45.264  37.610  1.00 18.27           C
        ATOM   1316  CA  LYS A 188A     58.565  43.487  40.936  1.00 15.22           C
        ATOM   1344  N  AGLN A 192      53.226  33.157  38.087  0.66 12.97           N

            PDB Format (Columns, Data, Jusitification, Data Type)
                1-4   ATOM , string
                7-11  ATOM Number (1), right, integer
                13-16 Atom Name (Ca), left, string
                17    Alternate location indicator (A,B,..Z), character
                18-20 Residue Name (ALA), right, string
                22    Chain Identifier (A), character
                23-26 Residue Sequence Number, right, integer
                27    Code for insertion of residues, character
                31-38 X Coordinate, right, real(8.3) #Three Decimal Places
                39-46 Y Coordinate, right, real(8.3) #Three Decimal Places
                47-54 Z Coordinate, right, real(8.3) #Three Decimal Places
                55-60 Occupancy,    right, real(6.2)
                61-66 Temperature Factor, right, real
                73-76 Segment identifier, left, string
                77-78 Element symbol, right, string
                1-3   TER
            Template PDB Format (For Reading/Writing PDB Templates)
                1-4   ATOM , string
                7-11  ATOM Number (1), right, integer
                13-16 Atom Name (Ca), left, string
                18-20 Residue Name (ALA), right, string
                22    Chain Identifier (A), character
                23-26 Residue Sequence Number, right, integer
                31-38 X Coordinate, right, real(8.3) #Three Decimal Places
                39-46 Y Coordinate, right, real(8.3) #Three Decimal Places
                47-54 Z Coordinate, right, real(8.3) #Three Decimal Places
                56-59 Template Residue Number, right, integer
                61-63 Template Residue Name, string 
                1-3   TER
            init PDB Format
                init header
                template header 1
                template pdb 1
                template header 2
                template pdb 2
                etc. 
    Methods:
        Read:  Used to read a PDB File.  Several inputs allow different 
            information to be stored.
        Write: Used to write a PDB File or a Template PDB File for a 
            selected chain(s).
        MakeModel: Use query to template alignment to generate threading 
            template.
        Transformation: Applies 4x4 Rotation Translation Matrix to 
            Specified Chain.
        Append: Append a copy of a chain or a set of chains from another 
            pdbchains object to this one.
        Extend: Extend one chain with a copy of another.  Usefull for making 
            artificial monomers.  
        Renumber: Used to renumber atom and residue number of selected chain(s)
        ReadFasta:  Reads Fasta File, outputs headers and sequences
        WriteFasta: Takes a list of Sequeces and Headers and Writes 
            Fasta File
        Delete: Delete all attribute data.
        IsStandardAminos: Checks sequence for any non standrad residues
        HasBackbone: Checks if protein chain has backbone atoms
        HasCa:  Checks if data stored is Protein.  Could be DNA or RNA
        Slice: Returns a pdbchains object containg selected chain and data 
               positions.
        SliceResNum: Returns a pdbchains object containing selected chain and 
                     residue numbers
        SliceAtomNum: Returns a pdbchains object containing selected atom 
                      numbers.
        ReferenceLink: Make sure that the numpy coord and atom info
            data structures share the same reference to the 
            coordinates.  The coord data structure is used for
            fast numeric caculations.  Whereas atom_info is used
            for organization of pdb data and reading and writing
            functions.
        ChangeChainID: Change the chain id of a given chain.
        CalculateSideChainCOM: Calculates side chain center of mass
    """
    def __init__(self):
        self.num_chains = 0 
        self.fin  = []
        self.fname = []
        self.coord = []
        self.ca_pos = []
        self.atom_info = []
        self.sequence = []
        self.isProtein = []
        self.abrevTocode = { "GLY": "G", "ALA": "A", "VAL": "V", 
            "LEU": "L", "ILE": "I", "MET": "M", "PHE": "F",
            "TRP": "W", "PRO": "P", "SER": "S", "THR": "T", 
            "CYS": "C", "TYR": "Y", "ASN": "N", "GLN": "Q", 
            "ASP": "D", "GLU": "E", "LYS": "K", "ARG": "R", 
            "HIS": "H", "UNK":"X","ASX":"B","GLX":"Z","MSE":"M",
            "SEC":"U","PYL":"K"}
        self.codeToabrev = {"G": "GLY", "A": "ALA", "V": "VAL",
                        "L": "LEU", "I": "ILE", "M": "MET", "F": "PHE",
                        "W": "TRP", "P": "PRO", "S": "SER", "T": "THR", 
                        "C": "CYS", "Y": "TYR", "N": "ASN", "Q": "GLN", 
                        "D": "ASP", "E": "GLU", "K": "LYS", "R": "ARG",
                        "H": "HIS", "X":"UNK","B":"ASX","Z":"GLX",
                        "U":"SEC"}
        self.standardResidues = {"G": "GLY", "A": "ALA", "V": "VAL",
            "L": "LEU", "I": "ILE", "M": "MET", "F": "PHE",
            "W": "TRP", "P": "PRO", "S": "SER", "T": "THR",
            "C": "CYS", "Y": "TYR", "N": "ASN", "Q": "GLN",
            "D": "ASP", "E": "GLU", "K": "LYS", "R": "ARG",
            "H": "HIS"}
        self.init_header = ''
        self.pdbHeader = ''
        self.template_headers = []

    def __repr__(self):
        return "PDBchains()"
    def __str__(self):
        return "member of PDBchains Class"

    def Read(self,fin, isTemplate = False, isInit = False,
            ca_only = False, backbone_only = False, 
            ter_split = True, allow_alt_loc = True, 
            allow_insert = True,isgziped = False):
        """  
        This method is used to read and store a pdb file or a 
        template/init file (Yang Zhang Format).  The default behavior assumes
        that a pdb file is being read.  Each chain seperated by TER in the
        pdb file is stored in a seperate list here.  The arguments to this
        function fine tunes what information from the file is stored.

        Arguments:
            fin (str):  Name of file to be read. ie 1h4w.pdb
            isTemplate (bool,optional): If set to true reads in a file in 
                template format.  All other arguments are ignored. 
            ca_only (bool,optional): If true only C alpha atoms are stored.
                DNA/RNA chains are ignored.
            backbone_only (bool,optional): If true only backbone atoms are 
                stored. ie. (N, Ca, O).  If ca_only is also true, 
                N and O atoms are skipped.
            ter_split (bool,optional): If true chains are seperated by the 
                TER line in the pdb file. If false, the TER line is ignored 
                and the pdbfile is stored as an artificial monomer. 
            allow_alt_loc (bool, optional): If set to False, only one 
                conformation is allowed per residue.  The first conformation 
                found in the pdb file is the one that is stored.  
                Otherwise if True this information is also stored.
            allow_insert (bool, optional): If set to False all inserted 
                residues are ignored. If True, they are kept.
            isInit: For reading in init.dat files.  Stores regular data
                and header information.  These files contain the same 
                structure as a regular cleaned pdb file.  Except
                the beggining of the file has a header for the whole line.
                And each template has a header before the ATOM data.  Each
                Template is still TER seperated. 
            isgziped: (bool): Default pdb file format is expected to be a text 
                file sometimes downloaded files are gziped compressed.  If this 
                flag is true the program reads the gzipped file.
        """ 
        #store all lines beginning with ATOM or TER
        tmp = fin.split('/')
        filename = tmp[ len(tmp) - 1]
        if isgziped:
            f = gzip.open(fin,'r:gz')
        else:
            f = open(fin)
        pdbin=[line.strip() for line in f.readlines() if line[0:4] == "ATOM" or 
                        line[0:3] == "TER"]
        f.close()

        #count number of chains based on TER.  
        #Store coresponding file name for each chain
        TER_count = 0
        for line in pdbin:
            if line[0:3] == "TER":
                TER_count +=1
        if TER_count == 0: 
            TER_count+=1
        for i in xrange(0,TER_count):
            self.fin.append(fin)
            self.fname.append(filename)
    
        NumChains = self.num_chains
        #Read in pdb or template
        if not isTemplate and not isInit: #read standard pdb
            #if default settings true.  Store Pdb
            if( ca_only==False and backbone_only==False and
                ter_split==True and allow_alt_loc==True and
                allow_insert==True):
                self._ReadPDB(pdbin)
            else:
                #else parse information in accorance with arguments then store 
                #information
                pdbin = self._ParsePDB(pdbin, ca_only, backbone_only, ter_split,
                                         allow_alt_loc, allow_insert)
                self._ReadPDB(pdbin)
        #read template pdbs
        else:
            self._ReadTemplate(pdbin)

        if isInit: #store header information
            f = open(fin)
            headers = [line.rstrip() for line in f if line[0:4] != "ATOM" and 
                            line[0:3] != "TER" and line[0:3] != "END"]
            f.close()
            #check number of headers match number of chains or number chains +1 
            #for when file header included.
            if len(headers) == TER_count:
                headers.insert(0,'')
            elif len(headers) == TER_count + 1:
                self.init_header = headers[0]
            else:
                sys.stderr.write("Number of headers is different than number of"
                        +" chains in file:" + fin +'\n')
                raise IndexError
            for i in xrange(NumChains, NumChains + len(headers)-1):
                self.template_headers[i] = headers[i+1]      

    def Write(self,fout, write_template = False, 
        specific_chains = None, write_details = False, append = False, 
        ter_split = True, write_init = False,write_pdbHeader = False,
        write_connect=False):
        """
        Write writes pdb data stored to PDB chains out to a file

        The Write function has several parameters for writing files
        in pdb or template based format.  The default behavior writes
        a pdb file out with the last piece information being the
        coordinates, ie no occupancy, temperature factor etc.
        The arguments to this function can fine tune what is being
        written

        Arguments:
            fout (str): Name of file to be written
            write_template (bool): If true, writes pdb in template
                format.  Requires that data stored has template
                information.
            write_init (bool): If true, writes file in init format.
                Data must be in template formate.  
                write_template must be set to False 
            specific_chains ([int]): List of integers specifing chains
                to be written.  If None writes all chains.
                specific_chains = [0,1,3] Writes the first, second
                and fourth chain to file.
            write_details (bool): If True, occupancy, temperature
                factor, segment factor and element symbol are
                include in file written.  This data must be
                stored to be written out.
            append (bool): If true appends data to existing file.
            ter_split (bool): If False no TER seperating seperate
                chains in written file.
            write_connect(bool): Append atom connect information to CA.
                trace. Requires that each CA atom has a unique atom number.
                if atom number repeated a flag is raised. If your not sure
                if atoms numbers are unique use the Renumber function.
        """
        #choose writing mode
        write_mode = 'w'
        if append:
            write_mode = 'a'
        outfile = open(fout, write_mode)
        #determine which chains to write to file
        chain_list = []
        if specific_chains is None:
            chain_list = range(0,self.num_chains)
        else:
            chain_list = specific_chains
            for num in chain_list:
                if not isinstance(num,int):
                    err_msg = "Specific chains must be integer list\n"
                    sys.stderr.write(err_msg)
                    raise TypeError 
                if num <0 or num >= self.num_chains:
                    err_msg =''.join(["Integer list must be between 0 and ",
                                "number of ","chains\n"])
                    sys.stderr.write(err_msg)
                    raise IndexError
        #write pdbHeader
        if write_pdbHeader:
            if self.pdbHeader:
                outfile.write(self.pdbHeader)
                if '\n' != self.pdbHeader[len(self.pdbHeader)-1]:
                    outfile.write('\n') 
        #write file
        if write_init:
            if self.init_header:
                outfile.write(self.init_header +'\n')
            write_template = True
        for chain_num in chain_list:
            if write_init:
                if not self.template_headers[chain_num]:
                    sys.stderr.write(''.join([self.fin[chain_num],
                                    " Does not have ","template headers\n"]))
                    raise IndexError
                outfile.write(self.template_headers[chain_num]+'\n')
            for ainfo in self.atom_info[chain_num]:
                outfile.write('%-6s' % "ATOM")         #1-6
                outfile.write('%5s' % ainfo.atom_num)  #7-11
                if len(ainfo.atom_name) <= 3:
                    outfile.write("  ")                    #12 why two spaces?
                    outfile.write('%-3s' % ainfo.atom_name) #13-16    
                else:
                    outfile.write(" ")
                    outfile.write('%-4s' % ainfo.atom_name)
                outfile.write('%1s' % ainfo.alt_loc)           #17
                outfile.write('%3s' % ainfo.res_name)          #18-20
                outfile.write(" ")                     #21 empty space 
                outfile.write('%1s' % ainfo.chain_id)          #22
                outfile.write('%4s' % ainfo.res_num)   #23-26
                outfile.write('%1s' % ainfo.res_insert)              #27
                outfile.write('%3s' % '' )             #28-30 empty space
                outfile.write('%8s' % ("%.3f" % ainfo.coord[0]))         #31-38
                outfile.write('%8s' % ("%.3f" % ainfo.coord[1]))         #39-46
                outfile.write('%8s' % ("%.3f" % ainfo.coord[2]))         #47-54
                if write_template:
                    if ainfo.atom_name != "CA":
                        err_msg = "While writing template atom besides "
                        err_msg += "CA found.  Make sure to read in file "
                        err_msg += "with ca_only set to True"
                        sys.stderr.write(err_msg)
                        raise ValueError
                    outfile.write('%5s' % ainfo.temp_num) #55-59
                    outfile.write('%4s' % ainfo.temp_res) #60-63
                elif write_details:
                    if (ainfo.occup is None or ainfo.temperature is None or
                        ainfo.seg is None or ainfo.ele_symb is None):
                        err_msg = ''.join(["Writing Details Fails for file: ",
                                    self.fin[chain_num],
                                    "\n File has uninitialiezed detail ",
                                     "values check atom "+str(atom_num),
                                    " in chain" +str(chain_num)])
                        sys.stderr.write(err_msg)
                        raise ValueError
                    else:   
                        outfile.write('%6s' % ainfo.occup) #55-60
                        outfile.write('%6s' % ainfo.temperature) #61-66
                        outfile.write('%-4s' % ainfo.seg) #73-76
                        outfile.write('%2s' % ainfo.ele_symb) #77-78
                outfile.write('\n')         

            if ter_split:
                outfile.write("TER\n")
        if not ter_split:
            outfile.write("TER\n")
        if write_connect:
            connectString = self.GenerateConnectData(chain_list)
            outfile.write(connectString)
        outfile.close()


    def GenerateConnectData(self,chainList):
        """
        Creates the connect string for creating artificial bonds between Ca
        atoms in a Calpha trace model.

        Arguments:
            chainList ([int]): Integer list of chains that are used to create
                bond connections. Chains are not connected by a bond.
        """
        uniqueInts = {}
        connectList = []
        for chainI in chainList:
            if len(self.ca_pos[chainI]) < 2:
                continue
            for i in xrange(0,len(self.ca_pos[chainI])-1):
                atomNum1 = self.atom_info[chainI][i].atom_num
                atomNum2 = self.atom_info[chainI][i+1].atom_num
                outStr = ''.join(['CONECT','%5s'%atomNum1,'%5s'%atomNum2,'\n'])
                connectList.append(outStr)
        connectStr = ''.join(connectList)
        return connectStr

    def MakeModel(self, chain_num, query_alignment, template_alignment, 
                    chainID = ''):
        """
        MakeModel creates a template model using the query to template sequence 
        alignment inorder to use the template pdb file coordinates for 
        representing the approximated folded structure of the query sequence.

        Requires:
            The template PDB file must first be read in.  The template_alignment
            and query_alignment must be same length.  They can be obtained from 
            any alignment program.  The expected format of the alignments:
            Example:
                template sequence: GGAAPAKLSR
                query sequence:    GAVAGRS
                After Aligment program.
                template_alignment: GGAAPAK-LSR-
                query_alignment:    -GAV-A-G--RS
            
            Additonaly no alternate locations or insertions should be stored 
            when reading in the PDB.  If alternate locations or insertions are 
            present in the pdb file.  You can prevent them from being stored by 
            setting allow_alt_loc = False and allow_insert = False 
            when calling the Read function.
        Arguments:
            chain_num (int): The chain in PDBchains object being used to 
                generate template model.
            query_alignment (str): query portion of aligment.
            template_alignment (str): template portion of alignment.  This 
                refers to the pdb that has been read in, and will be altered to 
                the template model.
            chainID (optional, char): Chain identification character.   
        Modifies: 
            Selected chain is converted into a template model.  In order to 
            write as a template in the Write function set write_template = True.
            Only CA positions are kept.
        Example:
            pdb = PDBchains()
            #Read in pdb.  No alternate locations or inserts allowed
            pdb.Read('PDB1.pdb',allow_alt_loc = False, allow_insert)
            template_seq = pdb.sequence[0]
            template_align, query_align = NWalign(template_seq,query_seq)
            pdb.MakeModel(0,query_align, template_align)        
        """
        chain = self.atom_info[chain_num]
        ca_pos = self.ca_pos[chain_num]
        sequence = []
        atom_info = [] 
        queryResnum = 0
        tempResnum = 0
        if not chainID:
            chainID = chain[0].chain_id

        coord_list = []
        for i in xrange(0, len(query_alignment) ):
            if query_alignment[i] != '-' and template_alignment[i] != '-':
                canum = ca_pos[ tempResnum ]
                coord_list.append(tempResnum)
                sequence.append(query_alignment[i])
                atom_info.append(
                    PDBstruct( atom_num = chain[canum].atom_num,
                        atom_name = chain[canum].atom_name,
                        alt_loc = '',
                        res_name = self.codeToabrev[query_alignment[i]],
                        chain_id = chainID,
                        res_num = queryResnum + 1,
                        res_insert = '',
                        coord = chain[canum].coord,
                        temp_num = chain[canum].res_num,
                        temp_res = chain[canum].res_name
                    )
                )
                queryResnum += 1
                tempResnum += 1
            elif template_alignment[i] != '-':
                tempResnum += 1
            else:
                queryResnum += 1

        ca_pos = np.array(range(0, len(sequence) ))
        self.atom_info[chain_num] = atom_info
        self.sequence[chain_num] = ''.join( sequence )
        self.ca_pos[chain_num] = ca_pos
        self.coord[chain_num] = self.coord[chain_num][ coord_list, :]
        #Make sure atom_info and coord share same values and reference
        self.ReferenceLink(chains = [chain_num])

    def ReadFasta(self, fin):
        """
        Reads and stores sequence and header information from
        a file in fasta format.

        Arguments:
            fin (str): fasta file
        Returns:
            header_list ([str]): list of headers for each fasta sequence
            sequence_list ([str]): list of sequences.
        """
        return ReadFasta(fin)

    def WriteFasta(self,fout, sequences, headers,resPerLine = 60, 
                    append = False):
        """
        Writes a set of sequences and headers out to file in fasta format.

        Requires the list of sequences and headers to be of the same size.

        Arguments:
            fout (str): Name of output file.
            sequences ([str]): List of sequences.
            headers ([str]): List of sequence headers.
            resPerLine (optional, int): Maximum portion of sequence writter per 
                line.  Default write 60 per line.
            append (optional, bool): If true appends data to an already existing
                file.
        """
        WriteFasta(fout, sequences, headers,resPerLine = resPerLine, 
                append = append)

    def Transformation(self,transformationMatrix,chains = None):
        """
        Applies Transformation Matrix to Protein Chain Coordinates

        This function applies the 4x4 Rotation Translation Matrix
        to a slected protein chain.  The default behavior applies
        the transformation to each chain in the current data
        structure.  Either or both the rotation and translation matricies 
        can be used.
        Default argument skips operation.

        Arguments
            tranformationMatrix (4x4 float numpy array)
            chains ([int]): An interger list of the selected chains to 
                apply the transormation matrix to.  Default
                applies transform to all chains.  
        """
        chain_list = []
        if chains is None:
            chain_list = range(0,self.num_chains)
        else:
            chain_list = chains
        for chain_num in chain_list:
            rotation = transformationMatrix[:3,:3]
            translation = transformationMatrix[:3,3].reshape((1,3))
            tmp_matrix = np.dot(self.coord[chain_num],
                            np.transpose(rotation))+translation
            #copy new coordinates to coord.  I want the reference from
            #atom x,y,z to still have the same values as in coord.
            #This wont happen if I create a new object.  This makes sure
            #Atom info is still referenced to coord and changes accordingly
            atom = 0
            x = 0
            y = 1
            z = 2
            for coordinate in tmp_matrix:
                self.coord[chain_num][atom][x] = coordinate[x]
                self.coord[chain_num][atom][y] = coordinate[y]
                self.coord[chain_num][atom][z] = coordinate[z]
                atom +=1

    def Extend(self, SelfChain, pdbchainsObject, chain, useRenumber = True):
        """
        Extends a pdb chain with another pdbchain.  Usefull for making
        aritificial monomers.

        Arguments:
            SelfChain (int): The chain in this PDBchains object that will
                be extended
            pdbchainsObject (PDBchains object): PDBchains object containing
                chain that will extend a chain in this object.
            chain (int): The chain in pdbchainsObject used for extending
            useRenumber (bool): Renumber artificial monomer from 1 to length of 
                chain.      
        Modifies:
            Extends SelfChain with the pdb data from the chain in pdbchainObject
        """

        ca_pos = pdbchainsObject.ca_pos[chain] + len(self.coord[SelfChain])
        self.ca_pos[SelfChain] =np.concatenate((self.ca_pos[SelfChain], ca_pos))
        self.coord[SelfChain] = np.concatenate((self.coord[SelfChain],
                                    pdbchainsObject.coord[chain]),axis=0)
        self.atom_info[SelfChain].extend(copy.deepcopy(
                                        pdbchainsObject.atom_info[chain]))   
        self.sequence[SelfChain]=self.sequence[SelfChain]+pdbchainsObject.sequence[chain]
        self.ReferenceLink(chains = [SelfChain])
        if useRenumber:
            self.Renumber(chains = [SelfChain])

    def Append(self,pdbchainsObject, chains = None):
        """
        Append a copy of the pdbchains object data from pdbchains object to 
        this object
        
        WARNING: Does not append init_header. Must be done seperately

        Arguments:
            pdbchainsObject (PDBchains object): An initizalized PDBchains
                object.
            chains (optional, [int]): Integer list of chains in pdbchainsObject 
                that will be copied to this object.  The default behavior
                copies all chains.  
        """
        chainList = []
        if chains is None:
            chainList = range(0, pdbchainsObject.num_chains)
        else:
            chainList = chains
        for i in chainList:
            self.num_chains += 1
            self.fin.append(copy.copy(pdbchainsObject.fin[i]))
            self.fname.append(copy.copy(pdbchainsObject.fname[i]))
            self.coord.append(np.copy(pdbchainsObject.coord[i]))
            self.ca_pos.append(copy.copy(pdbchainsObject.ca_pos[i]))
            self.atom_info.append(copy.deepcopy(pdbchainsObject.atom_info[i]))
            self.sequence.append(copy.copy(pdbchainsObject.sequence[i]))
            self.template_headers.append(copy.copy(
                                          pdbchainsObject.template_headers[i]))
            self.isProtein.append(copy.copy(pdbchainsObject.isProtein[i]))  
            #make coord and atom_info share same reference.
            self.ReferenceLink(chains = [self.num_chains-1])

    def Renumber(self, start_res_num = 1, start_atom_num = 1, 
        chains = None, renumber_res = True, renumber_atom = True,
        start_reset = True):
        """
        Renumbers the residue numbers and/or atoms

        Renumbers the residues numbers of a given chain.  
        The default method renumbers each residue and atom in each
        chain starting from 1.
        Arguments
            start_res_num (int): The first number used to start
                renumbering the residue numbers. Default is 1
            start_atom_num(int): The first numbers used to start
                renumbering the atom numbers.  The defaults is
                set to 1.
            renumber_res (bool): If true the residue numbers are
                renumbered.
            renumber_atom (bool): If true the atom numbers are 
                renumbered.
            chains ([int]): Integers list corresponding to chains
                to be renumbred.  Default behavior causes all 
                chains to get renumbered starting from 
                start_res_num.
            start_reset (bool): If True the residue renumber starts
                with start_res_num and start_atom_num for each chain
                Else it starts with the last res_num and atom_num in
                the previous chain          
        """
        chain_list = []
        if chains is None:
            chain_list = range(0,self.num_chains)
        else:
            chain_list = chains
        res_num = start_res_num
        atom_num = start_atom_num

        firstChain = True
        for chain_num in chain_list:
            current_res = self.atom_info[chain_num][0].res_num
            if start_reset:
                res_num = start_res_num
                atom_num = start_atom_num

            else:
                if not firstChain:
                    res_num += 1
            firstChain = False

            for i in xrange(0,len(self.atom_info[chain_num])):
                if self.atom_info[chain_num][i].res_num == current_res:
                    if renumber_res:
                        self.atom_info[chain_num][i].res_num = res_num
                    if renumber_atom:
                        self.atom_info[chain_num][i].atom_num = atom_num
                        atom_num +=1
                else:
                    current_res = self.atom_info[chain_num][i].res_num
                    res_num += 1
                    if renumber_res:
                        self.atom_info[chain_num][i].res_num = res_num
                    if renumber_atom:
                        self.atom_info[chain_num][i].atom_num = atom_num
                        atom_num +=1
    def Delete(self):
        """Delete resets all the member data."""
        self.num_chains = 0
        self.fin  = []
        self.sequence = []
        self.fname = []
        self.coord = []
        self.ca_pos = []
        self.atom_info = []
        self.init_header = ''
        self.template_headers = []
        self.isProtein = []
    def HasBackbone(self,chain_num):
        """
        Checks pdb chain for backbone atoms
        N CA C O
        Arguments:
            chain_num (int): Chain number to check
                for backbone atoms.
        Returns:
            True if all residues have backbone atoms.
            False if missing backbone atoms.    
        """
        if len(self.ca_pos[chain_num]) == 0:
            return False
        AtomsInChain = len(self.atom_info[chain_num])
        atom_info = self.atom_info[chain_num]

        atomDic = {}
        for atom in atom_info:
            if atom.res_num in atomDic:
                atomDic[atom.res_num].append(atom.atom_name)
            else:
                atomDic[atom.res_num] = [atom.atom_name]

        for key in atomDic:
            atomList = atomDic[key]
            if 'N' not in atomList:
                return False
            if 'CA' not in atomList:
                return False
            if 'C' not in atomList:
                return False
            if 'O' not in atomList:
                return False
        return True
        
    def HasCa(self,chain_num):
        """
        Checks to determine if chain has CA.
        
        Used to identify RNA and DNA.  They
        will not have a CA position thus
        this function would return False.

        Arguments
            chain_num (int): The chain to check for
                CA positions.
        Returns True if CA present else False.
        """

        if len(self.ca_pos[chain_num]) == 0:
            return False
        ca_pos = self.ca_pos[chain_num][0]
        if self.atom_info[chain_num][ca_pos] != "CA":
            return False
        else: 
            return True

    def Slice(self, chainNum, positionList):
        """     
        Given an integer list of positions ranging from 0 to the total number of
        atoms stored for the given chain number.  Returns a slice of the 
        pdbchain object containing the data correpsonding to the given chain and
        positions in the list.  Sequence information only stored for CA 
        positions

        Init Header are not copied over.

        Requires: 
            The positionList must be in numerical order ie go from low to high.

        Arguments:
            chain_num (int): The PDB chain to take the slice from.
            positionList ([int]): Integer list of positions to obtain data from 
                the given chain. positionList ranges from 0 to the number of 
                atoms is a chain.
        Returns:
            pdbslice (PDBchains object): A pdbchains object returning the slice 
        """
        if not positionList:
            return PDBchains()

        CaDic = {}
        for i in xrange(0,len(self.ca_pos[chainNum])):
            position = self.ca_pos[chainNum][i]
            CaDic[position] = None

        SeqHits = {}
        sequence = []
        ca_pos_list = []
        positionCount = 0
        #Get sequence slice and new positons for ca_pos
        for i in positionList:
            res_num = self.atom_info[chainNum][i].res_num
            if res_num not in SeqHits:
                SeqHits[res_num] = None
                res_name = self.atom_info[chainNum][i].res_name
                if self.isProtein[chainNum]: #For protein sequence
                    resid = self.abrevTocode[res_name]
                    sequence.append(resid)
                else: #For DNA/RNA sequence
                    resid = res_name[1] if len(res_name) == 2 else res_name[0]
                    sequence.append(resid)
            if i in CaDic:
                ca_pos_list.append(positionCount)

            positionCount += 1

        sequence = ''.join(sequence)
        pdbslice = PDBchains()
        pdbslice.num_chains = 1
        pdbslice.fin.append(self.fin[chainNum])
        pdbslice.fname.append(self.fname[chainNum]) 
        pdbslice.coord.append(  self.coord[chainNum][positionList,:] ) 
        pdbslice.ca_pos.append(np.array(ca_pos_list))
        pdbslice.sequence.append(sequence)
        pdbslice.isProtein.append(self.isProtein[chainNum])
        pdbslice.atom_info.append( [self.atom_info[chainNum][i] for i in positionList] )
        pdbslice.template_headers.append( [self.template_headers[chainNum]] )
        pdbslice.ReferenceLink() 
        return pdbslice

    def SliceResNum(self, chainNum, resNumList,ca_only = False,noRepeat=True):
        """     
        Given an integer list of residue numbers. Returns a slice of the 
        pdbchain object containing the data correpsonding to the given chain and
        residue numbers in the list.

        Notes:
            Init Header are not copied over

        Arguments:
            chain_num (int): The PDB chain to take the slice from.
            resNumList ([int]): Integer list of residue numbers to obtain data 
                from the given chain.
            ca_only (bool): If set to true only CA atoms are returned.
            noRepeat (bool): One one atom type per residue number is allowed if
                set to true
        Returns:
            pdbslice (PDBchains object): A pdbchains object returning the slice 
        """
        lookupTable = {}
        for resnum in resNumList:
            lookupTable[resnum] = None

        positionList = []
        if noRepeat:
            repeats = {}
            for i in xrange(0,len(self.atom_info[chainNum])):
                resNum = self.atom_info[chainNum][i].res_num
                atomName = self.atom_info[chainNum][i].atom_name
                keyRA = (resNum,atomName)
                if resNum in lookupTable and keyRA not in repeats: 
                    positionList.append(i)
                    repeats[keyRA] = None
        else:
            for i in xrange(0,len(self.atom_info[chainNum])):
                if self.atom_info[chainNum][i].res_num in lookupTable:
                    positionList.append(i)

        if ca_only:
            positionTmp = []
            for i in positionList:
                if self.atom_info[chainNum][i].atom_name == 'CA':
                    positionTmp.append(i)
            positionList = positionTmp

        return self.Slice(chainNum, positionList)

    def SliceAtomNum(self, chainNum, atomNumList):
        """     
        Given an integer list of atom numbers. Returns a slice of the pdbchain 
        object containing the data correpsonding to the given chain and residue 
        numbers in the list.

        Notes:
            Init Header are not copied over

        Arguments:
            chain_num (int): The PDB chain to take the slice from.
            atomNumList ([int]): Integer list of residue numbers to obtain data 
                from the given chain.
        Returns:
            pdbslice (PDBchains object): A pdbchains object returning the slice 
        """
        lookupTable = {}
        for atomNum in atomNumList:
            lookupTable[atomNum] = None

        positionList = []
        for i in xrange(0,len(self.atom_info[chainNum])):
            if self.atom_info[chainNum][i].atom_num in lookupTable:
                positionList.append(i)  

        return self.Slice(chainNum, positionList)


    def IsStandardAminos(self,sequence):
        """ 
        Check for non standard residue
        Arguments:
            sequence (string): sequnce string of protein chain
        Returns:
            False, if sequence has non standard residue. True Otherwise
        """
        for res in sequence:
            if res.upper() not in self.standardResidues:
                return False
        return True 

    def ReferenceLink(self,chains = None):
        """
        Gives the values in atom_info.coord the same values as those in 
        self.coord.

        The coord data structure is used for fast numeric caculations.  
        Whereas atom_info is used for organization of pdb data and 
        reading and writing functions.

        Also checks to that atom_info and coord have the same size (ie
        number of atoms stored).

        Arguments:
            chains ([int]): Which chains to restore the reference to.
                If chains = None ie Default behavior,
                reference restored to all chains applied to all 
                chains.

        Modifies:
            Coordinates in atom info and coord will be forced to 
            share the same object.
        """ 
        if chains is None:
            chains = range(0,self.num_chains)
        for chain in chains:
            numberAtoms = self.coord[chain].shape[0]
            if len(self.atom_info[chain]) != numberAtoms:
                err_msg = ''.join(["Atom_Info and Coord have different number ",
                                "of atoms stored\n", "Coord Number of atoms: ", 
                                str(numberAtoms) + '\n',
                                "atom_info number of atoms: ",
                                str(len(self.atom_info[chain]))+'\n'])
                sys.stderr.write(err_msg)
                raise IndexError
            for atomNum in xrange(0, numberAtoms):
                self.atom_info[chain][atomNum].coord =self.coord[chain][atomNum]

    def ChangeChainID(self,chainNumber,chainID):
        """
        Changes the chainID for a given chain.
    
        Arguments:
            chainNumber (int): The chain for which the chain id will be changed
            chainID (char): Character of chain id (one letter string) or empty 
                            string

        """
        chainID = str(chainID)
        if len(chainID) > 1:
            raise ValueError("chainID must be one character or an empty string")
 
        for i in xrange(0,len(self.atom_info[chainNumber]) ):
            self.atom_info[chainNumber][i].chain_id = str(chainID)


    def CalculateSideChainCOM(self,chainNumber):
        """
        Returns a dictionary containing all sidechain center of masses
        for each residue number in the given chain.

        Arguments:
            chainNumber (int): Which chain to extract sidechain COM from
        Returns:
            chain1sgCOM (dic): Each key is a residue number.  Each value
                contains the x,y,z coordinates for the side chain
                center of mass and the residue type ie
                chain1sgCOM[1] = [x,y,z,ResidueCode]
        """
        chain1sg = {}
        countDic = {}
        for i in xrange(0,len(self.atom_info[chainNumber]) ):
            atomName = self.atom_info[chainNumber][i].atom_name
            resNum   = self.atom_info[chainNumber][i].res_num 
            resName  = self.atom_info[chainNumber][i].res_name
            if resNum not in chain1sg:
                chain1sg[resNum] = [0,0,0,self.abrevTocode[resName]]
                countDic[resNum] = 0

            if resName == "GLY" and atomName == 'CA':
                x,y,z = self.atom_info[chainNumber][i].coord
                chain1sg[resNum][0] = x
                chain1sg[resNum][1] = y
                chain1sg[resNum][2] = z
                countDic[resNum]+=1 
    
            elif atomName not in ['N','CA','C','O']:    
                x,y,z = self.atom_info[chainNumber][i].coord
                chain1sg[resNum][0] += x
                chain1sg[resNum][1] += y
                chain1sg[resNum][2] += z
                countDic[resNum]+=1
        chain1sgCOM = {}
        for key in chain1sg:
            count = 1.0*countDic[key]
            if count == 0.0:
                continue
            vals = chain1sg[key]
            xavg = vals[0]/count
            yavg = vals[1]/count
            zavg = vals[2]/count
            resName = vals[3]
            chain1sgCOM[key] = [xavg,yavg,zavg,resName]

        return chain1sgCOM



#Private Helper Functions Below
    def _ReadPDB(self,pdbin):
        """
        Used for pdbin coming from a file in PDB format. Splits pdbin 
        into seperate chains based on TER.  Each ATOM line is then stored
        into the PDB_struct class.  All of the PDBchains attributes are stored
        here or in _ReadTemplate.  Except for init and template headers.  
        """
        pdbchains = self._TERsplit(pdbin)
        if not pdbchains: #error if read stores no data
                    err_msg = ''.join(["_TERsplit does not store pdbfile ",
                                    self.fin[0],
                                    "\n may not contain PDB coordinates ATOM ",
                                    "check file\n"])
                    sys.stderr.write(err_msg)
                    raise IOError   
        self.num_chains = self.num_chains + len(pdbchains)
        for chain in pdbchains:
            atomNum = 0
            ca_pos = []
            atom_info = []
            sequence = []
            xyz_array = np.zeros( (len(chain),3 ) )
            for line in chain:
                if line[13:16].strip() == "CA": #store C alpha locations
                    ca_pos.append(atomNum)
                    res_name = line[17:20].strip()
                    try:
                        sequence.append( self.abrevTocode[res_name] )
                    except:
                        sequence.append('A')
                if line[13:16].strip() == "C3'": #for storing DNA/RNA representative atom
                    ca_pos.append(atomNum)

                xyz_array[atomNum][0] = float( line[30:38] )
                xyz_array[atomNum][1] = float( line[38:46] )
                xyz_array[atomNum][2] = float( line[46:54] )
                if len(line) >= 78:
                    atom_info.append( 
                        PDBstruct( atom_num = int(line[6:11]), 
                            atom_name = line[12:16].strip(), 
                            alt_loc = line[16].strip(), 
                            res_name = line[17:20].strip(), 
                            chain_id = line[21],
                            res_num = int(line[22:26]), 
                            res_insert = line[26], 
                            coord = xyz_array[atomNum],
                            occup = line[54:60],
                            temperature = line[60:66],
                            seg = line[72:76],
                            ele_symb = line[76:78]
                        )
                    )
                else:   
                    atom_info.append( 
                        PDBstruct( atom_num = int(line[6:11]), 
                            atom_name = line[12:16].strip(), 
                            alt_loc = line[16].strip(), 
                            res_name = line[17:20].strip(), 
                            chain_id = line[21],
                            res_num = int(line[22:26]), 
                            res_insert = line[26], 
                            coord = xyz_array[atomNum],
                        )
                    )   
                atomNum +=1
            self.coord.append(xyz_array)
            self.ca_pos.append(np.array(ca_pos))
            self.atom_info.append(atom_info)
            self.template_headers.append('')
            if sequence: #the chain is a protein
                self.isProtein.append( True )
                self.sequence.append( ''.join(sequence) )
            else: #the chain is a nucleotide/store nucleotide sequence
                self.isProtein.append( False )
                res_name = atom_info[0].res_name
                res_code = res_name[1] if len(res_name) == 2 else res_name[0]   
                sequence.append(res_code)
                res_num = atom_info[0].res_num
                for ainfo in atom_info:
                    if res_num != ainfo.res_num:
                        res_name = ainfo.res_name
                        res_code = res_name[1] if len(res_name)==2 else res_name[0]
                        sequence.append(res_code)
                        res_num = ainfo.res_num
                self.sequence.append( ''.join(sequence))
                    
    def _TERsplit(self,pdbin):
        """
        Split pdbin list into list of list each with seperate chain. ie 
        pdbchain[0] is chain1 pdbchain[0][0] is the first line in pdbchain ie 
        ATOM 1 N ASN 1 ... 
        """
        pdbchains = []
        chain = []
        for line in pdbin:
            if line[0:4] == "ATOM":
                chain.append(line)

            elif line[0:3] == "TER":
                if chain: #if chain not empty
                    pdbchains.append(chain)
                    chain = []
            else:
                    err_file = self.inpath + '/' + fin
                    err_msg = "Line that doesnt begin with ATOM or TER\n"
                    err_msg += "Passed filter.  Error file is: " 
                    err_msg += err_file + '\n'
                    sys.stderr.write(err_msg)
                    raise IOError
        if chain: #for pdbs that dont terminate with TER at end
            pdbchains.append(chain)

        return pdbchains        



    #fin passed only to associate each file with chain name
    def _ReadTemplate(self, pdbin):
        """
        Used for pdbin coming from a file in Template or init format
        (yang zhang format). Splits pdbin into seperate chains based on TER.  
        Each ATOM line is then stored into the PDB_struct class.  All of the 
        PDBchains attributes except for init and template headers are 
        created here or in _ReadPDB.  
        """
        pdbchains = self._TERsplit(pdbin)
        if not pdbchains: #error if read stores no data
                    err_msg = "_TERsplit does not store pdbfile" + fin +'\n'
                    sys.stderr.write(err_msg)
                    raise IOError
        self.num_chains = self.num_chains + len(pdbchains)
        for chain in pdbchains:
            atomNum = 0
            ca_pos = []
            atom_info = []
            sequence = []
            xyz_array = np.zeros( (len(chain),3 ) )
            for line in chain:
                if line[13:16].strip() == "CA": #store C alpha locations
                    ca_pos.append(atomNum)
                    res_name = line[17:20].strip()
                    try:
                        sequence.append( self.abrevTocode[res_name] )
                    except:
                        sequence.append('A')
                xyz_array[atomNum][0] = float( line[30:38] )
                xyz_array[atomNum][1] = float( line[38:46] )
                xyz_array[atomNum][2] = float( line[46:54] )
                TNum = -1
                TRes ='XXX'
                try:
                    tnum = int(line[55:59])
                    tres = line[60:63].strip()
                    TNum = tnum
                    TRes = tres
                except:
                    pass
                atom_info.append( 
                    PDBstruct( atom_num = int(line[6:11]), 
                        atom_name = line[12:16].strip(), 
                        alt_loc = line[16].strip(), 
                        res_name = line[17:20].strip(), 
                        chain_id = line[21],
                        res_num = int(line[22:26]), 
                        res_insert = line[26], 
                        coord = xyz_array[atomNum],
                        temp_num = TNum,
                        temp_res = TRes
                        )
                )   
                atomNum +=1
            self.coord.append(xyz_array)
            self.ca_pos.append(np.array(ca_pos))
            self.atom_info.append(atom_info)
            self.sequence.append( ''.join(sequence))
            self.template_headers.append('')
            self.isProtein.append(True)
    def _ParsePDB(self,pdbin, ca_only, backbone_only, ter_split, allow_alt_loc, 
                    allow_insert):
        """Removes data from pdbin list based on argument parameters"""
        exclude = []
        #atom exclusion
        if backbone_only and not ca_only:
            atom_list = ["N","CA", "C", "O"]
            self._AtomParse(pdbin, exclude, atom_list)
        elif ca_only:
            atom_list = ["CA"]
            self._AtomParse(pdbin, exclude, atom_list)
        #if no split on TER
        if not ter_split:
            count = 0
            for line in pdbin:
                if line[0:3] == "TER":
                    exclude.append(count)
                count +=1

        #exclude alternate locations
        if not allow_alt_loc:
            CurrentResNum = -99999.5
            count = 0
            atom_hit = {}
            for line in pdbin:
                if line[0:3] == "TER":
                    count +=1
                    atom_hit = {}
                    CurrentResNum = -99999.5
                    continue

                if line[0:4] == "ATOM":
                    res_num = int(line[22:26])
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    if res_num != CurrentResNum:
                        CurrentResNum = res_num
                        atom_hit = {}

                    if (atom_name,res_name,res_num) in atom_hit:
                        exclude.append(count)
                        count += 1
                    else:
                        atom_hit[atom_name,res_name,res_num] = None
                        count += 1

        #exlused inserts
        if not allow_insert:
            #Column 27 is for insertion of residues.
            count = 0
            for line in pdbin:
                if line[0:3] == "TER" or line == '':
                    pass
                elif line[0:4] == "ATOM":
                    if line[26].strip() != '':
                        exclude.append(count)
                count +=1   


        #make exclude list non redundant and ordered in ascending order
        #slice pdbin according to exclude list
        exclude = set(exclude)
        include = list(set( range(0,len(pdbin)) ) - exclude) 
        include.sort()
        new_pdbin = [pdbin[i] for i in include]         
        return new_pdbin

            
    def _AtomParse(self,pdbin,exclude, atom_list):
        """Searches for atoms not in atom list.  This list is then
        used to prevent those atoms from being stored.  Mainly
        used to read in only calpha positions or backbone positions
        """
        count = 0
        for line in pdbin:
            if line == '' or line[0:3] == "TER":
                pass
            else:
                tmp_atom = line[13:16].strip()
                if tmp_atom not in atom_list:
                    exclude.append(count)
            count +=1

class ReplicaChains:
    """
    Module for reading and writing replica trajectory files from
    I-Tasser

    MODEL
    TER
    ENDMDL

    Need connect.
    dont connect chain split
    """
    def __init__(self):
        self.coord = []
        self.sequence = ''
        self.sequences = []
        self.num_chains = []
        self.num_trajectories = 0
        self.chainEnds = []
        pass

    def ReadTrajectory(self,inFile):
        pass

    def WriteTrajectory(self,outFile):
        pass
