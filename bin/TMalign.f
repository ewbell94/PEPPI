
**************************************************************************
*     This program is to identify the best alignment of two protein 
*     structures that gives the highest TM-score. Input structures must 
*     be in the PDB format. By default, TM-score is normalized by the 
*     second protein. Users can obtain a brief instruction by simply 
*     running the program without arguments. For comments/suggestions,
*     please contact email: zhng@umich.edu.
*     
*     Reference to cite:
*     Yang Zhang, Jeffrey Skolnick, Nucl. Acid Res. 2005 33: 2302-9
*     
*     Permission to use, copy, modify, and distribute this program for 
*     any purpose, with or without fee, is hereby granted, provided that
*     the notices on the head, the reference information, and this
*     copyright notice appear in all copies or substantial portions of 
*     the Software. It is provided "as is" without express or implied 
*     warranty.
************************ updating history ********************************
*     2005/06/01: A bug of two-point superposition was fixed.
*     2005/10/19: the program was reformed so that the alignment
*                 results are not dependent on the specific compilers.
*     2006/06/20: select 'A' if there is altLoc when reading PDB file.
*     2007/02/27: rotation matrix from Chain-1 to Chain-2 was added.
*     2007/04/18: added options with TM-score normalized by average
*                 length, shorter length, or longer length of two 
*                 structures.
*     2007/05/23: added additional output file 'TM.sup_all' for showing
*                 full-chain C-alpha traces while 'TM.sup' is only for 
*                 aligned regions.
*     2007/09/19: added a new feature alignment to deal with the problem
*                 of aligning fractional structures (e.g. protein
*                 interfaces).
*     2007/10/16: A bug for irregular bond-length models was fixed.
*     2009/03/14: A new initial alignment was added and previous initial
*                 alignments are further enhanced. This change increased
*                 accuracy by 4% but increases CPU cost by 40%.
*     2009/08/20: A bug for asymmetry alignment result was fixed.
*     2010/08/02: A new RMSD matrix was used to remove obsolete statements.
*                 Staled subroutines were deleted.
*     2011/01/03: The length of pdb file names were extended to 500.
*     2011/01/24: Fixed a bug on output file name created on 2011/01/03.
*     2011/01/30: An open source license is attached to the program.
*     2011/09/03: A new option "-d" is added to allow users to change
*                 TM-score normalization scale. A version number is attached 
*                 to the program from now on.
*     2011/10/11: A new scale (d0) was introduced for alignment search. This
*                 is to mainly improve alignment selection for small proteins 
*                 (e.g. L<50 residues) but also increase alignment coverage
*                 of larger proteins. Second, TM-align output format is changed
*                 and two TM-scores normalized by both chains are reported.
*     2011/10/12: Distance cutoff for gap is increased from 3.85A to 4.25A.
*                 Added 'TMalign -v' to allow user to check version number.
*     2012/01/24: Fix a bug for secondary structure definition
*     2012/04/16: Add an option to allow user to specify seed alignments, e.g.
*                 '-i align.txt'. This is used together with other inherent
*                 TM-align seeds. An example of the fasta file can be seen at
*                 http://zhanglab.ccmb.med.umich.edu/TM-align/align.txt.
*     2012/04/17: Add an option '-m matrix.txt' to output the rotation matrix
*                 in separate file, drop-off secondary-structure smooth 
*                 procedure, and add one iteration in initial5. This change 
*                 increases the alignment accuracy (TM-score) by 2%.
*     2012/04/19: Add additional output file 'TM.sup_atm' 'TM.sup_all_atm' for 
*                 showing all-atom superposition while 'TM.sup' and 'TM.sup_all' 
*                 are only for C-alpha traces.
*     2012/05/07: Improved RMSD calculation subroutine which speeds up TM-algin
*                 program by 10%.
*     2012/07/07: Add an option '-I align.txt' to allow user to STICK TO the
*                 inital alignment. This is different from '-i align.txt' where
*                 initial alignment can be optimized.
*     2013/05/08: Update TM-align so that it can read all alternate location
*                 indicators and residue insertions.
*     2013/05/11: Fixed a bug in array overflow.
*     2014/06/01: Added 'TM.sup_all_atm_lig' to display ligand structures
*     2015/09/14: Optimized I/O which increased speed by ~100%
*     2016/05/21: Fixed a bug on conformation output
*     2017/07/08: Added one iteration in initial4 to avoid asymmetric alignment
*     2019/07/08: Enable TM-align to support both PDB and mmCIF formats, and
*                 fixed a bug on file output
*     2019/08/18: Fixed multiple bugs associated with mmCIF formats
*     2019/08/22: added output scripts for pymol, C++ version was included.
**************************************************************************
      
c        1         2         3         4         5         6         7 !
c 3456789012345678901234567890123456789012345678901234567890123456789012345678

      program TMalign
      PARAMETER(nmax=5000)      !maximum length of the sequence
      PARAMETER(nmax2=10000)    !for alignment output
      
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/length/nseq1,nseq2
      common/d0/d0,anseq
      common/d0min/d0_min
      common/d00/d00,d002
      
      common/alignment/m_alignment,sequence(10),TM_ali,L_ali,rmsd_ali
      common/alignment1/m_alignment_stick
      character*10000 sequence
      
      character seq1(0:nmax),seq2(0:nmax),du1,du3,du4*2
      character*500 fnam,pdb(100),outname,falign,fmatrix
      character*3 aa(-1:20),aanam,ss(2,nmax)
      character*500 s,du,dum1,dum2
      character*504 xxx,xxx_p,xxx_pdb
      character aseq1(nmax2),aseq2(nmax2),aseq3(nmax2)
      character*8 version
      character*5000 s1         !maximum length of protein is 5000

ccc   mmCIF
      character*500 ctmp(1000),ch(nmax),ent(nmax),mn(nmax)
      character*500 ch_t,ent_t,mn_t
      character*20 Aatom(2,90000)
      character Agroup(2,90000)*6
      character Ares(2,90000)*3,Aalt(2,90000)
      character Ains(2,90000),Aent(2,90000)
      character Cins(2,nmax),Cch(2)
      character*2 Ach(2,90000),chid(2,90000)
      
      integer Aatomi(2,90000),Aresi(2,90000),n_cut(2)
      dimension iform(10),xx(2,90000),yy(2,90000),zz(2,90000)
      dimension nL(2),nLCA(2)
ccc
      
      dimension m1(nmax),m2(nmax),m12(2,nmax)
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/init/invmap_i(nmax)

      common/TM/TM,TMmax
      common/d8/d8
      common/initial4/mm(2,nmax)
      
      character*10 aa1,ra1,aa2,ra2,du2
      dimension ia1(90000),aa1(90000),ra1(90000),ir1(90000)
      dimension xa1(90000),ya1(90000),za1(90000)
      dimension ia2(90000),aa2(90000),ra2(90000),ir2(90000)
      dimension xa2(90000),ya2(90000),za2(90000)
      
      dimension ma1(nmax),ma2(nmax)
      dimension nc1(nmax),nc2(nmax)
      
      dimension nres1(2,nmax,32:122),nres2(2,nmax,32:122,32:122) !number of atoms
      character*5 atom1(50)     !atom name

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms !armsd is real
      data w /nmax*1.0/
ccc   

      data aa/ 'BCK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP','CYX'/
      character*1 slc(-1:20)
      data slc/'X','G','A','S','C',
     &     'V','T','I','P','M',
     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W','C'/

      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
         write(*,*)
         write(*,*)'Brief instruction for running TM-align program:'
         write(*,*)'(For detail: Zhang & Skolnick, Nucl. Acid. Res.',
     &     ' 33: 2302-9, 2005)'
         write(*,*)
         write(*,*)'1. Align ''chain_1.pdb'' and ''chain_2.pdb'':'
         write(*,*)'   >TMalign chain_1.pdb chain_2.pdb'
         write(*,*)
         write(*,*)'2. Ask TM-align to start with an alignment',
     &        ' specified in fasta file ''align.txt'':'
         write(*,*)'   >TMalign chain_1.pdb chain_2.pdb -i align.txt'
         write(*,*)'   or to stick the alignment to ''align.txt'':'
         write(*,*)'   >TMalign chain_1.pdb chain_2.pdb -I align.txt'
         write(*,*)
         write(*,*)'3. Output the superposition to ''TM.sup'', ',
     &        '''TM.sup_all'' and ''TM.sup_atm'':'
         write(*,*)'   >TMalign chain_1.pdb chain_2.pdb -o TM.sup'
         write(*,*)'      To view superimposed C-alpha traces of',
     &        ' aligned regions by rasmol or pymol:'
         write(*,*)'        >rasmol -script TM.sup'
         write(*,*)'        >pymol -d @TM.sup.pml'
         write(*,*)'      To view superimposed C-alpha traces of',
     &        ' all regions:'
         write(*,*)'        >rasmol -script TM.sup_all'
         write(*,*)'        >pymol -d @TM.sup_all.pml'
         write(*,*)'      To view superimposed full-atom structures of',
     &        ' aligned regions:'
         write(*,*)'        >rasmol -script TM.sup_atm'
         write(*,*)'        >pymol -d @TM.sup_atm.pml'
         write(*,*)'      To view superimposed full-atom structures of',
     &        ' all regions:'
         write(*,*)'        >rasmol -script TM.sup_all_atm'
         write(*,*)'        >pymol -d @TM.sup_all_atm.pml'
         write(*,*)'      To view superimposed full-atom structures of',
     &        ' all regions with ligands:'
         write(*,*)'        >rasmol -script TM.sup_all_atm_lig'
         write(*,*)'        >pymol -d @TM.sup_all_atm_lig.pml'
         write(*,*)
         write(*,*)'4. There are two TM-scores reported. You ',
     &        'should use the one normalized by'
         write(*,*)'   the length of the protein you ',
     &        'are interested in.'
         write(*,*)'   If you want TM-score normalized by the ',
     &        'average length of two proteins:'
         write(*,*)'      >TMalign chain_1.pdb chain_2.pdb -a'
         write(*,*)'   or TM-score normalized by an ',
     &        'assigned length (>L_min), e.g. 100 AA:'
         write(*,*)'      >TMalign chain_1.pdb chain_2.pdb -L 100'
         write(*,*)'   If you want TM-score scaled by an assigned d0,',
     &        ' e.g. 5 A:'
         write(*,*)'      >TMalign chain_1.pdb chain_2.pdb -d 5'
         write(*,*)
         write(*,*)'5. Output TM-align rotation matrix:'
         write(*,*)'   >TMalign chain_1.pdb chain_2.pdb -m matrix.txt'
         write(*,*)
         goto 9999
      endif
      
      version='20190822'
      if(fnam.eq.'-v')then
         write(*,*)'TM-align Version ',version
         goto 9999
      endif
      
******* options ----------->
      m_out=-1                  !decided output
      m_fix=-1                  !fixed length-scale only for output
      m_ave=-1                  !using average length
      m_d0_min=-1               !diminum d0 for search
      m_d0=-1                   !given d0 for output
      m_alignment=-1            !without initial alignment
      m_alignment_stick=-1      !without initial alignment
      m_matrix=-1               !no output of matrix
      narg=iargc()
      i=0
      j=0
 115  continue
      i=i+1
      call getarg(i,fnam)
      if(fnam.eq.'-o')then
         m_out=1
         i=i+1
         call getarg(i,outname)
      elseif(fnam.eq.'-a')then  !change superposed output but not the alignment
         m_ave=1
         i=i+1
      elseif(fnam.eq.'-L')then  !change both L_all and d0
         m_fix=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)L_fix
      elseif(fnam.eq.'-d')then
         m_d0=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)d0_fix
      elseif(fnam.eq.'-i')then
         m_alignment=1
         i=i+1
         call getarg(i,fnam)
         falign=fnam
      elseif(fnam.eq.'-I')then
         m_alignment_stick=1
         i=i+1
         call getarg(i,fnam)
         falign=fnam
      elseif(fnam.eq.'-m')then
         m_matrix=1
         i=i+1
         call getarg(i,fnam)
         fmatrix=fnam
      else
         j=j+1
         pdb(j)=fnam
      endif
      if(i.lt.narg)goto 115

      irmx=0                    !no conformation output
      if(m_matrix.eq.1.or.m_out.eq.1)then !we need to extract rotation matrix
         irmx=1
      endif
      
**********decide file format (PDB or mmCIF) ------------------->
      do j=1,2
         iform(j)=0             !format of pdb(j)
         open(unit=10,file=pdb(j),status='old')
         do while(iform(j).eq.0)
            read(10,*)s
            if(s(1:6).eq.'HEADER'.or.s(1:4).eq.'ATOM'.or.
     &           s(1:6).eq.'REMARK')then
               iform(j)=1       !PDB format
            elseif(s(1:4).eq.'data'.or.s(1:1).eq.'#'.or.
     &           s(1:5).eq.'loop_')then
               iform(j)=2       !mmCIF format
            endif
         enddo
         if(iform(j).eq.0)then
            write(*,*)'error: file must in PDB or mmCIF format!'
         endif
         close(10)
      enddo
*******^^^^ format is decided ^^^^^^^^^^^^^^^^^^^^^^^^^^
      
cccccccccRead data from CA file ---------------->
c     we only need to read following (keep first chain, keep only one altLoc):
c     xa(i,j,k)----(x,y,z)
c     mm(2,nmax)---residue order number, for gapless threading
c     ss(2,nmax)---residue name ('GLY') for seq_ID calculation and output
c     seq1(0,nmax),seq2(0,nmax)----single characters ('G', 'N')
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 1001 ic=1,2            !ic=1,2 for file1 and file2
         i=0
         open(unit=10,file=pdb(ic),status='old')
         if(iform(ic).eq.1)then !file in PDB format------->
            do while (.true.)   !start do-while ---------->
               read(10,'(A500)',end=1013) s
               if(i.gt.0)then
                  if(s(1:3).eq.'TER'.or.s(1:3).eq.'MOD'.
     &                 or.s(1:3).eq.'END')then
                     goto 1013  !only read the first chain
                  endif
               endif
               if(s(1:4).eq.'ATOM')then !read 'ATOM' =======>
                  read(s(13:16),*)du4 !will remove space before 'CA'
                  if(du4.eq.'CA')then !read 'CA' ---------->
                     du1=s(27:27) !residue insertion tag
*******remove repeated altLoc atoms -------->
                     mk=1
                     if(s(17:17).ne.' ')then !with Alternate atom
                        read(s(23:26),*)i8 !res ID
                        if(nres1(ic,i8,ichar(du1)).ge.1)then !since CA, cannot be >1
                           mk=-1 !this residue was already read, skip it
                        endif
                     endif
c^^^^^^^^^altLoc checked (mk.ne.1 will be skipped) ^^^^^^^^^
                     if(mk.eq.1)then !mk=1 ------------>
                        i=i+1
                        read(s,'(a6,I5,a6,A3,A1,A1,i4,A1,a3,3F8.3)')
     &                       du,itmp,du,aanam,du,chid(ic,i),mm(ic,i),
     &                       Cins(ic,i),du,
     &                       xa(1,i,ic-1),xa(2,i,ic-1),xa(3,i,ic-1)
***   
                        i8=mm(ic,i)
                        nres1(ic,i8,ichar(du1))=
     &                       nres1(ic,i8,ichar(du1))+1 !nres1 only for check altLoc
***   
                        ss(ic,i)=aanam !residue name, 'GLY', for seq_ID
                        if(i.ge.nmax)goto 1013
                     endif      !<-----mk=1
                  endif         !<---- read 'CA'
               endif            !<==== read 'ATOM'
            enddo               !<----end do-while
         else                   !<----mmCIF format,if(iform(1).eq.2)
            in=0                !number of entries
            do while (.true.)   !start do-while ---------->
               read(10,'(A500)',end=1013) s
               if(i.gt.0)then   !skip unuseful read to save time
                  if(s(1:1).eq.'#'.or.s(1:5).eq.'loop_'.or.
     &                 s(1:6).eq.'HETATM')goto 1013
               endif
               if(s(1:11).eq.'_atom_site.')then
                  in=in+1
                  read(s,*)du
                  if(du.eq.'_atom_site.label_atom_id')i_atom=in !CA,O, no need
                  if(du.eq.'_atom_site.label_alt_id')i_alt=in !'.',A,B
                  if(du.eq.'_atom_site.label_comp_id')i_res=in !GLY,LEU
                  if(du.eq.'_atom_site.label_asym_id')i_ch=in !A,B
                  if(du.eq.'_atom_site.auth_asym_id')i_ch=in !A,B, using later one
                  if(du.eq.'_atom_site.label_entity_id')i_ent=in !1,2,a,b
                  if(du.eq.'_atom_site.label_seq_id')i_resi=in !1,2,3
                  if(du.eq.'_atom_site.auth_seq_id')i_resi=in !1,2,3, identical to PDB res
                  if(du.eq.'_atom_site.pdbx_PDB_ins_code')i_ins=in !A,?
                  if(du.eq.'_atom_site.Cartn_x')i_x=in !x, 1.234
                  if(du.eq.'_atom_site.Cartn_y')i_y=in !y, 1.234
                  if(du.eq.'_atom_site.Cartn_z')i_z=in !z, 1.234
                  if(du.eq.'_atom_site.pdbx_PDB_model_num')i_mn=in !model number, 1,2,3
               endif
               if(s(1:4).eq.'ATOM')then !read 'ATOM' =======>
                  read(s,*)(ctmp(j),j=1,in) !no space before characters
                  if(ctmp(i_atom).eq.'CA')then !read 'CA' ---------->
                     if(i.gt.0)then
                        if(i_mn.gt.1)then !sometimes it may have no i_mn
                           read(ctmp(i_mn),*)mn_t
                           if(mn_t.ne.mn(i)) goto 1013 !only read first model
                        endif
                        read(ctmp(i_ch),*)ch_t
                        if(ch_t.ne.ch(i)) goto 1013 !only read first chain
                        read(ctmp(i_ent),*)ent_t
                        if(ent_t.ne.ent(i)) goto 1013 !only read first entity
                     endif
                     
*******remove repeated altLoc atoms (this does not care inserted residues) -------->
                     mk=1
                     du1=ctmp(i_ins) !read insertion tag
                     if(ctmp(i_alt).ne.'.')then !with Alternate atom
                        read(ctmp(i_resi),*)i8 !res ID
                        if(nres1(ic,i8,ichar(du1)).ge.1)then !since CA, cannot be >1
                           mk=-1 !this residue was already read, skip it
                        endif
                     endif
c^^^^^^^^^altLoc checked (mk.ne.1 will be skipped) ^^^^^^^^^
                     if(mk.eq.1)then !mk=1 ------------>
                        i=i+1
                        read(ctmp(i_ch),*)ch(i) !for check other chain
                        read(ctmp(i_ent),*)ent(i) !for check other entity
                        read(ctmp(i_mn),*)mn(i) !for check model number
***   
                        read(ctmp(i_x),*)xa(1,i,ic-1)
                        read(ctmp(i_y),*)xa(2,i,ic-1)
                        read(ctmp(i_z),*)xa(3,i,ic-1)
                        read(ctmp(i_resi),*)mm(ic,i) !residue order, 4,5,6
                        read(ctmp(i_ins),*)Cins(ic,i) !residue insertion, A, ?
                        if(Cins(ic,i).eq.'?')Cins(ic,i)=' '
                        ss(ic,i)=ctmp(i_res) !residue name, 'GLY', for seq_ID
***   
                        i8=mm(ic,i)
                        chid(ic,i)=ch(i)
                        nres1(ic,i8,ichar(du1))=
     &                       nres1(ic,i8,ichar(du1))+1 !nres1 only for check altLoc
***   
                        if(i.ge.nmax)goto 1013
                     endif      !<-----mk=1
                  endif         !<-----read 'CA'
               endif            !<==== read 'ATOM'
            enddo               !<----end do-while
         endif                  !if(iform(ic).eq.2)
 1013    continue
         close(10)

c-------convert 'GLY' to 'G' ------------>
         if(ic.eq.1)then
            nseq1=i
            do i=1,nseq1
               do j=-1,20
                  if(ss(1,i).eq.aa(j))then
                     seq1(i)=slc(j)
                     goto 121
                  endif
               enddo
               seq1(i)=slc(-1)
 121           continue
            enddo
         else
            nseq2=i
            do i=1,nseq2
               do j=-1,20
                  if(ss(2,i).eq.aa(j))then
                     seq2(i)=slc(j)
                     goto 122
                  endif
               enddo
               seq2(i)=slc(-1)
 122           continue
            enddo
         endif
 1001 continue                  !ic=1,2 for file1 and file2

c^^^^^^^^^^^^^read input files completed ^^^^^^^^^^^^^^^^^^^^^^

ccccccccc read initial alignment file from 'alignment.txt':
      if(m_alignment.eq.1.or.m_alignment_stick.eq.1)then
         open(unit=10,file=falign,status='old')
         n_p=0
         do while (.true.)
            read(10,'(A5000)',end=1012)s1
            if(s1(1:1).eq.">")then
               n_p=n_p+1
               sequence(n_p)=''
               if(n_p.gt.2)goto 1012
            else
               if(n_p.gt.0)then
                  sequence(n_p)=sequence(n_p)(1:len_trim(sequence(n_p)))
     &                 //s1(1:len_trim(s1))
               endif
            endif
         enddo
 1012    continue
         close(10)
         if(n_p.lt.2)then
            write(*,*)'ERROR: FASTA format is wrong, two proteins',
     &           ' should be included'
            stop
         endif
         if(len_trim(sequence(1)).ne.len_trim(sequence(2)))then
            write(*,*)'Warning: FASTA format may be wrong, the',
     &           ' length in alignment should be equal. But we',
     &           ' run it anyway'
         endif
      endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
*     Scale of TM-score in search based on length of smaller protein --------->
      anseq_min=min(nseq1,nseq2) !both search and d8_cut use nseq_min
      anseq=anseq_min           !length for defining TMscore in search
      d8=1.5*anseq_min**0.3+3.5 !remove pairs with dis>d8 during search & final
      if(anseq.gt.19)then       !L=19, d0=0.168
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=0.168
      endif
      d0_min=d0+0.8             !best for search, this should be changed when calculate real TM-score
      if(d0.lt.d0_min)d0=d0_min !min d0 in search=0.968, min d0 in output=0.5
      d00=d0                    !for quickly calculate TM-score in searching
      if(d00.gt.8)d00=8
      if(d00.lt.4.5)d00=4.5
      d002=d00**2
      nseq=max(nseq1,nseq2)
      
***** do alignment **************************
      CALL TM_align             !to find invmap(j)
      
************************************************************
***   Refine alignment by cutting dis>d8 ------------------------>
      n_al=0
      do j=1,nseq2
         if(invmap0(j).gt.0)then
            i=invmap0(j)
            n_al=n_al+1
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
            m1(n_al)=i          !for recording residue order
            m2(n_al)=j
         endif
      enddo
      d0_input=d0               !scaled by seq_min
      call TMscore8(d0_input,n_al,xtm1,ytm1,ztm1,n_al,
     &     xtm2,ytm2,ztm2,TM,Rcomm,Lcomm) !TM-score with dis<d8 only
      j=0
      n_eq=0
      do i=1,n_al
         dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &        (ztm1(i)-ztm2(i))**2)
         if(dis2.le.d8.or.m_alignment_stick.eq.1)then
            j=j+1
            xtm1(j)=xtm1(i)
            ytm1(j)=ytm1(i)
            ztm1(j)=ztm1(i)
            xtm2(j)=xtm2(i)
            ytm2(j)=ytm2(i)
            ztm2(j)=ztm2(i)

            r_1(1,j)=xtm1(i)
            r_1(2,j)=ytm1(i)
            r_1(3,j)=ztm1(i)
            r_2(1,j)=xtm2(i)
            r_2(2,j)=ytm2(i)
            r_2(3,j)=ztm2(i)
            
            m1(j)=m1(i)         !record alignment
            m2(j)=m2(i)
            if(ss(1,m1(i)).eq.ss(2,m2(i)))then
               n_eq=n_eq+1
            endif
         endif
      enddo
      n8_al=j
      seq_id=float(n_eq)/(n8_al+0.00000001)
      call u3b(w,r_1,r_2,n8_al,0,rms,u,t,ier)
      rmsd=dsqrt(rms/n8_al)
      
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*^^^^^^^ alignment is done, all cutoffs were based on shorter chain^^^^^^^^
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

************************************************************
***   Output TM-score -------------------------->
      d0_out=5                  !only for showing residue-pair distance
      d0_min=0.5                !for TM-score output, consistent stdrd TM-score
*     Based on Chain_1===>
      anseq=nseq1
      if(anseq.gt.21)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      d0_input=d0
      call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n8_al,
     &     xtm2,ytm2,ztm2,TM8,Rcomm,Lcomm,0) !normal TMscore
      TM1=TM8*n8_al/anseq
*     Based on Chain_2===>
      anseq=nseq2
      if(anseq.gt.21)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      d0_input=d0
      call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n8_al,
     &     xtm2,ytm2,ztm2,TM8,Rcomm,Lcomm,irmx) !normal TMscore
      TM2=TM8*n8_al/anseq
*     Based on Average length===>
      if(m_ave.eq.1)then
         anseq=(nseq1+nseq2)/2.0
         if(anseq.gt.21)then
            d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
         else
            d0=d0_min
         endif
         if(d0.lt.d0_min)d0=d0_min
         d0_input=d0
         call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n8_al,
     &        xtm2,ytm2,ztm2,TM8,Rcomm,Lcomm,irmx) !normal TMscore
         TM12=TM8*n8_al/anseq
      endif
*     Based on assigned length===>
      if(m_fix.eq.1)then
         anseq=L_fix            !input length
         if(anseq.gt.21)then
            d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
         else
            d0=d0_min
         endif
         if(d0.lt.d0_min)d0=d0_min
         d0_input=d0
         call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n8_al,
     &        xtm2,ytm2,ztm2,TM8,Rcomm,Lcomm,irmx) !normal TMscore
         TML=TM8*n8_al/anseq
      endif
*     Based on user-specified d0===>
      if(m_d0.eq.1)then
         d0=d0_fix
         d0_out=d0_fix
         d0_input=d0
         call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n8_al,
     &        xtm2,ytm2,ztm2,TM8,Rcomm,Lcomm,irmx) !normal TMscore
         TMfix=TM8*n8_al/nseq2
      endif
      
*****************************************
*     screen output of TM-align results *
*****************************************
      write(*,*)
      write(*,*)'*****************************************************',
     &     '*********************'
      write(*,*)'*                        TM-align (Version ',version,
     &     ')                     *'
      write(*,*)'* An algorithm for protein structure alignment and co',
     &     'mparison            *'
      write(*,*)'* Based on statistics:                               ',
     &     '                    *'
      write(*,*)'*       0.0 < TM-score < 0.30, random structural simi',
     &     'larity              *'
      write(*,*)'*       0.5 < TM-score < 1.00, in about the same fold',
     &     '                    *'
      write(*,*)'* Reference: Y Zhang and J Skolnick, Nucl Acids Res 3',
     &     '3, 2302-9 (2005)    *'
      write(*,*)'* Please email your comments and suggestions to: zhng',
     &     '@umich.edu          *'
      write(*,*)'*****************************************************',
     &     '*********************'
      write(*,*)
      write(*,101)pdb(1)
 101  format('Name of Chain_1: ',A50)
      write(*,102)pdb(2)
 102  format('Name of Chain_2: ',A50)
      write(*,103)nseq1
 103  format('Length of Chain_1: ',I4,' residues')
      write(*,201)nseq2
 201  format('Length of Chain_2: ',I4,' residues')
      
      if(m_alignment.eq.1.or.m_alignment_stick.eq.1)then
         write(*,72)TM_ali,L_ali,rmsd_ali
 72      format('User-specified initial alignment: TM/Lali/rmsd= ',
     &        f7.5,', ',I4,', ',f6.3)
      endif
      
      write(*,*)
      write(*,203)n8_al,rmsd,seq_id
 203  format('Aligned length= ',I4,', RMSD= ',f6.2,
     &     ', Seq_ID=n_identical/n_aligned= ',f5.3)
      write(*,204)TM1
 204  format('TM-score= ',f7.5,' (if normalized by length of Chain_1)')
      write(*,205)TM2
 205  format('TM-score= ',f7.5,' (if normalized by length of Chain_2)')
      if(m_ave.eq.1)then
         write(*,206)TM12,(nseq1+nseq2)/2.0
 206     format('TM-score= ',f7.5,
     &        ' (if normalized by average length of chains =',f6.1,')')
      endif
      if(m_fix.eq.1)then
         write(*,207)TML,L_fix
 207     format('TM-score= ',f7.5,
     &        ' (if scaled by user-specified L=',I4,')')
      endif
      if(m_d0.eq.1)then
         write(*,208)TMfix,d0_fix
 208     format('TM-score= ',f7.5,
     &        ' (if scaled by user-specified d0=',f4.1,')')
      endif
      write(*,210)
 210  format('(You should use TM-score normalized by length',
     &     ' of the reference protein)')
      write(*,*)

********* extract rotation matrix based on TMscore8 ------------>
      if(m_matrix.eq.1.or.m_out.eq.1)then !we need to extract rotation matrix
         L=0
         do i=1,n8_al
            k=m1(i)
            L=L+1
            r_1(1,L)=xa(1,k,0)
            r_1(2,L)=xa(2,k,0)
            r_1(3,L)=xa(3,k,0)
            r_2(1,L)=xtm1(i)
            r_2(2,L)=ytm1(i)
            r_2(3,L)=ztm1(i)
         enddo
         if(L.le.3)then
            write(*,*)'Aligned length is too short,',
     &           ' no matrix outout'
            goto 211
         endif
         call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2, this will be used by whole-chain
 211     continue
      endif

********* output rotation matrix -------------------------->
      if(m_matrix.eq.1)then
         open(unit=1,file=fmatrix,status='unknown')
         write(1,*)'-------- Rotation matrix to rotate Chain_1 to ',
     &        'Chain_2 ------'
         write(1,*)'m          t(m)         u(m,1)         u(m,2) ',
     &        '        u(m,3)'
         do i=1,3
            write(1,209)i,t(i),u(i,1),u(i,2),u(i,3)
         enddo
         write(1,*)'Code for rotating Chain_1 from (x,y,z) to (X,Y,Z):'
         write(1,*)'   do i=1,L'
         write(1,*)'     X(i)=t(1)+u(1,1)*x(i)+u(1,2)*y(i)+u(1,3)*z(i)'
         write(1,*)'     Y(i)=t(2)+u(2,1)*x(i)+u(2,2)*y(i)+u(2,3)*z(i)'
         write(1,*)'     Z(i)=t(3)+u(3,1)*x(i)+u(3,2)*y(i)+u(3,3)*z(i)'
         write(1,*)'   enddo'
         write(1,*)
         close(1)
      endif
 209  format(I2,f18.10,f15.10,f15.10,f15.10)
      
************  output aligned sequences **************************
      ii=0
      i1_old=1
      i2_old=1
      do i=1,n8_al
         do j=i1_old,m1(i)-1
            ii=ii+1
            aseq1(ii)=seq1(j)
            aseq2(ii)='-'
            aseq3(ii)=' '
         enddo
         do j=i2_old,m2(i)-1
            ii=ii+1
            aseq1(ii)='-'
            aseq2(ii)=seq2(j)
            aseq3(ii)=' '
         enddo
         ii=ii+1
         aseq1(ii)=seq1(m1(i))
         aseq2(ii)=seq2(m2(i))
         dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &     (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
         if(dis2.le.d0_out)then
           aseq3(ii)=':'
         else
           aseq3(ii)='.'
         endif
         i1_old=m1(i)+1
         i2_old=m2(i)+1
      enddo
      do i=i1_old,nseq1
         ii=ii+1
         aseq1(ii)=seq1(i)
         aseq2(ii)='-'
         aseq3(ii)=' '
      enddo
      do i=i2_old,nseq2
         ii=ii+1
         aseq1(ii)='-'
         aseq2(ii)=seq2(i)
         aseq3(ii)=' '
      enddo
      write(*,50)d0_out
 50   format('(":" denotes aligned residue pairs of d < ',f3.1,
     &     ' A, "." denotes other aligned residues)')
      write(*,10)(aseq1(i),i=1,ii)
      write(*,10)(aseq3(i),i=1,ii)
      write(*,10)(aseq2(i),i=1,ii)
 10   format(10000A1)
      write(*,*)

c^^^^^^^^^^screen output is done ^^^^^^^^^^^^^^^^^^^^^^^^^^

      if(m_out.ne.1)then
         stop
      endif
     
*******************************************************
*     output alignment structure for Rasmal review    *
*******************************************************

c*****************************************************************
c************* Step-1: read all-atom structures ******************
c*****************************************************************
      do 1002 ic=1,2            !ic=1,2 for file1 and file2
         nL(ic)=0               !number of lines in PDB file
         nLCA(ic)=0             !number of lines with CA in PDB file
         n_cut(ic)=0            !end of first chain, decide to show aaa_all_atm_lig
         open(unit=10,file=pdb(ic),status='old')
         if(iform(ic).eq.1)then !file in PDB format %%%%%%%%%%------->
            do while (.true.)
               read(10,'(A500)',end=1015) s
               if(n_cut(ic).eq.0.and.nL(ic).gt.0)then !decide n_cut
                  if(s(1:3).eq.'TER'.or.s(1:3).eq.'MOD'.
     &                 or.s(1:3).eq.'END')then
                     n_cut(ic)=nL(ic)
                  endif
               endif
               
***   skip model_number >=2 ----------->
               if(s(1:5).eq.'MODEL')then
                  read(s,*)du,model_num
                  if(model_num.gt.1)then
                     do while(.true.)
                        read(10,'(A500)',end=1015) s
                        if(s(1:6).eq.'ENDMDL')goto 1014
                     enddo
                  endif
               endif
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
               
               if(s(1:6).eq.'ATOM  '.or.s(1:6).eq.'HETATM')then !read ATOM/HETATM ====>
*******remove repeated altLoc atoms -------->
                  mk=1
                  du1=s(27:27)  !read insertion tag
                  du3=s(22:22)  !chain ID
                  if(s(17:17).ne.' ')then !with Alternate atom
                     du2=s(13:16) !atom ID
c                     read(s(13:16),*)du2
                     read(s(23:26),*)i8 !res ID
                     do i1=1,nres2(ic,i8,ichar(du1),ichar(du3)) !#of atoms for res_insert
                        if(du2.eq.atom1(i1))then !such atom was already read 
                           mk=-1
                        endif
                     enddo
                  endif
c^^^^^^^^^altLoc checked (mk.ne.1) will be skipped) ^^^^^^^^^
                  
                  if(mk.eq.1)then !mk=1--------->
                     nL(ic)=nL(ic)+1
                     read(s,'(a6,I5,a1,A4,a1,A3,a1,A1,I4,A1,a3,3F8.3)')
     &                    Agroup(ic,nL(ic)),Aatomi(ic,nL(ic)),
     &                    du,Aatom(ic,nL(ic)),Aalt(ic,nL(ic)),
     &                    Ares(ic,nL(ic)),du,Ach(ic,nL(ic)),
     &                    Aresi(ic,nL(ic)),Ains(ic,nL(ic)),du,
     &                    xx(ic,nL(ic)),yy(ic,nL(ic)),zz(ic,nL(ic))
c                     read(Aatom(ic,nL(ic)),*)Aatom(ic,nL(ic)) !remove space before ' CA'
***   
                     i8=Aresi(ic,nL(ic))
                     if(nres2(ic,i8,ichar(du1),ichar(du3)).lt.50)then
                        nres2(ic,i8,ichar(du1),ichar(du3))=
     &                       nres2(ic,i8,ichar(du1),ichar(du3))+1 !#of atoms for res_ins
                        atom1(nres2(ic,i8,ichar(du1),ichar(du3)))=
     &                       Aatom(ic,nL(ic)) !atom ID
                     endif
***   
                     if(nL(ic).ge.90000)goto 1015
                  endif         !<------mk=1
               endif            !<======read ATOM/HETATM
 1014          continue         !skip model_num >1
            enddo               !do while
         else                   !<----mmCIF format,if(iform(1).eq.2) %%%%%%%%%%%%%%%%
            in=0                !number of entries
            do while (.true.)   !start do-while ---------->
               read(10,'(A500)',end=1015) s
               if(nL(ic).gt.0)then !skip unuseful read to save time
                  if(s(1:1).eq.'#'.or.s(1:5).eq.'loop_')goto 1015
               endif
               if(s(1:11).eq.'_atom_site.')then
                  in=in+1
                  read(s,*)du
                  if(du.eq.'_atom_site.group_PDB')i_group=in !'ATOM','HETATM'
                  if(du.eq.'_atom_site.id')i_atomi=in !1,2,3,4
                  if(du.eq.'_atom_site.type_symbol')i_type=in !N,C,O
                  if(du.eq.'_atom_site.label_atom_id')i_atom=in !CA,O
                  if(du.eq.'_atom_site.label_alt_id')i_alt=in !'.',A,B
                  
                  if(du.eq.'_atom_site.label_comp_id')i_res=in !GLY,LEU
                  if(du.eq.'_atom_site.label_asym_id')i_ch=in !A,B
                  if(du.eq.'_atom_site.auth_asym_id')i_ch=in !A,B, using later one
                  if(du.eq.'_atom_site.label_entity_id')i_ent=in !1,2,a,b
                  if(du.eq.'_atom_site.label_seq_id')i_resi=in !1,2,3
                  if(du.eq.'_atom_site.auth_seq_id')i_resi=in !1,2,3, identical to PDB res
                  if(du.eq.'_atom_site.pdbx_PDB_ins_code')i_ins=in !A,?
                  
                  if(du.eq.'_atom_site.Cartn_x')i_x=in !x, 1.234
                  if(du.eq.'_atom_site.Cartn_y')i_y=in !y, 1.234
                  if(du.eq.'_atom_site.Cartn_z')i_z=in !z, 1.234
                  if(du.eq.'_atom_site.pdbx_PDB_model_num')i_mn=in !model number, 1,2,3
               endif
               
               if(s(1:4).eq.'ATOM'.or.s(1:6).eq.'HETATM')then !read 'ATOM' =======>
                  read(s,*)(ctmp(j),j=1,in)
***   skip model_num >1 --------------------------------------->
                  if(i_mn>1)then
                     read(ctmp(i_mn),*)model_num
                     if(model_num.gt.1)goto 1016 !skip model_num >1
                  endif
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                  
                  if(n_cut(ic).eq.0.and.nL(ic).gt.0.and.
     &                 nLCA(ic).gt.0)then !decide n_cut for mmcif
                     read(ctmp(i_ch),*)ch_t
                     read(ctmp(i_ent),*)ent_t
                     if(ch_t.ne.Ach(ic,nL(ic)).or.ent_t.ne.
     &                    Aent(ic,nL(ic)))then
                        n_cut(ic)=nL(ic)
                     endif
                  endif
                  
*******remove repeated altLoc atoms (this does not care inserted residues) -------->
                  mk=1
                  du1=ctmp(i_ins) !read insertion tag
                  du3=ctmp(i_ch) !chain ID
                  if(ctmp(i_alt).ne.'.')then !with Alternate atom
                     du2=ctmp(i_atom) !atom ID
                     read(ctmp(i_resi),*)i8 !res ID
                     do i1=1,nres2(ic,i8,ichar(du1),ichar(du3)) !#of atoms for res_insert
                        if(du2.eq.atom1(i1))then !such atom was already read 
                           mk=-1
                        endif
                     enddo
                  endif
c^^^^^^^^^altLoc checked (mk.ne.1) will be skipped) ^^^^^^^^^

                  if(mk.eq.1)then !mk=1 -------------->
                     nL(ic)=nL(ic)+1
                     read(ctmp(i_group),*)Agroup(ic,nL(ic))
                     read(ctmp(i_atomi),*)Aatomi(ic,nL(ic))
                     read(ctmp(i_atom),*)Aatom(ic,nL(ic))
********add space before Aatom(ic,nL(ic)) for output ------>
                     if(Aatom(ic,nL(ic))(4:4).eq.' ')then
                        Aatom(ic,nL(ic))=' '//Aatom(ic,nL(ic))
                     endif
c     read(ctmp(i_alt),*)Aalt(ic,nL(ic)) ! not used, because we alway use 1 atom for altLoc
c     if(Aalt(ic,nL(ic)).eq.'.')Aalt(ic,nL(ic))=' '
                     if(Aatom(ic,nL(ic)).eq.' CA')then
                        nLCA(ic)=nLCA(ic)+1 !for deciding n_cut of mmCIF
                     endif
                     
                     read(ctmp(i_res),*)Ares(ic,nL(ic))
                     read(ctmp(i_ch),*)Ach(ic,nL(ic)) !for check other chain
                     read(ctmp(i_ent),*)Aent(ic,nL(ic)) !for check other entity
                     if(ctmp(i_resi)(1:1).eq.'.')ctmp(i_resi)='0' !for HETATM
                     read(ctmp(i_resi),*)Aresi(ic,nL(ic))
                     read(ctmp(i_ins),*)Ains(ic,nL(ic))
                     if(Ains(ic,nL(ic)).eq.'?')Ains(ic,nL(ic))=' '
***   
                     read(ctmp(i_x),*)xx(ic,nL(ic))
                     read(ctmp(i_y),*)yy(ic,nL(ic))
                     read(ctmp(i_z),*)zz(ic,nL(ic))
***   
                     i8=Aresi(ic,nL(ic))
                     if(nres2(ic,i8,ichar(du1),ichar(du3)).lt.50)then
                        nres2(ic,i8,ichar(du1),ichar(du3))=
     &                       nres2(ic,i8,ichar(du1),ichar(du3))+1 !#of atoms for res_ins
                        atom1(nres2(ic,i8,ichar(du1),ichar(du3)))=
     &                       Aatom(ic,nL(ic)) !atom ID
                     endif
***   
                     if(nL(ic).ge.90000)goto 1015
                  endif         !<---- mk=1
 1016             continue      !skip model_num >1
               endif            !<==== read 'ATOM'
            enddo               !<----end do-while
         endif                  !<----mmCIF format,if(iform(1).eq.2) %%%%%%%%%
 1015    continue
         close(10)
 1002 continue                  !ic=1,2 for file1 and file2
      
********** mark aligned residues --------------->
      do i=1,n8_al
         m12(1,i)=m1(i)
         m12(2,i)=m2(i)
      enddo
      Cch(1)='A'
      Cch(2)='B'
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c*********************************************************************
c************* Step-2: output structures for Rasmol ******************
c*********************************************************************
      do 1019 ip=1,5
         if(ip.eq.1) xxx=outname(1:len_trim(outname)) !'aaa', aligned CA
         if(ip.eq.2) xxx=outname(1:len_trim(outname))//'_all' !'aaa_all', all CA
         if(ip.eq.3) xxx=outname(1:len_trim(outname))//'_atm' !'aaa_atm', all aligned atoms
         if(ip.eq.4) xxx=outname(1:len_trim(outname))//'_all_atm' !'aaa_all_atm', all atoms
         if(ip.eq.5) xxx=outname(1:len_trim(outname))//'_all_atm_lig' !'aaa_all_atm_lig', atom+ligand
         OPEN(unit=10,file=xxx,status='unknown')
         
*************for pymol preset ------>
         xxx_p=xxx(1:len_trim(xxx))//'.pml' !for pymol script
         xxx_pdb=xxx(1:len_trim(xxx))//'.pdb' !for pymol pdb
         OPEN(unit=11,file=xxx_p,status='unknown')
         OPEN(unit=12,file=xxx_pdb,status='unknown')
         
         write(11,'(A,A)')'load ',xxx_pdb(1:len_trim(xxx_pdb))
         write(11,'(A)')'hide all'
         write(11,'(A)')'bg_color white'
         write(11,'(A)')'color blue, chain A'
         write(11,'(A)')'color red, chain B'
         write(11,'(A)')'set transparency=0.2'
         write(11,'(A)')'set sphere_scale, 0.25'
         write(11,'(A)')'set stick_radius, 0.3'
*^^^^^^^^^ pymol preset complete ^^^^^^^^^^^^^^^^^^^
         
***   script:
         if(ip.eq.1.or.ip.eq.2)then !'aaa','aaa_all'
            write(10,'(A)')'load inline'
            write(10,'(A)')'select *A'
            write(10,'(A)')'wireframe .45'
            write(10,'(A)')'select *B'
            write(10,'(A)')'wireframe .20'
            write(10,'(A)')'select all'
            write(10,'(A)')'color white'
            do i=1,n8_al
               dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &              (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
               if(dis2.le.d0_out)then
                  write(10,902)mm(1,m1(i)),mm(2,m2(i)) !select residue
                  write(10,'(A)')'color red'
               endif
            enddo
            write(10,'(A)')'select all'
            
            write(11,'(A)')'show sticks, chain A'
            write(11,'(A)')'show sticks, chain B'
         elseif(ip.eq.3.or.ip.eq.4)then !'aaa_atm', 'aaa_all_atm'
            write(10,'(A)')'load inline'
            write(10,'(A)')'select *A'
            write(10,'(A)')'color blue'
            write(10,'(A)')'select *B'
            write(10,'(A)')'color red'
            write(10,'(A)')'select all'
            write(10,'(A)')'cartoon'
            
            write(11,'(A)')'show cartoon, chain A'
            write(11,'(A)')'show cartoon, chain B'
         else                   !'aaa_all_atm_lig
            write(10,'(A)')'load inline'
            write(10,'(A)')'select all'
            write(10,'(A)')'cartoon'
            write(10,'(A)')'select *A'
            write(10,'(A)')'color blue'
            write(10,'(A)')'select *B'
            write(10,'(A)')'color red'
            write(10,'(A)')'select ligand'
            write(10,'(A)')'wireframe 0.25'
            write(10,'(A)')'select solvent'
            write(10,'(A)')'spacefill 0.25'
            write(10,'(A)')'select all'
            
            write(11,'(A)')'show cartoon, chain A'
            write(11,'(A)')'show cartoon, chain B'
            write(11,'(A)')'show stick, HETATM'
            write(11,'(A)')'show spheres, solvent'
         endif
         write(10,'(A)')'exit'
         write(10,903)version
         write(10,104)pdb(1),nseq1
         write(10,105)pdb(2),nseq2,int(anseq),d0
         write(10,106)n8_al,rmsd,TM2,seq_id

c------------output coordinates ----------->         
         do ic=1,2
            if(ic.eq.1)then     !na=1-5000
               na=0             !for CONECT
            else                !na=5001-10000
               na=5000          !for CONECT
            endif
            
            if(n_cut(ic).eq.0)then
               n_cut(ic)=nL(ic) !there is no cut
            endif
            
            do i=1,nL(ic)
c---------- decide if we should output this line ------------>
               mk=0
               if(ip.eq.1)then  !'aaa', aligned CA
                  if(i.le.n_cut(ic))then
                     if(Aatom(ic,i).eq.' CA')then
                        do j=1,n8_al
                           if(Aresi(ic,i).eq.mm(ic,m12(ic,j)))then !residue order
                              if(Ains(ic,i).eq.Cins(ic,m12(ic,j)))then !residue insertion
                                 mk=1
                                 goto 120 !each line output once
                              endif
                           endif
                        enddo
                     endif
                  endif
               elseif(ip.eq.2)then !'aaa_all', all CA
                  if(i.le.n_cut(ic))then
                     if(Aatom(ic,i).eq.' CA')then
                        mk=1
                     endif
                  endif
               elseif(ip.eq.3)then !'aaa_atm', all aligned protein atoms
                  if(i.le.n_cut(ic))then
                     do j=1,n8_al
                        if(Ach(ic,i).eq.chid(ic,m12(ic,j)))then ! residue order (e.g. 100)
                           if(Aresi(ic,i).eq.mm(ic,m12(ic,j)))then ! residue order (e.g. 100)
                              if(Ains(ic,i).eq.Cins(ic,m12(ic,j)))then !resi insertion (eg, A)
                                 mk=1
                                 goto 120 !each line output once
                              endif
                           endif
                        endif
                     enddo
                  endif
               elseif(ip.eq.4)then !'aaa_all_atm', all protein atoms
                  if(Ach(ic,i).eq.chid(ic,1))then ! residue order (e.g. 100)
                     if(i.le.n_cut(ic))then
                        mk=1
                     endif
                  endif
               elseif(ip.eq.5)then !'aaa_all_atm_lig', first chain and all ligands
                  if(i.le.n_cut(ic).or.Agroup(ic,i).eq.'HETATM')then
                     mk=1
                  endif
               endif
 120           continue

c^^^^^^^^^^
               if(mk.eq.1)then  !mk=1----------->
                  if(ic.eq.1)then
                     ax=t(1)+u(1,1)*xx(ic,i)+u(1,2)*yy(ic,i)+
     &                    u(1,3)*zz(ic,i)
                     ay=t(2)+u(2,1)*xx(ic,i)+u(2,2)*yy(ic,i)+
     &                    u(2,3)*zz(ic,i)
                     az=t(3)+u(3,1)*xx(ic,i)+u(3,2)*yy(ic,i)+
     &                    u(3,3)*zz(ic,i)
                  else
                     ax=xx(ic,i)
                     ay=yy(ic,i)
                     az=zz(ic,i)
                  endif
                  
c--------------output (x,y,z) --------------------------------->                  
                  na=na+1       !number of atoms for CONECT
                  if(ip.eq.1.or.ip.eq.2)then !need CONECT
                     write(10,1236)Agroup(ic,i),na,'',
     &                    Aatom(ic,i),'',Ares(ic,i),'',Cch(ic),
     &                    Aresi(ic,i),Ains(ic,i),'',ax,ay,az
                     write(12,1236)Agroup(ic,i),na,'',
     &                    Aatom(ic,i),'',Ares(ic,i),'',Cch(ic),
     &                    Aresi(ic,i),Ains(ic,i),'',ax,ay,az
                  else          !all atoms
                     write(10,1236)Agroup(ic,i),Aatomi(ic,i),'',
     &                    Aatom(ic,i),'',Ares(ic,i),'',Cch(ic),
     &                    Aresi(ic,i),Ains(ic,i),'',ax,ay,az
                     write(12,1236)Agroup(ic,i),Aatomi(ic,i),'',
     &                    Aatom(ic,i),'',Ares(ic,i),'',Cch(ic),
     &                    Aresi(ic,i),Ains(ic,i),'',ax,ay,az
                  endif
               endif            !<----mk=1
            enddo               !do i=1,nL(ic)

            write(10,'(A)')'TER' !TER
            write(12,'(A)')'TER' !TER
            if(ip.eq.1.or.ip.eq.2)then
               if(ic.eq.1)then
                  do i=2,na
                     write(10,'(A6,I5,I5)')'CONECT',i-1,i !CONECT atom numbers
                     write(12,'(A6,I5,I5)')'CONECT',i-1,i !CONECT atom numbers
                  enddo
               else
                  do i=5002,na
                     write(10,'(A6,I5,I5)')'CONECT',i-1,i !CONECT atom numbers
                     write(12,'(A6,I5,I5)')'CONECT',i-1,i !CONECT atom numbers
                  enddo
               endif
            endif
         enddo                  !do ic=1,2
         close(10)
         close(11)
         close(12)
 1019 continue                  !ip=1,5

 902     format('select ',I4,':A,',I4,':B')
 903     format('REMARK TM-align Version ',A8,'')
 104     format('REMARK Chain 1:',A10,'  Size=',I4)
 105     format('REMARK Chain 2:',A10,'  Size=',I4,
     &        ' (TM-score is normalized by ',I4,', d0=',f6.2,')')
 106     format('REMARK Aligned length=',I4,', RMSD=',f6.2,
     &        ', TM-score=',f7.5,', ID=',f5.3)
 1236    format(A6,I5,a1,A4,a1,A3,a1,A1,I4,A1,a3,3F8.3)

*^^^^^^^^^^^^^^^^^^ output finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 9999 END


***********************************************************************
***********************************************************************
*     Structure superposition
***********************************************************************
***********************************************************************
***********************************************************************
      SUBROUTINE TM_align
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      common/d0/d0,anseq
      dimension gapp(100)
      
      common/alignment/m_alignment,sequence(10),TM_ali,L_ali,rmsd_ali
      common/alignment1/m_alignment_stick
      character*10000 sequence
      
      TMmax=0
      n_gapp=2
      gapp(1)=-0.6
      gapp(2)=0
      
      ddcc=0.4
      if(anseq.le.40)then
         ddcc=0.1
      endif
      
ccc   stick to the initial alignment -------------->
      if(m_alignment_stick.eq.1)then
         do j=1,nseq2
            invmap(j)=-1
         enddo
         i1=0
         i2=0
         L1=LEN_TRIM(sequence(1))
         L2=LEN_TRIM(sequence(2))
         L=L1
         if(L2.lt.L)L=L2
         do 6661 i=1,L
            if(sequence(1)(i:i).ne.'-')i1=i1+1
            if(i1.gt.nseq1)then
               goto 6661
            endif
            if(sequence(2)(i:i).ne.'-')then
               i2=i2+1
               if(i2.gt.nseq2)then
                  goto 6661
               endif
               if(sequence(1)(i:i).ne.'-')then
                  invmap(i2)=i1
               endif
            endif
 6661    enddo
         call standard_TMscore(TM_ali,L_ali,rmsd_ali) !calc TM-score from invmap, nmlzd by nseq2
ccc   
         call get_score         !TM, matrix score(i,j)
         if(TM.gt.TMmax)then
            TMmax=TM
            do j=1,nseq2
               invmap0(j)=invmap(j)
            enddo
         endif
         return
      endif
      
*11111111111111111111111111111111111111111111111111111111
*     get initial alignment from global gapless threading
**********************************************************
      call get_initial1          !gapless threading
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, matrix score(i,j)
      
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif
      
*****************************************************************
*       initerative alignment, for different gap_open:
*****************************************************************
      DO 1 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 11 id=1,30          !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
            
            call get_score      !calculate TM-score, score(i,j)
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 111
            endif
            TM_old=TM
 11      continue
 111     continue
 1    continue
      
*222222222222222222222222222222222222222222222222222222222
*     get initial alignment from secondary structure alignment
**********************************************************
      call get_initial2         !DP for secondary structure
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif
      if(TM.le.TMmax*0.2)goto 2222
*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 2 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 22 id=1,30          !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
            
            call get_score      !calculate TM-score, score(i,j)
c     write(*,21)gap_open,rmsd_al,n_al,TM
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 222
            endif
            TM_old=TM
 22      continue
 222     continue
 2    continue
 2222 continue
      
*555555555555555555555555555555555555555555555555555555555555555555
*     get initial alignment of local structure superposition
*******************************************************************
      call get_initial5
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, matrix score(i,j)
      
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif
      if(TM.le.TMmax*ddcc)goto 5555
*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 5 i_gapp=1,n_gapp      !different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 55 id=1,2           !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
            
            call get_score      !calculate TM-score, score(i,j)
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 555
            endif
            TM_old=TM
 55      continue
 555     continue
 5    continue
 5555 continue
      
*333333333333333333333333333333333333333333333333333333333333
*     get initial alignment from invmap0+SS
*************************************************************
      call get_initial3         !invmap0+SS
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif
      if(TM.le.TMmax*ddcc)goto 3333
*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 3 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 33 id=1,30          !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)

            call get_score      !calculate TM-score, score(i,j)
c     write(*,21)gap_open,rmsd_al,n_al,TM
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 333
            endif
            TM_old=TM
 33      continue
 333     continue
 3    continue
 3333 continue

*444444444444444444444444444444444444444444444444444444444
*     initial alignment from gapless threading on largest continous fragments
**********************************************************
      call get_initial4         !gapless threading
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, matrix score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif
      if(TM.le.TMmax*ddcc)goto 4444
*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 4 i_gapp=2,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 44 id=1,2           !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
            
            call get_score      !calculate TM-score, score(i,j)
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
 44      continue
 4    continue
 4444 continue
      
*666666666666666666666666666666666666666666666666666666666666
*     get initial alignment from user's input:
*************************************************************
      if(m_alignment.ne.1)goto 6666
      do j=1,nseq2
         invmap(j)=-1
      enddo
      i1=0
      i2=0
      L1=LEN_TRIM(sequence(1))
      L2=LEN_TRIM(sequence(2))
c     write(*,*)'seq1= ',trim(sequence(1))
c     write(*,*)'seq2= ',trim(sequence(2))
      L=L1
      if(L2.lt.L)L=L2
      do 666 i=1,L
c         write(*,*)i,sequence(1)(i:i),sequence(2)(i:i)
         if(sequence(1)(i:i).ne.'-')i1=i1+1
         if(sequence(2)(i:i).ne.'-')then
            i2=i2+1
            if(i2.gt.nseq2)then
               goto 666
            endif
            if(sequence(1)(i:i).ne.'-')then
               invmap(i2)=i1
c               write(*,*)i2,i1
            endif
         endif
 666  enddo
ccc   
c      L_ali=0
c      do j=1,nseq2
c     write(*,*)j,invmap(j)
c         if(invmap(j).gt.0)then
c            L_ali=L_ali+1
c         endif
c      enddo
c     write(*,*)'L_ali=',L_ali
      call standard_TMscore(TM_ali,L_ali,rmsd_ali) !calc TM-score from invmap, nmlzd by nseq2
ccc   
      call get_score            !TM, matrix score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif
*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 6 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 66 id=1,30          !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
            
            call get_score      !calculate TM-score, score(i,j)
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
 66      continue
 6    continue
 6666 continue
      
c^^^^^^^^^^^^^^^ best alignment invmap0(j) found ^^^^^^^^^^^^^^^^^^
      RETURN
      END

**************************************************************
*     get initial alignment invmap0(i) from gapless threading
**************************************************************
      subroutine get_initial1
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      
      aL=min(nseq1,nseq2)
      idel=aL/2.0               !minimum size of considered fragment
      if(idel.le.5)idel=5
      n1=-nseq2+idel
      n2=nseq1-idel
      GL_max=0
      do ishift=n1,n2
         L=0
         do j=1,nseq2
            i=j+ishift
            if(i.ge.1.and.i.le.nseq1)then
               L=L+1
               invmap(j)=i
            else
               invmap(j)=-1
            endif
         enddo
         if(L.ge.idel)then
            call get_GL(GL)
            if(GL.gt.GL_max)then
               GL_max=GL
               do i=1,nseq2
                  invmap_i(i)=invmap(i)
               enddo
            endif
         endif
      enddo
      
      return
      end

**************************************************************
*     get initial alignment invmap0(i) from secondary structure
**************************************************************
      subroutine get_initial2
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/sec/isec(nmax),jsec(nmax)
      common/init/invmap_i(nmax)

********** assign secondary structures ***************
c     1->coil, 2->helix, 3->turn, 4->strand
      do i=1,nseq1
         isec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq1)then
            dis13=diszy(0,j1,j3)
            dis14=diszy(0,j1,j4)
            dis15=diszy(0,j1,j5)
            dis24=diszy(0,j2,j4)
            dis25=diszy(0,j2,j5)
            dis35=diszy(0,j3,j5)
            isec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo
      do i=1,nseq2
         jsec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq2)then
            dis13=diszy(1,j1,j3)
            dis14=diszy(1,j1,j4)
            dis15=diszy(1,j1,j5)
            dis24=diszy(1,j2,j4)
            dis25=diszy(1,j2,j5)
            dis35=diszy(1,j3,j5)
            jsec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo
      
********** score matrix **************************
      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               score(i,j)=1
            else
               score(i,j)=0
            endif
         enddo
      enddo
      
********** find initial alignment: invmap(j) ************
      gap_open=-1.0             !should be -1
      call DP(NSEQ1,NSEQ2)      !produce alignment invmap(j)
      do i=1,nseq2
         invmap_i(i)=invmap(i)
      enddo

*^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

**************************************************************
*     get initial alignment invmap0(i) from secondary structure 
*     and previous alignments
**************************************************************
      subroutine get_initial3
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/sec/isec(nmax),jsec(nmax)
      common/init/invmap_i(nmax)

********** score matrix **************************
      do i=1,nseq2
         invmap(i)=invmap0(i)
      enddo
      call get_score1           !get score(i,j) using RMSD martix
      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               score(i,j)=0.5+score(i,j)
            else
               score(i,j)=score(i,j)
            endif
         enddo
      enddo

********** find initial alignment: invmap(j) ************
      gap_open=-1.0             !should be -1
      call DP(NSEQ1,NSEQ2)      !produce alignment invmap(j)
      do i=1,nseq2
         invmap_i(i)=invmap(i)
      enddo

*^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

**************************************************************
*     get initial alignment invmap0(i) from fragment gapless threading
**************************************************************
      subroutine get_initial4
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      common/initial4/mm(2,nmax)
      logical contin

      dimension ifr2(2,nmax,nmax),Lfr2(2,nmax),i_fr2(2)
      dimension ifr(nmax)
      
      fra_min=4                 !>=4,minimum fragment for search
      fra_min1=fra_min-1        !cutoff for shift, save time
      dcu0=4.25
      
      GL_max=0
      
c      do k=1,2                  !k=1, fragment from protein1; k=2, fragment from protein2
      do k=2,1,-1               !k=1, fragment from protein1; k=2, fragment from protein2
ccc   Find the smallest continuous fragments on protein-k -------->
         dcu=dcu0               !breaking bond-length
         if(k.eq.1)then
            nseq0=nseq1
            r_min=nseq1/3.0     !minimum fragment, in case too small protein
         else
            nseq0=nseq2
            r_min=nseq2/3.0     !minimum fragment, in case too small protein
         endif
         if(r_min.gt.fra_min)r_min=fra_min
         
 20      nfr=1                  !number of fragments
         j=1                    !number of residues at nfr-fragment
         ifr2(k,nfr,j)=1        !residue ID of nfr-fragment
         Lfr2(k,nfr)=j          !length of the fragment
         
         do i=2,nseq0
            dis=diszy(k-1,i-1,i) !str,res,res
            contin=.false.
            if(dcu.gt.dcu0)then
               if(dis.lt.dcu)then
                  contin=.true.
               endif
            elseif(mm(k,i).eq.(mm(k,i-1)+1))then
               if(dis.lt.dcu)then
                  contin=.true.
               endif
            endif
            if(contin)then
               j=j+1
               ifr2(k,nfr,j)=i
               Lfr2(k,nfr)=j
            else
               nfr=nfr+1
               j=1
               ifr2(k,nfr,j)=i
               Lfr2(k,nfr)=j
            endif
         enddo
         Lfr_max=0
         i_fr2(k)=1             !ID of the maximum piece
         do i=1,nfr
            if(Lfr_max.lt.Lfr2(k,i))then
               Lfr_max=Lfr2(k,i)
               i_fr2(k)=i
            endif
         enddo
         if(Lfr_max.lt.r_min)then
            dcu=dcu+0.01
            goto 20
         endif
         
         L_fr=Lfr2(k,i_fr2(k))  !length of the maximum fragment
         do i=1,L_fr
            ifr(i)=ifr2(k,i_fr2(k),i)
         enddo
         
ccc   find the best initial alignment------------>
         if(k.eq.1)then         !using fragment from protein-1
            if(L_fr.eq.nseq1)then !to make it different from initial1
               n1=int(nseq1*0.1) !0
               n2=int(nseq1*0.89) !2
               j=0
               do i=n1,n2
                  j=j+1
                  ifr(j)=ifr(n1+j)
               enddo
               L_fr=j
            endif
            
            nseq1_=L_fr
            aL=min(nseq1_,nseq2)
            idel=aL/2.5         !minimum size of considered fragment
            if(idel.le.fra_min1)idel=fra_min1
            n1=-nseq2+idel      !shift1
            n2=nseq1_-idel      !shift2
c            write(*,*)idel,aL,n1,n2,'aaa---'
            do ishift=n1,n2
               L=0
               do j=1,nseq2
                  i=j+ishift
                  if(i.ge.1.and.i.le.nseq1_)then
                     L=L+1
                     invmap(j)=ifr(i)
                  else
                     invmap(j)=-1
                  endif
               enddo
               if(L.ge.idel)then
                  call get_GL(GL)
                  if(GL.gt.GL_max)then
                     GL_max=GL
                     do i=1,nseq2
                        invmap_i(i)=invmap(i)
                     enddo
                  endif
               endif
            enddo
            
c            write(*,*)'GL_max=',GL_max
c            do i=1,nseq2
c               write(*,*)i,invmap_i(i),'111'
c            enddo
            
         else
            if(L_fr.eq.nseq2)then !to make it different from initial1
               n1=int(nseq2*0.1) !0
               n2=int(nseq2*0.89) !2
               j=0
               do i=n1,n2
                  j=j+1
                  ifr(j)=ifr(n1+j)
               enddo
               L_fr=j
            endif
            
            nseq2_=L_fr
            aL=min(nseq1,nseq2_)
            idel=aL/2.5         !minimum size of considered fragment
            if(idel.le.fra_min1)idel=fra_min1
            n1=-nseq2_+idel
            n2=nseq1-idel
c            write(*,*)idel,aL,n1,n2,'bbb----'
            do ishift=n1,n2
               L=0
               do j=1,nseq2
                  invmap(j)=-1
               enddo
               do j=1,nseq2_
                  i=j+ishift
                  if(i.ge.1.and.i.le.nseq1)then
                     L=L+1
                     invmap(ifr(j))=i
                  endif
               enddo
               if(L.ge.idel)then
                  call get_GL(GL)
                  if(GL.gt.GL_max)then
                     GL_max=GL
                     do i=1,nseq2
                        invmap_i(i)=invmap(i)
                     enddo
                  endif
               endif
            enddo

c            write(*,*)'GL_max=',GL_max
c            do i=1,nseq2
c               write(*,*)i,invmap_i(i),'1112222'
c            enddo
         endif
      enddo
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      return
      end

**************************************************************
*    fifth initial alignement. Using local structure super   *
*    position.                                               *
**************************************************************
      subroutine get_initial5
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/d0/d0,anseq
      common/d0min/d0_min
      common/init/invmap_i(nmax)
      common/sec/isec(nmax),jsec(nmax)
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms
      data w /nmax*1.0/
      common/inv/invmap_a(nmax)
      integer aL,m1,m2,n_frag
      dimension n_frag(10)
      
***** setting parameters ************************************
      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01
      GLmaxA=0
      aL=min(nseq1,nseq2)
      
c     jump on sequence1 -------------->
      if(nseq1.gt.250)then
         n_jump1=45
      elseif(nseq1.gt.200)then
         n_jump1=35
      elseif(nseq1.gt.150)then
         n_jump1=25
      else
         n_jump1=15
      endif
      if(n_jump1.gt.nseq1/3)n_jump1=nseq1/3
      
c     jump on sequence2 -------------->
      if(nseq2.gt.250)then
         n_jump2=45
      elseif(nseq2.gt.200)then
         n_jump2=35
      elseif(nseq2.gt.150)then
         n_jump2=25
      else
         n_jump2=15
      endif
      if(n_jump2.gt.nseq2/3)n_jump2=nseq2/3
      
c     fragment to superimpose -------------->
      n_frag(1)=20
      n_frag(2)=100
      if(n_frag(1).gt.aL/3)n_frag(1)=aL/3
      if(n_frag(2).gt.aL/2)n_frag(2)=aL/2

c     start superimpose search -------------->
      do i_frag=1,2
         m1=nseq1-n_frag(i_frag)+1
         m2=nseq2-n_frag(i_frag)+1
         do ii=1,m1,n_jump1
            do jj=1,m2,n_jump2
               do k=1,n_frag(i_frag)
                  iii=ii+k-1
                  jjj=jj+k-1
                  r_1(1,k)=xa(1,iii,0)
                  r_1(2,k)=xa(2,iii,0)
                  r_1(3,k)=xa(3,iii,0)
                  r_2(1,k)=xa(1,jjj,1)
                  r_2(2,k)=xa(2,jjj,1)
                  r_2(3,k)=xa(3,jjj,1)
               enddo
               
*********superpose the two structures and rotate it *****************
               call u3b(w,r_1,r_2,n_frag(i_frag),1,rms,u,t,ier) !u rotate r_1 to r_2
               do i1=1,nseq1
                  xx=t(1)+u(1,1)*xa(1,i1,0)+u(1,2)*xa(2,i1,0)+u(1,3)
     &                 *xa(3,i1,0)  
                  yy=t(2)+u(2,1)*xa(1,i1,0)+u(2,2)*xa(2,i1,0)+u(2,3)
     &                 *xa(3,i1,0)
                  zz=t(3)+u(3,1)*xa(1,i1,0)+u(3,2)*xa(2,i1,0)+u(3,3)
     &                 *xa(3,i1,0)
                  do j1=1,nseq2
                     dd=(xx-xa(1,j1,1))**2+(yy-xa(2,j1,1))**2+
     &                    (zz-xa(3,j1,1))**2
                     score(i1,j1)=1/(1+dd/d02) ! changing
                  enddo
               enddo
               
*********extract alignement with score(i,j) *****************
               call DP(NSEQ1,NSEQ2)
               call get_GL(GL)
               if(GL.gt.GLmaxA)then
                  GLmaxA=GL
                  do j1=1,nseq2
                     invmap_i(j1)=invmap(j1)
                  enddo
               endif

            enddo
         enddo
      enddo
      
      return
      end
      
*************************************************************
*     assign secondary structure:
*************************************************************
      function diszy(i,i1,i2)
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      diszy=sqrt((xa(1,i1,i)-xa(1,i2,i))**2
     &     +(xa(2,i1,i)-xa(2,i2,i))**2
     &     +(xa(3,i1,i)-xa(3,i2,i))**2)
      return
      end

*************************************************************
*     assign secondary structure:
*************************************************************
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
      make_sec=1
      delta=2.1
      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
                        make_sec=2 !helix
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
                        make_sec=4 !strand
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end

****************************************************************
*     quickly calculate TM-score with given invmap(i) in 3 iterations
****************************************************************
      subroutine get_GL(GL)
      PARAMETER(nmax=5000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      common/TM/TM,TMmax
      common/d00/d00,d002

      dimension xo1(nmax),yo1(nmax),zo1(nmax)
      dimension xo2(nmax),yo2(nmax),zo2(nmax)
      dimension dis2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
            xo1(n_al)=xa(1,i,0)
            yo1(n_al)=xa(2,i,0)
            zo1(n_al)=xa(3,i,0)
            xo2(n_al)=xa(1,j,1)
            yo2(n_al)=xa(2,j,1)
            zo2(n_al)=xa(3,j,1)
         endif
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      GL=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         GL=GL+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration------------->
      d002t=d002
 21   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
      if(j.lt.3.and.n_al.gt.3)then
         d002t=d002t+.5
         goto 21
      endif
      L=j
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
      G2=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G2=G2+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration------------->
      d002t=d002+1
 22   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
      if(j.lt.3.and.n_al.gt.3)then
         d002t=d002t+.5
         goto 22
      endif
      L=j
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
      G3=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G3=G3+1/(1+dis2(i)/(d0**2))
      enddo
      if(G2.gt.GL)GL=G2
      if(G3.gt.GL)GL=G3

c^^^^^^^^^^^^^^^^ GL done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

****************************************************************
*     with invmap(i) calculate TM-score and martix score(i,j) for rotation 
****************************************************************
      subroutine get_score
      PARAMETER(nmax=5000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/ut/u,t

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for TM-score:
            xtm1(n_al)=xa(1,i,0) !for TM-score
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
         endif
      enddo
***   calculate TM-score for the given alignment----------->
      d0_input=d0
      call TMscore8_search(d0_input,n_al,xtm1,ytm1,ztm1,
     &     n_al,xtm2,ytm2,ztm2,TM,Rcomm,Lcomm) !simplified search engine
      TM=TM*n_al/anseq          !TM-score
***   calculate score matrix score(i,j)------------------>
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d0**2)
         enddo
      enddo
      
c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

****************************************************************
*     with invmap(i) calculate score(i,j) using RMSD rotation 
****************************************************************
      subroutine get_score1
      PARAMETER(nmax=5000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      common/d0min/d0_min
      common/TM/TM,TMmax

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo
***   calculate score matrix score(i,j)------------------>
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d02)
         enddo
      enddo

c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate TM-score for a given alignment specified by invmap (Based on Chain_2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine standard_TMscore(TM2,L_ali,RMSD)
      PARAMETER(nmax=5000)      !maximum length of the sequence
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/d0min/d0_min
      
ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms !armsd is real
      data w /nmax*1.0/
ccc   

      d0_min=0.5                !for TM-score output, consistent stdrd TM-score
      anseq=nseq2
      if(anseq.gt.21)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      d0_input=d0               !scaled by seq_min
      
cccc  collect aligned residues from invmap ------->
      n_al=0
      do j=1,nseq2
         if(invmap(j).gt.0)then
            i=invmap(j)
            n_al=n_al+1
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)

            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo
      call u3b(w,r_1,r_2,n_al,0,rms,u,t,ier)
      L_ali=n_al
      RMSD=dsqrt(rms/n_al)
c      write(*,*)'---------',rms,n_al,RMSD,u(1,1),t(1)
      
      call TMscore(d0_input,n_al,xtm1,ytm1,ztm1,n_al,
     &     xtm2,ytm2,ztm2,TM8,Rcomm,Lcomm,0) !normal TMscore
      TM2=TM8*n_al/anseq
      
c^^^^^^^^^^ TM-score calculation is done ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end
      
*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     TM--TM-score of the comparison
*     Rcomm--RMSD of two structures in the common aligned residues
*     Lcomm--Length of the common aligned regions
*
*     Note: 
*     1, Always put native as the second structure, by which TM-score
*        is normalized.
*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
*        TM-score superposition.
*************************************************************************
*************************************************************************
*** dis<8, simplified search engine
      subroutine TMscore8_search(dx,L1,x1,y1,z1,L2,x2,y2,z2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=5000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)
      dimension iL0(nmax)
      common/ut/u,t

      dimension x1(nmax),y1(nmax),z1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ***************************
*     because L1=L2 in this special case---------->
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

************/////
*     parameters:
*****************
***   d0------------->
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      
      n_it=20                   !maximum number of iterations
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        k=0
        do i=1,iL_max,40        !this is the simplification!
          k=k+1
          iL0(k)=i
        enddo
        if(iL0(k).lt.iL_max)then
          k=k+1
          iL0(k)=iL_max
        endif
        n_shift=k
        do 300 i_shift=1,n_shift
           iL=iL0(i_shift)
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           Rcomm=0              !not used
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun8       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun8    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M
      
******** return the final rotation ****************
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max
      
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     TM--TM-score of the comparison
*     Rcomm--RMSD of two structures in the common aligned residues
*     Lcomm--Length of the common aligned regions
*
*     Note: 
*     1, Always put native as the second structure, by which TM-score
*        is normalized.
*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
*        TM-score superposition.
*************************************************************************
*************************************************************************
***   dis<8, but same search engine
      subroutine TMscore8(dx,L1,x1,y1,z1,L2,x2,y2,z2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=5000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ***************************
*     because L1=L2 in this special case---------->
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

************/////
*     parameters:
*****************
***   d0------------->
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max    !on aligned residues, [1,nseqA]
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun8       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun8    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

******** return the final rotation ****************
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun8
      PARAMETER(nmax=5000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score
      common/d8/d8

      d_tmp=d
      d_tmp2=d*d
      d82=d8*d8
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis2=(xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2
         if(dis2.lt.d_tmp2)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         if(dis2.le.d82)then
            score_sum=score_sum+1/(1+dis2/d0/d0)
         endif
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
         d_tmp2=d_tmp*d_tmp
         goto 21
      endif
      score=score_sum/float(nseqB) !TM-score

      return
      end

*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     TM--TM-score of the comparison
*     Rcomm--RMSD of two structures in the common aligned residues
*     Lcomm--Length of the common aligned regions
*
*     Note: 
*     1, Always put native as the second structure, by which TM-score
*        is normalized.
*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
*        TM-score superposition.
*************************************************************************
*************************************************************************
***  normal TM-score:
      subroutine TMscore(dx,L1,x1,y1,z1,L2,x2,y2,z2,
     &     TM,Rcomm,Lcomm,irmx)
      PARAMETER(nmax=5000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ***************************
*     because L1=L2 in this special case---------->
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

************/////
*     parameters:
*****************
***   d0------------->
c      d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max      !on aligned residues, [1,nseqA]
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M
      
******** return the final rotation ****************
      if(irmx.eq.1)then          !we need coordinates for output structure
         LL=0
         do i=1,ka0
            m=k_ali0(i)         !record of the best alignment
            r_1(1,i)=xa(iA(m))
            r_1(2,i)=ya(iA(m))
            r_1(3,i)=za(iA(m))
            r_2(1,i)=xb(iB(m))
            r_2(2,i)=yb(iB(m))
            r_2(3,i)=zb(iB(m))
            LL=LL+1
         enddo
         call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
         do j=1,nseqA
            x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
            y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
            z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
         enddo
      endif
      
      TM=score_max

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun
      PARAMETER(nmax=5000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score
      
      d_tmp=d
      d_tmp2=d_tmp*d_tmp
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis2=(xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2
         if(dis2.lt.d_tmp2)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         score_sum=score_sum+1/(1+dis2/d0/d0)
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
         d_tmp2=d_tmp*d_tmp
         goto 21
      endif
      score=score_sum/float(nseqB) !TM-score

      return
      end

********************************************************************
*     Dynamic programming for alignment.
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
*     
*     Please note this subroutine is not a correct implementation of 
*     the N-W dynamic programming because the score tracks back only 
*     one layer of the matrix. This code was exploited in TM-align 
*     because it is about 1.5 times faster than a complete N-W code
*     and does not influence much the final structure alignment result.
*     In 1/1000 case, it may result in asymmetry, i.e. A_to_B!=B_to_A
*     For example, '1se9A.pdb' and '2edpA.pdb' (with score(i,j)=0 or 1 in SS)
********************************************************************
      SUBROUTINE DP(NSEQ1,NSEQ2)
      PARAMETER(nmax=5000)
      LOGICAL*1 DIR
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      dimension DIR(0:nmax,0:nmax),VAL(0:nmax,0:nmax)
      REAL H,V
***   initialize the matrix:
      val(0,0)=0
      do i=1,nseq1
        dir(i,0)=.false.
        val(i,0)=0
      enddo
      do j=1,nseq2
        dir(0,j)=.false.
        val(0,j)=0
        invmap(j)=-1
      enddo

***   decide matrix and path:
      DO j=1,NSEQ2
        DO i=1,NSEQ1
          D=VAL(i-1,j-1)+SCORE(i,j)
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN
          
          IF((D.GE.H).AND.(D.GE.V)) THEN
            DIR(I,J)=.true.
            VAL(i,j)=D
          ELSE
            DIR(I,J)=.false.
            if(V.GE.H)then
              val(i,j)=v
            else
              val(i,j)=h
            end if
          ENDIF
        ENDDO
      ENDDO
      
***   extract the alignment:
      i=NSEQ1
      j=NSEQ2
      DO WHILE((i.GT.0).AND.(j.GT.0))
        IF(DIR(i,j))THEN
          invmap(j)=i
          i=i-1
          j=j-1
        ELSE
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN
          IF(V.GE.H) THEN
            j=j-1
          ELSE
            i=i-1
          ENDIF
        ENDIF
      ENDDO
     
c^^^^^^^^^^^^^^^Dynamical programming done ^^^^^^^^^^^^^^^^^^^
      return
      END
      
cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m                    (given)
c  x    - x(i,m) are coordinates of atom c m in set x          (given)
c  y    - y(i,m) are coordinates of atom c m in set y          (given)
c  n    - n is number of atom pairs                            (given)
c  mode  - 0:calculate rms     only                            (given,short)
c          1:calculate     u,t only                            (given,medium)
c          2:calculate rms,u,t                                 (given,longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs            (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      double precision w(*), x(3,*), y(3,*)
      integer n, mode
      
      double precision rms, u(3,3), t(3)
      integer ier
      
      integer i, j, k, l, m1, m
      integer ip(9), ip2312(4)
      double precision r(3,3), xc(3), yc(3), wc
      double precision a(3,3), b(3,3), e(3), rr(6), ss(6)
      double precision e0, d, spur, det, cof, h, g
      double precision cth, sth, sqrth, p, sigma
      double precision c1x, c1y, c1z, c2x, c2y, c2z
      double precision s1x, s1y, s1z, s2x, s2y, s2z
      double precision sxx, sxy, sxz, syx, syy, syz, szx, szy, szz
      
      double precision sqrt3, tol, zero
      
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
      
      wc  = zero
      rms = zero
      e0  = zero
      s1x = zero
      s1y = zero
      s1z = zero
      s2x = zero
      s2y = zero
      s2z = zero
      sxx = zero
      sxy = zero
      sxz = zero
      syx = zero
      syy = zero
      syz = zero
      szx = zero 
      szy = zero
      szz = zero
      
      do i=1, 3
         xc(i) = zero
         yc(i) = zero
         t(i) = zero
         do j=1, 3
            r(i,j) = zero
            u(i,j) = zero
            a(i,j) = zero
            if( i .eq. j ) then
               u(i,j) = 1.0
               a(i,j) = 1.0
            end if
         end do
      end do
      
      ier = -1
      if( n .lt. 1 ) return
      ier = -2
      
      do m=1, n
         c1x=x(1, m)
         c1y=x(2, m)
         c1z=x(3, m)
         
         c2x=y(1, m)
         c2y=y(2, m)
         c2z=y(3, m)
         
         s1x = s1x + c1x
         s1y = s1y + c1y;
         s1z = s1z + c1z;
         
         s2x = s2x + c2x;
         s2y = s2y + c2y;
         s2z = s2z + c2z;
         
         sxx = sxx + c1x*c2x; 
         sxy = sxy + c1x*c2y; 
         sxz = sxz + c1x*c2z; 
         
         syx = syx + c1y*c2x; 
         syy = syy + c1y*c2y; 
         syz = syz + c1y*c2z;
         
         szx = szx + c1z*c2x; 
         szy = szy + c1z*c2y; 
         szz = szz + c1z*c2z;
      end do
      
      xc(1) = s1x/n;
      xc(2) = s1y/n; 	
      xc(3) = s1z/n;
      
      yc(1) = s2x/n;
      yc(2) = s2y/n; 	
      yc(3) = s2z/n;
      if(mode.eq.2.or.mode.eq.0) then ! need rmsd
         do m=1, n		
            do i=1, 3
               e0 = e0+ (x(i, m)-xc(i))**2 + (y(i, m)-yc(i))**2
            end do				
         end do
      endif
      
      r(1, 1) = sxx-s1x*s2x/n;
      r(2, 1) = sxy-s1x*s2y/n;
      r(3, 1) = sxz-s1x*s2z/n;
      r(1, 2) = syx-s1y*s2x/n;
      r(2, 2) = syy-s1y*s2y/n;
      r(3, 2) = syz-s1y*s2z/n;
      r(1, 3) = szx-s1z*s2x/n;
      r(2, 3) = szy-s1z*s2y/n;
      r(3, 3) = szz-s1z*s2z/n;
      
      det = r(1,1) * ( (r(2,2)*r(3,3)) - (r(2,3)*r(3,2)) )
     &     - r(1,2) * ( (r(2,1)*r(3,3)) - (r(2,3)*r(3,1)) )
     &     + r(1,3) * ( (r(2,1)*r(3,2)) - (r(2,2)*r(3,1)) )
      
      sigma = det
      
      m = 0
      do j=1, 3
         do i=1, j
            m = m+1
            rr(m) = r(1,i)*r(1,j) + r(2,i)*r(2,j) + r(3,i)*r(3,j)
         end do
      end do
      
      spur = (rr(1)+rr(3)+rr(6)) / 3.0
      cof = (((((rr(3)*rr(6) - rr(5)*rr(5)) + rr(1)*rr(6))
     &     - rr(4)*rr(4)) + rr(1)*rr(3)) - rr(2)*rr(2)) / 3.0
      det = det*det
      
      do i=1, 3
         e(i) = spur
      end do
      if( spur .le. zero ) goto 40
      d = spur*spur
      h = d - cof
      g = (spur*cof - det)/2.0 - spur*h
      if( h .le. zero ) then
         if( mode .eq. 0 ) then
            goto 50
         else
            goto 30
         end if
      end if
      sqrth = dsqrt(h)
      d = h*h*h - g*g
      if( d .lt. zero ) d = zero
      d = datan2( dsqrt(d), -g ) / 3.0
      cth = sqrth * dcos(d)
      sth = sqrth*sqrt3*dsin(d)
      e(1) = (spur + cth) + cth
      e(2) = (spur - cth) + sth
      e(3) = (spur - cth) - sth
      
      if( mode .eq. 0 ) then
         goto 50
      end if
      
      do l=1, 3, 2
         d = e(l)
         ss(1) = (d-rr(3)) * (d-rr(6))  - rr(5)*rr(5)
         ss(2) = (d-rr(6)) * rr(2)      + rr(4)*rr(5)
         ss(3) = (d-rr(1)) * (d-rr(6))  - rr(4)*rr(4)
         ss(4) = (d-rr(3)) * rr(4)      + rr(2)*rr(5)
         ss(5) = (d-rr(1)) * rr(5)      + rr(2)*rr(4)
         ss(6) = (d-rr(1)) * (d-rr(3))  - rr(2)*rr(2)
         
         if( dabs(ss(1)) .ge. dabs(ss(3)) ) then
            j=1
            if( dabs(ss(1)) .lt. dabs(ss(6)) ) j = 3
         else if( dabs(ss(3)) .ge. dabs(ss(6)) ) then
            j = 2
         else
            j = 3
         end if
         
         d = zero
         j = 3 * (j - 1)
         
         do i=1, 3
            k = ip(i+j)
            a(i,l) = ss(k)
            d = d + ss(k)*ss(k)
         end do
         if( d .gt. zero ) d = 1.0 / dsqrt(d)
         do i=1, 3
            a(i,l) = a(i,l) * d
         end do
      end do
      
      d = a(1,1)*a(1,3) + a(2,1)*a(2,3) + a(3,1)*a(3,3)
      if ((e(1) - e(2)) .gt. (e(2) - e(3))) then
         m1 = 3
         m = 1
      else
         m1 = 1
         m = 3
      endif
      
      p = zero
      do i=1, 3
         a(i,m1) = a(i,m1) - d*a(i,m)
         p = p + a(i,m1)**2
      end do
      if( p .le. tol ) then
         p = 1.0
         do 21 i=1, 3
            if (p .lt. dabs(a(i,m))) goto 21
            p = dabs( a(i,m) )
            j = i
 21      continue
         k = ip2312(j)
         l = ip2312(j+1)
         p = dsqrt( a(k,m)**2 + a(l,m)**2 )
         if( p .le. tol ) goto 40
         a(j,m1) = zero
         a(k,m1) = -a(l,m)/p
         a(l,m1) =  a(k,m)/p
      else
         p = 1.0 / dsqrt(p)
         do i=1, 3
            a(i,m1) = a(i,m1)*p
         end do
      end if
      
      a(1,2) = a(2,3)*a(3,1) - a(2,1)*a(3,3)
      a(2,2) = a(3,3)*a(1,1) - a(3,1)*a(1,3)
      a(3,2) = a(1,3)*a(2,1) - a(1,1)*a(2,3)
      
 30   do l=1, 2
         d = zero
         do i=1, 3
            b(i,l) = r(i,1)*a(1,l) + r(i,2)*a(2,l) + r(i,3)*a(3,l)
            d = d + b(i,l)**2
         end do
         if( d .gt. zero ) d = 1.0 / dsqrt(d)
         do i=1, 3
            b(i,l) = b(i,l)*d
         end do
      end do
      d = b(1,1)*b(1,2) + b(2,1)*b(2,2) + b(3,1)*b(3,2)
      p = zero
      
      do i=1, 3
         b(i,2) = b(i,2) - d*b(i,1)
         p = p + b(i,2)**2
      end do
      if( p .le. tol ) then
         p = 1.0
         do 22 i=1, 3
            if(p.lt.dabs(b(i,1)))goto 22
            p = dabs( b(i,1) )
            j = i
 22      continue
         k = ip2312(j)
         l = ip2312(j+1)
         p = dsqrt( b(k,1)**2 + b(l,1)**2 )
         if( p .le. tol ) goto 40
         b(j,2) = zero
         b(k,2) = -b(l,1)/p
         b(l,2) =  b(k,1)/p
      else
         p = 1.0 / dsqrt(p)
         do i=1, 3
            b(i,2) = b(i,2)*p
         end do
      end if
      
      b(1,3) = b(2,1)*b(3,2) - b(2,2)*b(3,1)
      b(2,3) = b(3,1)*b(1,2) - b(3,2)*b(1,1)
      b(3,3) = b(1,1)*b(2,2) - b(1,2)*b(2,1)
      
      do i=1, 3
         do j=1, 3
            u(i,j) = b(i,1)*a(j,1) + b(i,2)*a(j,2) + b(i,3)*a(j,3)
         end do
      end do
      
 40   do i=1, 3
         t(i) = ((yc(i) - u(i,1)*xc(1)) - u(i,2)*xc(2)) - u(i,3)*xc(3)
      end do
 50   do i=1, 3
         if( e(i) .lt. zero ) e(i) = zero
         e(i) = dsqrt( e(i) )
      end do
      
      ier = 0
      if( e(2) .le. (e(1) * 1.0d-05) ) ier = -1
      
      d = e(3)
      if( sigma .lt. 0.0 ) then
         d = - d
         if( (e(2) - e(3)) .le. (e(1) * 1.0d-05) ) ier = -1
      end if
      d = (d + e(2)) + e(1)
      
      if(mode .eq. 2.or.mode.eq.0) then ! need rmsd
         rms = (e0 - d) - d
         if( rms .lt. 0.0 ) rms = 0.0
      endif
      
      return
      end
      
