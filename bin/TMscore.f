*************************************************************************
*     This program is to compare two protein structures and identify the 
*     best superposition that has the highest TM-score. Input structures 
*     must be in the PDB format. By default, TM-score is normalized by 
*     the second protein. Users can obtain a brief instruction by simply
*     running the program without arguments. For comments/suggestions,
*     please contact email: zhng@umich.edu.
*     
*     Reference: 
*     Yang Zhang, Jeffrey Skolnick, Proteins, 2004 57:702-10.
*     
*     Permission to use, copy, modify, and distribute this program for 
*     any purpose, with or without fee, is hereby granted, provided that
*     the notices on the head, the reference information, and this
*     copyright notice appear in all copies or substantial portions of 
*     the Software. It is provided "as is" without express or implied 
*     warranty.
******************* Updating history ************************************
*     2005/10/19: the program was reformed so that the score values.
*                 are not dependent on the specific compilers.
*     2006/06/20: selected 'A' if there is altLoc when reading PDB file.
*     2007/02/05: fixed a bug with length<15 in TMscore_32.
*     2007/02/27: rotation matrix from Chain-1 to Chain-2 was added.
*     2007/12/06: GDT-HA score was added, fixed a bug for reading PDB.
*     2010/08/02: A new RMSD matrix was used and obsolete statement removed.
*     2011/01/03: The length of pdb file names were extended to 500.
*     2011/01/30: An open source license is attached to the program.
*     2012/05/07: Improved RMSD calculation subroutine which speeds up 
*                 TM-score program by 30%.
*     2012/06/05: Added option '-l L' which calculates TM-score (and maxsub
*                 and GDT scores) normalized by a specific length 'L'.
*     2012/12/17: Added 'TM.sup_atm' to superpose full-atom structures.
*                 The former superposition is for CA-trace only.
*     2013/05/08: Update TM-score so that it can read all alternate location
*                 indicators and residue insertions.
*     2013/05/11: Fix a bug in array overflow.
*     2016/03/23: Extended the program to allow calculating TM-score for for 
*                 complex structure comparisons, where multiple-chains are
*                 merged into a single chain. Chain ID is now included in
*                 the output files.
*     2019/07/08: Enabled TM-score to support both PDB and mmCIF formats,
*                 and updated structure reading which makes program faster.
*     2019/08/18: Fixed multiple bugs associated with mmCIF formats.
*     2019/08/22: added output scripts for pymol, C++ version was included.
*************************************************************************
      
c        1         2         3         4         5         6         7 !
c 3456789012345678901234567890123456789012345678901234567890123456789012345678
      
      program TMscore
      PARAMETER(nmax=5000)      !maximum number of residues
      PARAMETER(namax=50000)    !maximum of atoms
      
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)

      character*500 fnam,pdb(100),outname,xxx_p,xxx_pdb
      character*3 aa(-1:20),seqA(nmax),seqB(nmax),aanam,ss(2,nmax)
      character*500 s,du
      character*1 chA(nmax),chB(nmax),ch(nmax)
      character*1 chain1(namax),chain2(namax)
      character seq1A(nmax),seq1B(nmax),ali(nmax),du1,du3,du4*2
      character sequenceA(nmax),sequenceB(nmax),sequenceM(nmax)
      character ins1(nmax),ins2(nmax),ains1(namax),ains2(namax)
      
      dimension L_ini(100),iq(nmax)
      common/scores/score,score_maxsub,score_fix,score10
      common/GDT/n_GDT05,n_GDT1,n_GDT2,n_GDT4,n_GDT8
      double precision score,score_max,score_fix,score_fix_max
      double precision score_maxsub,score10
      dimension xa(nmax),ya(nmax),za(nmax)
      dimension x2(2,nmax),y2(2,nmax),z2(2,nmax)
      
      character*10 aa1,ra1,aa2,ra2,du2
      integer mm(2,nmax)

      dimension nres1(2,nmax,32:122),nres2(2,nmax,32:122,32:122) !number of atoms
      character*5 atom1(50)     !atom name
      dimension m12(2,nmax)
      
ccc   mmCIF
      character*500 ctmp(1000),ent(nmax),mn(nmax)
      character*500 ch_t,ent_t,mn_t
      character*20 Aatom(2,namax)
      character Agroup(2,namax)*6
      character Ares(2,namax)*3,Aalt(2,namax)
      character Ains(2,namax),Ach(2,namax),Aent(2,namax)
      character Cins(2,nmax)

      integer Aatomi(2,namax),Aresi(2,namax)
      dimension iform(10),xx(2,namax),yy(2,namax),zz(2,namax),nL(2)
ccc

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   
      
      double precision u2(3,3),t2(3) !for output

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

*****instructions ----------------->
      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
         write(*,*)
         write(*,*)'Brief instruction for running TM-score program:'
         write(*,*)'(For detail: Zhang & Skolnick,  Proteins, 2004',
     &        ' 57:702-10)'
         write(*,*)
         write(*,*)'1. Run TM-score to compare ''model'' and ',
     &        '''native'':'
         write(*,*)'   >TMscore model native'
         write(*,*)
         write(*,*)'2. Run TM-score to compare two complex structures',
     &        ' with multiple chains'
         write(*,*)'   (Compare all chains with the same chain',
     &        ' identifier)'
         write(*,*)'   >TMscore -c model native'
         write(*,*)
         write(*,*)'3. TM-score normalized with an assigned scale d0',
     &        ' e.g. 5 A:'
         write(*,*)'   >TMscore model native -d 5'
         write(*,*)
         write(*,*)'4. TM-score normalized by a specific length, ',
     &        'e.g. 120 AA:'
         write(*,*)'   >TMscore model native -l 120'
         write(*,*)
         write(*,*)'5. TM-score with superposition output, e.g. ',
     &        '''TM.sup'':'
         write(*,*)'   >TMscore model native -o TM.sup'
         write(*,*)'   To view superimposed CA-traces by rasmol ',
     &        'or pymol:'
         write(*,*)'      >rasmol -script TM.sup'
         write(*,*)'      >pymol -d @TM.sup.pml'
         write(*,*)'   To view superimposed atomic models:'
         write(*,*)'      >rasmol -script TM.sup_atm'
         write(*,*)'      >pymol -d @TM.sup_atm.pml'
         write(*,*)
         goto 9999
      endif
      
******* options ----------->
      m_out=-1
      m_fix=-1
      m_len=-1
      narg=iargc()
      i=0
      j=0
      m_complex=0
 115  continue
      i=i+1
      call getarg(i,fnam)
      if(fnam.eq.'-o')then
         m_out=1
         i=i+1
         call getarg(i,outname)
      elseif(fnam.eq.'-d')then
         m_fix=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)d0_fix
      elseif(fnam.eq.'-c')then
         m_complex=1
      elseif(fnam.eq.'-l')then
         m_len=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)l0_fix
      else
         j=j+1
         pdb(j)=fnam
      endif
      if(i.lt.narg)goto 115

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
     &              s(1:5).eq.'loop_')then
               iform(j)=2       !mmCIF format
            endif
         enddo
         if(iform(j).eq.0)then
            write(*,*)'error: file must in PDB or mmCIF format!'
         endif
         close(10)
      enddo
*******^^^^ format is decided ^^^^^^^^^^^^^^^^^^^^^^^^^^
      
c     write(*,*)'m_complex=',m_complex

cccccccccRead data from CA file ---------------->
c     we only need to read following (keep first chain, keep only one altLoc):
c     x,y,z2(ic,i)----(x,y,z)
c     mm(2,nmax)---residue order number, for gapless threading
c     ss(2,nmax)---residue name ('GLY') for seq_ID calculation and output
c     Cins(ic,i)---insertion
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 1001 ic=1,2            !ic=1,2 for file1 and file2
         i=0
         open(unit=10,file=pdb(ic),status='old')
         if(iform(ic).eq.1)then !file in PDB format------->
            do while (.true.)   !start do-while ---------->
               read(10,'(A500)',end=1013) s
               if(i.gt.0)then
                  if(m_complex.eq.0)then
                     if(s(1:3).eq.'TER'.or.s(1:3).eq.'MOD'.
     &                    or.s(1:3).eq.'END')then
                        goto 1013 !only read the first chain
                     endif
                  else
***   skip model_number >=2 ----------->
                     if(s(1:5).eq.'MODEL')then
                        read(s,*)du,model_num
                        if(model_num.gt.1)then
                           do while(.true.)
                              read(10,'(A500)',end=1013) s
                              if(s(1:6).eq.'ENDMDL')goto 1014
                           enddo
                        endif
                     endif
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                  endif
               endif
               if(s(1:4).eq.'ATOM')then !read 'ATOM' =======>
                  read(s(13:16),*)du4 !will remove space before 'CA'
                  if(du4.eq.'CA')then !read 'CA' ---------->
                     du1=s(27:27) !residue insertion tag
*******remove repeated altLoc atoms (this does not care inserted residues) -------->
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
                        read(s,'(a6,I5,a6,A3,a1,A1,I4,A1,a3,3F8.3)') !variable include space, different from read(*,*)
     &                       du,itmp,du,aanam,du,ch(i),mm(ic,i),
     &                       Cins(ic,i),du,x2(ic,i),y2(ic,i),z2(ic,i)
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
 1014          continue
            enddo               !<----end do-while
         else                   !<----mmCIF format,if(iform(1).eq.2)
            in=0                !number of entries
            do while (.true.)   !start do-while ---------->
               read(10,'(A500)',end=1013) s
               if(i.gt.0)then   !skip unuseful read to save time, suppose '#' appear before complex completed
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
                  if(du.eq.'_atom_site.label_seq_id')i_resi=in !1,2,3, 
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
                        if(m_complex.eq.0)then !monomer
                           if(i_mn.gt.1)then !sometimes it may have no i_mn
                              read(ctmp(i_mn),*)mn_t
                              if(mn_t.ne.mn(i)) goto 1013 !only read first model
                           endif
                           read(ctmp(i_ch),*)ch_t
                           if(ch_t.ne.ch(i)) goto 1013 !only read first chain
                           read(ctmp(i_ent),*)ent_t
                           if(ent_t.ne.ent(i)) goto 1013 !only read first entity
                        else    !dimer
                           if(i_mn.gt.1)then !sometimes it may have no i_mn
                              read(ctmp(i_mn),*)mn_t
                              if(mn_t.ne.mn(i)) goto 1016 !only read first model
                           endif
                        endif
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
                        read(ctmp(i_x),*)x2(ic,i)
                        read(ctmp(i_y),*)y2(ic,i)
                        read(ctmp(i_z),*)z2(ic,i)
                        read(ctmp(i_resi),*)mm(ic,i) !residue order, 4,5,6
                        read(ctmp(i_ins),*)Cins(ic,i) !residue insertion, A, ?
                        if(Cins(ic,i).eq.'?')Cins(ic,i)=' '
                        ss(ic,i)=ctmp(i_res) !residue name, 'GLY', for seq_ID
***   
                        i8=mm(ic,i)
                        nres1(ic,i8,ichar(du1))=
     &                       nres1(ic,i8,ichar(du1))+1 !nres1 only for check altLoc
***   
                        if(i.ge.nmax)goto 1013
                     endif      !<-----mk=1
                  endif         !<-----read 'CA'
               endif            !<==== read 'ATOM'
 1016          continue         !skip model_num >1
            enddo               !<----end do-while
         endif                  !if(iform(ic).eq.2)
 1013    continue
         close(10)
         
c-------convert 'GLY' to 'G', etc ------------>
         if(ic.eq.1)then
            nseqA=i
            do i=1,nseqA
               chA(i)=ch(i)
               nresA(i)=mm(ic,i)
               seqA(i)=ss(ic,i) !GLY
               ins1(i)=Cins(1,i)
               xa(i)=x2(1,i)
               ya(i)=y2(1,i)
               za(i)=z2(1,i)
               do j=-1,20
                  if(ss(1,i).eq.aa(j))then
                     seq1A(i)=slc(j) !G
                     goto 121
                  endif
               enddo
               seq1A(i)=slc(-1)
 121           continue
            enddo
         else
            nseqB=i
            do i=1,nseqB
               chB(i)=ch(i)
               nresB(i)=mm(ic,i)
               seqB(i)=ss(ic,i) !'GLY'
               ins2(i)=Cins(2,i)
               xb(i)=x2(2,i)
               yb(i)=y2(2,i)
               zb(i)=z2(2,i)
               do j=-1,20
                  if(ss(2,i).eq.aa(j))then
                     seq1B(i)=slc(j) !'G'
                     goto 122
                  endif
               enddo
               seq1B(i)=slc(-1)
 122           continue
            enddo
         endif
 1001 continue                  !ic=1,2 for file1 and file2
c^^^^^^^^^^^^^read input files completed ^^^^^^^^^^^^^^^^^^^^^^

c      do i=1,nseqA
c         write(*,*)i,xa(i),seqA(i)
c      enddo
c      do i=1,nseqB
c         write(*,*)i,xb(i),seqB(i)
c      enddo

******************************************************************
*     pickup the aligned residues:
******************************************************************
      if(m_complex.eq.0)then    !monomer
         k=0
         do i=1,nseqA
            do j=1,nseqB
               if(nresA(i).eq.nresB(j))then
                  if(ins1(i).eq.ins2(j))then
                     k=k+1
                     iA(k)=i
                     iB(k)=j
                     goto 205
                  endif
               endif
            enddo
 205        continue
         enddo
      else                      !complex
         k=0
         do i=1,nseqA
            do j=1,nseqB
               if(nresA(i).eq.nresB(j).and.chA(i).eq.chB(j))then
                  if(ins1(i).eq.ins2(j))then !Residue_insertion_code
                     k=k+1
                     iA(k)=i
                     iB(k)=j
                     goto 206
                  endif
               endif
            enddo
 206        continue
         enddo
      endif
      n_ali=k                   !number of aligned residues
      if(n_ali.lt.1)then
         write(*,*)'There is no common residues in the input structures'
         goto 9999
      endif
      
c      do i=1,n_ali
c         write(*,*)i,iA(i),iB(i)
c      enddo

******************************************************************
*     check the residue serial numbers ------------->
      if(m_complex.eq.0)then
         nA_repeat=0
         do i=1,nseqA
            do j=i+1,nseqA
               if(nresA(i).eq.nresA(j))then
                  if(ins1(i).eq.ins1(j))then
                     nA_repeat=nA_repeat+1
                  endif
               endif
            enddo
         enddo
         if(nA_repeat.gt.0)then
            write(*,380)nA_repeat
         endif
 380     format('Warning: TMscore calculation can be wrong because ',
     &        'there are ',I3,' residues with the serial number same ',
     &        'as others in first Chain. Please modify the PDB file ',
     &        'and rerun the program!!')
         nB_repeat=0
         do i=1,nseqB
            do j=i+1,nseqB
               if(nresB(i).eq.nresB(j))then
                  if(ins2(i).eq.ins2(j))then
                     nB_repeat=nB_repeat+1
                  endif
               endif
            enddo
         enddo
         if(nB_repeat.gt.0)then
            write(*,381)nB_repeat
         endif
 381     format('Warning: TMscore calculation can be wrong because ',
     &        'there are ',I3,' residues with the serial number same ',
     &        'as others in second Chain. Please modify the PDB file ',
     &        'and rerun the program!!')
      endif
      
************/////
*     parameters:
*****************
***   d0------------->
      if(nseqB.gt.15)then
         d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      else
         d0=0.5
      endif
      if(m_len.eq.1)then
         d0=1.24*(l0_fix-15)**(1.0/3.0)-1.8
      endif
      if(d0.lt.0.5)d0=0.5
      if(m_fix.eq.1)d0=d0_fix
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      if(m_fix.eq.1)d_output=d0_fix
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
      score_maxsub_max=-1       !MaxSub-score
      score10_max=-1            !TM-score10
      n_GDT05_max=-1            !number of residues<0.5
      n_GDT1_max=-1             !number of residues<1
      n_GDT2_max=-1             !number of residues<2
      n_GDT4_max=-1             !number of residues<4
      n_GDT8_max=-1             !number of residues<8
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
              ka=ka+1
              k_ali(ka)=k
              LL=LL+1
           enddo
           if(i_init.eq.1)then  !global superposition
              call u3b(w,r_1,r_2,LL,2,rms,u,t,ier) !0:rmsd; 1:u,t; 2:rmsd,u,t
              armsd=dsqrt(rms/LL)
              rmsd_ali=armsd
           else
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           endif
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
           if(score10_max.lt.score10)score10_max=score10
           if(score_maxsub_max.lt.score_maxsub)score_maxsub_max=
     &          score_maxsub
           if(n_GDT05_max.lt.n_GDT05)n_GDT05_max=n_GDT05
           if(n_GDT1_max.lt.n_GDT1)n_GDT1_max=n_GDT1
           if(n_GDT2_max.lt.n_GDT2)n_GDT2_max=n_GDT2
           if(n_GDT4_max.lt.n_GDT4)n_GDT4_max=n_GDT4
           if(n_GDT8_max.lt.n_GDT8)n_GDT8_max=n_GDT8
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
              if(score10_max.lt.score10)score10_max=score10
              if(score_maxsub_max.lt.score_maxsub)score_maxsub_max
     &             =score_maxsub
              if(n_GDT05_max.lt.n_GDT05)n_GDT05_max=n_GDT05
              if(n_GDT1_max.lt.n_GDT1)n_GDT1_max=n_GDT1
              if(n_GDT2_max.lt.n_GDT2)n_GDT2_max=n_GDT2
              if(n_GDT4_max.lt.n_GDT4)n_GDT4_max=n_GDT4
              if(n_GDT8_max.lt.n_GDT8)n_GDT8_max=n_GDT8
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
      
      ratio=1
      if(m_len.gt.0)then
         ratio=float(nseqB)/float(l0_fix)
      endif
      
******************************************************************
*     Output
******************************************************************
***   output TM-scale ---------------------------->
      write(*,*)
      write(*,*)'*****************************************************',
     &     '************************'
      write(*,*)'*                                 TM-SCORE           ',
     &     '                       *'
      write(*,*)'* A scoring function to assess the similarity of prot',
     &     'ein structures         *'
      write(*,*)'* Based on statistics:                               ',
     &     '                       *'
      write(*,*)'*       0.0 < TM-score < 0.17, random structural simi',
     &     'larity                 *'
      write(*,*)'*       0.5 < TM-score < 1.00, in about the same fold',
     &     '                       *'
      write(*,*)'* Reference: Yang Zhang and Jeffrey Skolnick, ',
     &     'Proteins 2004 57: 702-710     *'
      write(*,*)'* For comments, please email to: zhng@umich.edu      ',
     &     '                       *'
      write(*,*)'*****************************************************',
     &     '************************'
      write(*,*)
      write(*,501)pdb(1),nseqA
 501  format('Structure1: ',A10,'  Length= ',I4)
      if(m_len.eq.1)then
         write(*,411)pdb(2),nseqB
         write(*,412)l0_fix
      else
         write(*,502)pdb(2),nseqB
      endif
 411  format('Structure2: ',A10,'  Length= ',I4)
 412  format('TM-score is notmalized by ',I4)
 502  format('Structure2: ',A10,'  Length= ',I4,
     &     ' (by which all scores are normalized)')
      write(*,503)n_ali
 503  format('Number of residues in common= ',I4)
      write(*,513)rmsd_ali
 513  format('RMSD of  the common residues= ',F8.3)
      write(*,*)
      if(m_len.eq.1)then
         score_max=score_max*float(nseqB)/float(l0_fix)
      endif
      write(*,504)score_max,d0
 504  format('TM-score    = ',f6.4,'  (d0=',f5.2,')')
      write(*,505)score_maxsub_max*ratio
 505  format('MaxSub-score= ',f6.4,'  (d0= 3.50)')
      score_GDT=(n_GDT1_max+n_GDT2_max+n_GDT4_max+n_GDT8_max)
     &     /float(4*nseqB)
      write(*,506)score_GDT*ratio,n_GDT1_max/float(nseqB)*ratio,
     &     n_GDT2_max/float(nseqB)*ratio,n_GDT4_max/float(nseqB)*ratio,
     &     n_GDT8_max/float(nseqB)*ratio
 506  format('GDT-TS-score= ',f6.4,' %(d<1)=',f6.4,' %(d<2)=',f6.4,
     $     ' %(d<4)=',f6.4,' %(d<8)=',f6.4)
      score_GDT_HA=(n_GDT05_max+n_GDT1_max+n_GDT2_max+n_GDT4_max)
     &     /float(4*nseqB)
      write(*,507)score_GDT_HA*ratio,n_GDT05_max/float(nseqB)*ratio,
     &     n_GDT1_max/float(nseqB)*ratio,n_GDT2_max/float(nseqB)*ratio,
     &     n_GDT4_max/float(nseqB)*ratio
 507  format('GDT-HA-score= ',f6.4,' %(d<0.5)=',f6.4,' %(d<1)=',f6.4,
     $     ' %(d<2)=',f6.4,' %(d<4)=',f6.4)
      write(*,*)
      
***   recall and output the superposition of maxiumum TM-score:
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
         xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo

********* extract rotation matrix ------------>
      write(*,*)'-------- rotation matrix to rotate Chain-1 to ',
     &     'Chain-2 ------'
      write(*,*)'i          t(i)         u(i,1)         u(i,2) ',
     &     '        u(i,3)'
      do i=1,3
         write(*,304)i,t(i),u(i,1),u(i,2),u(i,3)
      enddo
c      do j=1,nseqA
c         xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
c         yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
c         zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
c         write(*,*)j,xt(j),yt(j),zt(j)
c      enddo
      write(*,*)
 304  format(I2,f18.10,f15.10,f15.10,f15.10)

********* rmsd in superposed regions --------------->
      d=d_output                !for output
      call score_fun()          !give i_ali(i), score_max=score now
      
******* for output ---------->
      if(m_out.eq.1)then
         do i=1,3
            t2(i)=t(i)
            do j=1,3
               u2(i,j)=u(i,j)
            enddo
         enddo
      endif
*******************************

***   record aligned residues by i=[1,nseqA], for sequenceM()------------>
      do i=1,nseqA
         iq(i)=0
      enddo
      do i=1,n_cut
         j=iA(i_ali(i))         ![1,nseqA]
         k=iB(i_ali(i))         ![1,nseqB]
         dis=sqrt((xt(j)-xb(k))**2+(yt(j)-yb(k))**2+(zt(j)-zb(k))**2)
c         write(*,*)i,j,k,dis,d_output,'1--'
         if(dis.lt.d_output)then
            iq(j)=1
         endif
      enddo
*******************************************************************
***   output aligned sequences
      if(m_complex.eq.0)then
         k=0
         i=1
         j=1
 800     continue
         if(i.gt.nseqA.and.j.gt.nseqB)goto 802
         if(i.gt.nseqA.and.j.le.nseqB)then
            k=k+1
            sequenceA(k)='-'
            sequenceB(k)=seq1B(j)
            sequenceM(k)=' '
            j=j+1
            goto 800
         endif
         if(i.le.nseqA.and.j.gt.nseqB)then
            k=k+1
            sequenceA(k)=seq1A(i)
            sequenceB(k)='-'
            sequenceM(k)=' '
            i=i+1
            goto 800
         endif
         if(nresA(i).eq.nresB(j))then
            k=k+1
            sequenceA(k)=seq1A(i)
            sequenceB(k)=seq1B(j)
            if(iq(i).eq.1)then
               sequenceM(k)=':'
            else
               sequenceM(k)=' '
            endif
            i=i+1
            j=j+1
            goto 800
         elseif(nresA(i).lt.nresB(j))then
            k=k+1
            sequenceA(k)=seq1A(i)
            sequenceB(k)='-'
            sequenceM(k)=' '
            i=i+1
            goto 800
         elseif(nresB(j).lt.nresA(i))then
            k=k+1
            sequenceA(k)='-'
            sequenceB(k)=seq1B(j)
            sequenceM(k)=' '
            j=j+1
            goto 800
         endif
 802     continue
      else
         k=0
         i=1                    !
         j=1
 803     continue
         if(i.gt.nseqA.and.j.gt.nseqB)goto 804
         if(i.gt.nseqA.and.j.le.nseqB)then
            k=k+1
            sequenceA(k)='-'
            sequenceB(k)=seq1B(j)
            sequenceM(k)=' '
            j=j+1
            goto 803
         endif
         if(i.le.nseqA.and.j.gt.nseqB)then
            k=k+1
            sequenceA(k)=seq1A(i)
            sequenceB(k)='-'
            sequenceM(k)=' '
            i=i+1
            goto 803
         endif
         if(chA(i).ne.chB(j))then
            kk=0                !check if chA(i) match later chains
            do i2=j,nseqB
               if(chA(i).eq.chB(i2))then
                  kk=kk+1
               endif
            enddo
            if(kk.eq.0)then     !chA(i) is not matched
               k=k+1
               sequenceA(k)=seq1A(i)
               sequenceB(k)='-'
               sequenceM(k)=' '
               i=i+1
               goto 803
            endif
            
            kk=0                !check if chB(j) match later chains
            do i1=i,nseqA
               if(chA(i1).eq.chB(j))then
                  kk=kk+1
               endif
            enddo
            if(kk.eq.0)then     !chB(j) is not matched
               k=k+1
               sequenceA(k)='-'
               sequenceB(k)=seq1B(j)
               sequenceM(k)=' '
               j=j+1
               goto 803
            endif
            write(*,*)'Chains in complexes have crossed order, ',
     &           'therefore, no alignment is output'
            stop
         else
            if(nresA(i).eq.nresB(j))then
               k=k+1
               sequenceA(k)=seq1A(i)
               sequenceB(k)=seq1B(j)
               if(iq(i).eq.1)then
                  sequenceM(k)=':'
               else
                  sequenceM(k)=' '
               endif
               i=i+1
               j=j+1
               goto 803
            elseif(nresA(i).lt.nresB(j))then
               k=k+1
               sequenceA(k)=seq1A(i)
               sequenceB(k)='-'
               sequenceM(k)=' '
               i=i+1
               goto 803
            elseif(nresB(j).lt.nresA(i))then
               k=k+1
               sequenceA(k)='-'
               sequenceB(k)=seq1B(j)
               sequenceM(k)=' '
               j=j+1
               goto 803
            endif
         endif
 804     continue
      endif
      
ccc   RMSD (d<5.0)-------->
      LL=0
      do i=1,n_cut
         m=i_ali(i)             ![1,nseqA]
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,0,rms,u,t,ier)
      armsd=dsqrt(rms/LL)
      rmsd=armsd

      write(*,600)d_output,n_cut,rmsd
 600  format('Superposition in the TM-score: Length(d<',f3.1,
     $     ')=',i3,'  RMSD=',f6.2)
      write(*,603)d_output
 603  format('(":" denotes the residue pairs of distance < ',f3.1,
     &     ' Angstrom)')
      write(*,601)(sequenceA(i),i=1,k)
      write(*,601)(sequenceM(i),i=1,k)
      write(*,601)(sequenceB(i),i=1,k)
      write(*,602)(mod(i,10),i=1,k)
 601  format(2000A1)
 602  format(2000I1)
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
         open(unit=10,file=pdb(ic),status='old')
         if(iform(ic).eq.1)then !file in PDB format %%%%%%%%%%------->
            do while (.true.)
               read(10,'(A500)',end=1015) s
               if(nL(ic).gt.0)then !decide n_cut
                  if(m_complex.eq.0)then
                     if(s(1:3).eq.'TER'.or.s(1:3).eq.'MOD'.
     &                    or.s(1:3).eq.'END')then
                        goto 1015 !no ligand allowed
                     endif
                  else
***   skip model_number >=2 ----------->
                     if(s(1:5).eq.'MODEL')then
                        read(s,*)du,model_num
                        if(model_num.gt.1)then
                           do while(.true.)
                              read(10,'(A500)',end=1015) s
                              if(s(1:6).eq.'ENDMDL')goto 1017
                           enddo
                        endif
                     endif
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                  endif
               endif
               if(s(1:6).eq.'ATOM  ')then !read ATOM ====>
*******remove repeated altLoc atoms -------->
                  mk=1
                  du1=s(27:27)  !read insertion tag
                  du3=s(22:22)  !chain ID
                  if(s(17:17).ne.' ')then !with Alternate atom
                     du2=s(13:16) !atom ID
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
                     read(s,8999)Agroup(ic,nL(ic)),Aatomi(ic,nL(ic)),
     &                    du,Aatom(ic,nL(ic)),Aalt(ic,nL(ic)),
     &                    Ares(ic,nL(ic)),du,Ach(ic,nL(ic)),
     &                    Aresi(ic,nL(ic)),Ains(ic,nL(ic)),du,
     &                    xx(ic,nL(ic)),yy(ic,nL(ic)),zz(ic,nL(ic))
***   
                     i8=Aresi(ic,nL(ic))
                     if(nres2(ic,i8,ichar(du1),ichar(du3)).lt.50)then
                        nres2(ic,i8,ichar(du1),ichar(du3))=
     &                       nres2(ic,i8,ichar(du1),ichar(du3))+1 !#of atoms for res_ins
                        atom1(nres2(ic,i8,ichar(du1),ichar(du3)))=
     &                       Aatom(ic,nL(ic)) !atom ID
                     endif
***   
                     if(nL(ic).ge.namax)goto 1015
                  endif         !<------mk=1
               endif            !<======read ATOM
 1017          continue
            enddo               !<--- do while
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
               if(s(1:4).eq.'ATOM')then !read 'ATOM' =======>
                  read(s,*)(ctmp(j),j=1,in)
                  if(nL(ic).gt.0)then
                     if(m_complex.eq.0)then !monomer
                        if(i_mn.gt.1)then !sometimes it may have no i_mn
                           read(ctmp(i_mn),*)mn_t
                           if(mn_t.ne.mn(i)) goto 1015 !only read first model
                        endif
                        read(ctmp(i_ch),*)ch_t
                        if(ch_t.ne.Ach(ic,nL(ic))) goto 1015
                        read(ctmp(i_ent),*)ent_t
                        if(ent_t.ne.Aent(ic,nL(ic))) goto 1015
                     else
                        if(i_mn.gt.1)then !sometimes it may have no i_mn
                           read(ctmp(i_mn),*)mn_t
                           if(mn_t.ne.mn(i)) goto 1018 !only read first model
                        endif
                     endif
                  endif
                  
*******remove repeated altLoc atoms -------->
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
                     
                     read(ctmp(i_res),*)Ares(ic,nL(ic))
                     read(ctmp(i_ch),*)Ach(ic,nL(ic)) !for check other chain
                     read(ctmp(i_ent),*)Aent(ic,nL(ic)) !for check other entity
                     read(ctmp(i_mn),*)mn(i) !for check model number
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
                     if(nL(ic).ge.namax)goto 1015
                  endif         !<---- mk=1
               endif            !<==== read 'ATOM'
 1018          continue
            enddo               !<----end do-while
         endif                  !<----mmCIF format,if(iform(1).eq.2) %%%%%%%%%
 1015    continue
         close(10)
 1002 continue                  !ic=1,2 for file1 and file2
      
 8999 format(a6,I5,a1,A4,a1,A3,a1,A1,I4,A1,a3,3F8.3)

c**************************************************************
c*************  Step-2: output 'aaa' (CA only)  ***************
c**************************************************************
      OPEN(unit=7,file=outname,status='unknown') !pdb1.aln + pdb2.aln

*************for pymol preset ------>
      xxx_p=outname(1:len_trim(outname))//'.pml' !for pymol script
      xxx_pdb=outname(1:len_trim(outname))//'.pdb' !for pymol script
      OPEN(unit=11,file=xxx_p,status='unknown')
      OPEN(unit=12,file=xxx_pdb,status='unknown')
      
      write(11,'(A,A)')'load ',xxx_pdb(1:len_trim(xxx_pdb))
      write(11,'(A)')'hide all'
      write(11,'(A)')'bg_color white'
      write(11,'(A)')'color blue, chain A'
      write(11,'(A)')'color red, chain B'
      write(11,'(A)')'set transparency=0.2'
      write(11,'(A)')'set stick_radius, 0.3'
      write(11,'(A)')'show sticks, chain A'
      write(11,'(A)')'show sticks, chain B'
*^^^^^^^^^ pymol preset complete ^^^^^^^^^^^^^^^^^^^
      
 900  format(A)
 901  format('select atomno= ',I5)
 902  format('select atomno= ',I5)
      write(7,900)'load inline'
      write(7,900)'select atomno<50001'
      write(7,900)'wireframe .45'
      write(7,900)'select atomno>50000'
      write(7,900)'wireframe .15'
      write(7,900)'select all'
      write(7,900)'color white' !all above atoms are white
      do i=1,n_cut
         write(7,901)iA(i_ali(i))
         write(7,900)'color red'
         write(7,902)50000+iB(i_ali(i))
         write(7,900)'color red'
      enddo
      write(7,900)'select all'
      write(7,900)'exit'
      write(7,514)rmsd_ali
 514  format('REMARK  RMSD of the common residues=',F8.3)
      write(7,515)score_max,d0
 515  format('REMARK  TM-score=',f6.4,' (d0=',f5.2,')')
      do i=1,nseqA
         write(7,1236)i,seqA(i),chA(i),nresA(i),ins1(i),
     &        xt(i),yt(i),zt(i)
         write(12,1236)i,seqA(i),'A',nresA(i),ins1(i),
     &        xt(i),yt(i),zt(i)
      enddo
      write(7,1238)
      write(12,1238)
      do i=2,nseqA
         write(7,1239)i-1,i
         write(12,1239)i-1,i
      enddo
      do i=1,nseqB
         write(7,1236)50000+i,seqB(i),chB(i),nresB(i),ins2(i),
     &        xb(i),yb(i),zb(i)
         write(12,1236)50000+i,seqB(i),'B',nresB(i),ins2(i),
     &        xb(i),yb(i),zb(i)
      enddo
      write(7,1238)
      write(12,1238)
      do i=2,nseqB
         write(7,1239)50000+i-1,50000+i
         write(12,1239)50000+i-1,50000+i
      enddo
      close(7)
      close(11)
      close(12)
 1236 format('ATOM  ',i5,'  CA  ',A3,' ',A1,I4,A1,3X,3F8.3)
 1238 format('TER')
 1239 format('CONECT',I5,I5)


c**************************************************************
c*************  Step-3: output 'aaa_atm' (all-atom) ***********
c**************************************************************
      outname=outname(1:len_trim(outname))//'_atm'
      OPEN(unit=7,file=outname,status='unknown') !pdb1.aln + pdb2.aln
      
*************for pymol preset ------>
      xxx_p=outname(1:len_trim(outname))//'.pml' !for pymol script
      xxx_pdb=outname(1:len_trim(outname))//'.pdb' !for pymol script
      OPEN(unit=11,file=xxx_p,status='unknown')
      OPEN(unit=12,file=xxx_pdb,status='unknown')
      
      write(11,'(A,A)')'load ',xxx_pdb(1:len_trim(xxx_pdb))
      write(11,'(A)')'hide all'
      write(11,'(A)')'bg_color white'
      write(11,'(A)')'color blue, chain A'
      write(11,'(A)')'color red, chain B'
      write(11,'(A)')'set transparency=0.2'
      write(11,'(A)')'set stick_radius, 0.3'
      write(11,'(A)')'show cartoon, chain A'
      write(11,'(A)')'show cartoon, chain B'
*^^^^^^^^^ pymol preset complete ^^^^^^^^^^^^^^^^^^^

      write(7,900)'load inline'
      write(7,900)'select atomno<50001'
      write(7,900)'color blue'
      write(7,900)'select atomno>50000'
      write(7,900)'color red'
      write(7,900)'select all'
      write(7,900)'cartoon'
      write(7,900)'exit'
      write(7,514)rmsd_ali
      write(7,515)score_max,d0
***   chain1:
      do i=1,nL(1)
         ax=t2(1)+u2(1,1)*xx(1,i)+u2(1,2)*yy(1,i)+u2(1,3)*zz(1,i)
         ay=t2(2)+u2(2,1)*xx(1,i)+u2(2,2)*yy(1,i)+u2(2,3)*zz(1,i)
         az=t2(3)+u2(3,1)*xx(1,i)+u2(3,2)*yy(1,i)+u2(3,3)*zz(1,i)

         write(7,8888)i,Aatom(1,i),Ares(1,i),Ach(1,i),
     &        Aresi(1,i),Ains(1,i),ax,ay,az
         write(12,8888)i,Aatom(1,i),Ares(1,i),'A',
     &        Aresi(1,i),Ains(1,i),ax,ay,az
      enddo
      write(7,1238)             !TER
      write(12,1238)             !TER
***   chain2:
      do i=1,nL(2)
         write(7,8888)50000+i,Aatom(2,i),Ares(2,i),Ach(2,i),
     &        Aresi(2,i),Ains(2,i),xx(2,i),yy(2,i),zz(2,i)
         write(12,8888)50000+i,Aatom(2,i),Ares(2,i),'B',
     &        Aresi(2,i),Ains(2,i),xx(2,i),yy(2,i),zz(2,i)
      enddo
      write(7,1238)             !TER
      write(12,1238)             !TER
      close(7)
      close(11)
      close(12)
 8888 format('ATOM  ',I5,1x,A4,1x,A3,' ',A1,I4,A1,3x,3F8.3)
      
*^^^^^^^^^^^^^^^^^^ output finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 9999 END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun
      PARAMETER(nmax=5000)

      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score,score_maxsub,score_fix,score10
      common/GDT/n_GDT05,n_GDT1,n_GDT2,n_GDT4,n_GDT8
      double precision score,score_max,score_fix,score_fix_max
      double precision score_maxsub,score10

      d_tmp=d
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      n_GDT05=0                 !for GDT-score, # of dis<0.5
      n_GDT1=0                  !for GDT-score, # of dis<1
      n_GDT2=0                  !for GDT-score, # of dis<2
      n_GDT4=0                  !for GDT-score, # of dis<4
      n_GDT8=0                  !for GDT-score, # of dis<8
      score_maxsub_sum=0        !Maxsub-score
      score_sum=0               !TMscore
      score_sum10=0             !TMscore10
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xt(i)-xb(j))**2+(yt(i)-yb(j))**2+(zt(i)-zb(j))**2)
***   for iteration:
         if(dis.lt.d_tmp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
***   for GDT-score:
         if(dis.le.8)then
            n_GDT8=n_GDT8+1
            if(dis.le.4)then
               n_GDT4=n_GDT4+1
               if(dis.le.2)then
                  n_GDT2=n_GDT2+1
                  if(dis.le.1)then
                     n_GDT1=n_GDT1+1
                     if(dis.le.0.5)then
                        n_GDT05=n_GDT05+1
                     endif
                  endif
               endif
            endif
         endif
***   for MAXsub-score:
         if(dis.lt.3.5)then
            score_maxsub_sum=score_maxsub_sum+1/(1+(dis/3.5)**2)
         endif
***   for TM-score:
         score_sum=score_sum+1/(1+(dis/d0)**2)
***   for TM-score10:
         if(dis.lt.10)then
            score_sum10=score_sum10+1/(1+(dis/d0)**2)
         endif
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
         goto 21
      endif
      score_maxsub=score_maxsub_sum/float(nseqB) !MAXsub-score
      score=score_sum/float(nseqB) !TM-score
      score10=score_sum10/float(nseqB) !TM-score10
      
      return
      end

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
      

