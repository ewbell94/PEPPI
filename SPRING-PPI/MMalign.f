*******************************************************************************
*     This program is to identify the best alignment of two protein complex
*     structures and output the alignment with best TM-score. By default,
*     TM-score is normalized by the second complex. For comments/suggestions,
*     please contact email: srayanta@umich.edu or zhng@umich.edu)    
*     
*     Reference to cite:
*     Srayanta Mukherjee, Yand Zhang. Nucleic Acids Research 2009; 37:e83 
*     
*     Permission to use, copy, modify, and distribute this program for 
*     any purpose, with or without fee, is hereby granted, provided that
*     the notices on the head, the reference information, and this
*     copyright notice appear in all copies or substantial portions of 
*     the Software. It is provided "as is" without express or implied 
*     warranty.
*********************** updating history **************************************
*    2009/01/27: Extension from dimer to multimer align
*    2009/09/30: A bug in the output display was fixed
*    2013/08/15: Fixed a number of compiler bugs
*******************************************************************************
      
      program compares
      PARAMETER(nmax=5000)
      PARAMETER(nmax2=10000)

      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/length/nseq1,nseq2
      common/length1/nseq11,nseq22,nseq11_p,nseq22_p
      common/d0/d0,anseq
      common/d0min/d0_min
      common/d00/d00,d002
      common/Tmsin/TMsin(nmax,nmax)
      character*100 fnam,pdb(100),outname,pdb1,pdb2
      character*3 aa(-1:20),aanam,ss1(nmax),ss2(nmax)
      character*100 s,du
      character*200 outnameall_tmp,outnameall
      character seq1(0:nmax),seq2(0:nmax),sequ
      character aseq1(nmax2),aseq2(nmax2),aseq3(nmax2)
      character*1 temp
      character*2 tmp
      dimension xx(nmax),yy(nmax),zz(nmax)
      dimension m1(nmax),m2(nmax)
      dimension mf1(nmax),mf2(nmax)
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      dimension xtmf1(nmax),ytmf1(nmax),ztmf1(nmax)
      dimension xtmf2(nmax),ytmf2(nmax),ztmf2(nmax)
      common/init/invmap_i(nmax)
      common/cutoff/int_cut
      common/TM/TM,TMmax,TM8temp
      common/n1n2/n1(nmax),n2(nmax)
      common/d8/d8
      common/initial4/mm1(nmax),mm2(nmax)
      common/initial5/inter1(nmax),inter2(nmax),inter3(nmax)
      common/initial5/inter4(nmax)
      common /lens/ len_first1,len_first2,len_second1,len_second2
      common/index/ind
ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
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
      common/choice/ichachoice,ichain,ii
      integer a(nmax2)
      !common npm
      common/arr/ib(1000000,50),ic(1000000,50),id(1000000,50)
      common/arr1/ia(1000000,50)
      common/chain/ibig,ismall
      COMMON/BACKBONE/iXA(3,nmax,0:1)   
      COMMON/BACKBONE/iXB(3,nmax,0:1)
      common/length/inseq1,inseq2
      common/initial4/imm1(nmax),imm2(nmax),imm3(nmax),imm4(nmax)        
      common/initial6/ifcounter1,ifcounter2
      common/intscore/n_int
      character*100 ipdb*(100)
      character*100 is,idu
      character*3 iaa(-1:20),iaanam,iss1(nmax),iss2(nmax),iss3(nmax)
      character*3 iss4(nmax)
      character iseq1(0:nmax),iseq2(0:nmax),iseq3(0:nmax),iseq4(0:nmax)
      dimension ixx(nmax),iyy(nmax),izz(nmax)
      dimension im1(nmax),im2(nmax),im3(nmax),im4(nmax)
      common/flag/intflag,domflag
      common/initial4/mempo(nmax)
      common/chains/ichainterm1(0:nmax),ichainterm2(0:nmax)
      common/nchains/itercount,itercount1,icombo1,icombo2 
      Logical mtc
      integer ii,jj,chk
      real*8 dist,xaa,yaa,zaa,cov
ccc

      data iaa/ 'BCK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP','CYX'/
      character*1 islc(-1:20)
      data islc/'X','G','A','S','C',
     &     'V','T','I','P','M',                                         

     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W','C'/

      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
         write(*,*)
         write(*,*)'Brief instruction for running MM-align program:'
         write(*,*)'(For detail: Mukherjee & Zhang, Nucl. Acid. Res.',
     &     ' 37:e83, 2009)'
         write(*,*)
         write(*,*)'Please make sure that the two chains of both '//
     &        'complexes are seperated by a "TER".'
         write(*,*)
         write(*,*)'1. Align ''complex1.pdb'' to ''complex2.pdb'''
         write(*,*)'  (By default, TM-score is normalized by the ',
     &        'length of ''complex2.pdb'')'
         write(*,*)'  >MMalign complex1.pdb complex2.pdb'
         write(*,*)
         write(*,*)'2. Run MM-align and output the superposition ',
     &        'to ''MM.sup'' and ''MM.sup_all'':'
         write(*,*)'  >MMalign complex1.pdb complex2.pdb -o MM.sup'
         write(*,*)'   To view the superimposed structures of the',
     &        ' aligned regions by rasmol:'
         write(*,*)'  >rasmol -script MM.sup'
         write(*,*)'   To view the superimposed structures of all',
     &        ' regions by rasmol:'
         write(*,*)'  >rasmol -script MM.sup_all'
         write(*,*)
         write(*,*)'3. If you want TM-score normalized by ',
     &        'an assigned length, e.g. 100 aa:'
         write(*,*)'  >MMalign complex1.pdb complex2.pdb -L 100'
         write(*,*)'   If you want TM-score normalized by the ',
     &        'average length of two structures:'
         write(*,*)'  >MMalign complex1.pdb complex2.pdb -a'
         write(*,*)'   If you want TM-score normalized by the ',
     &        'shorter length of two structures:'
         write(*,*)'  >MMalign complex1.pdb complex2.pdb -b'
         write(*,*)'   If you want TM-score normalized by the ',
     &        'longer length of two structures:'
         write(*,*)'  >MMalign complex1.pdb complex2.pdb -c'
         write(*,*)
         write(*,*)'4. If you want MM-align to enforce interface'//
     &         ' alignment defined by a cutoff e.g. 10A :'
         write(*,*)'  >MMalign complex1.pdb complex2.pdb -I 10' 
         write(*,*)'   (Distance cutoff is used to determine interface',
     &        ' resdiues. Default 8A.)'
         write(*,*)
         write(*,*)'5. If you want MM-align to allow cross'//
     &        ' alignment e.g: Domain swap, Gene Fusion :'
         write(*,*)'  >MMalign complex1.pdb complex2.pdb -DG ' 
  
         write(*,*)
         write(*,*)'(First three options,do not change the',
     &         ' final structure alignment result)'
         write(*,*)
         goto 9999
      endif
 
******* options ----------->
      m_out=-1                  !decided output
      m_fix=-1                  !fixed length-scale only for output
      m_ave=-1                  !using average length
      m_d0_min=-1               !diminum d0 for search
      m_d0=-1                   !given d0 for both search and output
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
      elseif(fnam.eq.'-L')then  !change both L_all and d0
         m_fix=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)L_fix
      elseif(fnam.eq.'-dmin')then
         m_d0_min=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)d0_min_input
      elseif(fnam.eq.'-d0')then
         m_d0=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)d0_fix
      elseif(fnam.eq.'-a')then ! this will change the superposed output but not the alignment
         m_ave=1
         i=i+1
      elseif(fnam.eq.'-b')then
        m_ave=2
         i=i+1
      elseif(fnam.eq.'-c')then
        m_ave=3
         i=i+1
      elseif(fnam.eq.'-I')then
         intflag=1
         i=i+1
c         m_out=1
         call getarg(i,fnam)
         if(fnam.eq. ' ')then
            int_cut=8
         else
            write(*,*)fnam
            read(fnam,*)int_cut
         endif
      elseif(fnam.eq.'-DG')then
         dom_flag=1         
      else
         j=j+1
         pdb(j)=fnam
      endif
      if(i.lt.narg)goto 115
      
ccccc Searching for interface residues *************      
      if (intflag.eq.1)then
         pdb1=pdb(1)
         pdb2=pdb(2)
         call interface(pdb1,pdb2,int_cut)
      endif     
ccccc Determine no of chains ccccccccc
      ichainterm1(0)=1
      ichainterm2(0)=1
      open(unit=10,file=pdb(1),status='old')
      i=0
      itercount=0
      do while (.true.)
         read (10,9001,end=1234) s
         if (s(1:3).eq.'TER')then
            itercount=itercount+1
            ichainterm1(itercount)=i+1 !-ichainterm1(itercount-1)
         endif
         if(s(1:3).eq.'ATO')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or
     &           .s(13:16).eq.'  CA')then
               if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
                  i=i+1
               endif
            endif
         endif
      enddo
 1234 continue
      
      open(unit=10,file=pdb(2),status='old')
      i=0
      itercount1=0
      do while (.true.)
         read (10,9001,end=1235) s
         if (s(1:3).eq.'TER')then
            itercount1=itercount1+1
            ichainterm2(itercount1)=i+1
         endif
         if(s(1:3).eq.'ATO')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or
     &           .s(13:16).eq.'  CA')then
               if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
                  i=i+1
               endif
            endif
         endif
      enddo
 1235 continue
ccccccccc read data from first CA file:
      temp=""
      open(unit=10,file=pdb(1),status='old')
      i=0
      l1=0
      flag=0
      terflag=0
      counter=1
      do while (.true.)
         read(10,9001,end=1010) s
         if(i.gt.0.and.s(1:3).eq.'TER'.and.terflag.eq.itercount)then
            goto 1010
         endif
         if(i.gt.0.and.s(1:3).eq.'TER'.and.terflag.lt.itercount)then
              terflag=terflag+1
              goto 1017 
         endif
         if(s(1:3).eq.'ATO')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or
     &           .s(13:16).eq.'  CA')then
               if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
                  i=i+1
                  read(s,9000)du,aanam,du,mm1(i),du,
     $                 xa(1,i,0),xa(2,i,0),xa(3,i,0)
                  do j=-1,20
                     if(aanam.eq.aa(j))then
                        seq1(i)=slc(j)
                        goto 21
                     endif
                  enddo
                  seq1(i)=slc(-1)
 21               continue
                  ss1(i)=aanam
                  if(i.ge.nmax)goto 1010
               endif
               if(s(10:11) .eq. " 2" .and. flag .eq. 0)then
                    temp=s(22:22)
                    flag=1
               endif
            endif
         endif
 1017  continue
      enddo
 1010 continue
 9000 format(A17,A3,A2,i4,A4,3F8.3)
 9001 format(A100)
      close(10)
      nseq1=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
ccccccccc read data from the second CA file:
      open(unit=10,file=pdb(2),status='old')
      i=0
      l2=0
      flag1=0
      terflag1=0
      do while (.true.)
         read(10,9001,end=1011) s
         if(i.gt.0.and.s(1:3).eq.'TER'.and.terflag1.eq.itercount1)then
            goto 1011
         endif
         if(i.gt.0.and.s(1:3).eq.'TER'.and.terflag1.lt.itercount1)then
            terflag1=terflag1+1
            goto 1018
         endif
         if(s(1:3).eq.'ATO')then
            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.
     &           s(13:16).eq.'  CA')then
               if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
                  i=i+1
                  read(s,9000)du,aanam,du,mm2(i),du,
     $                 xa(1,i,1),xa(2,i,1),xa(3,i,1)
                  do j=-1,20
                     if(aanam.eq.aa(j))then
                        seq2(i)=slc(j)
                        goto 22
                     endif
                  enddo
                  seq2(i)=slc(-1)
 22               continue
                  ss2(i)=aanam
                  if(i.ge.nmax)goto 1011
               endif
            if(s(10:11) .eq. " 2" .and. flag1 .eq. 0)then
                  temp=s(22:22)
                  flag1=1
               endif
            endif
         endif
 1018 continue
      enddo
 1011 continue
      close(10)
      nseq2=i
      
ccccccError message in case there is only 1 chain 
      if (terflag1.gt.50.and.terflag.gt.50)then
         write(*,*)
         write(*,*)'*       Error. Too many chains in Complex files.'//
     &        ' Please split chains and try again        *'
         write(*,*)
         goto 9999
      endif
      if (terflag.gt.50)then
         write(*,*)
         write(*,*)'*    Error: Too many chains in complex 1.'//
     &        ' Please split chains and try again       *'
         write(*,*)
         goto 9999
      endif
      if (terflag1.gt.50)then
         write(*,*)
         write(*,*)'*    Error: Too many chains in complex 2.'//
     &        ' Please split chains and try again       *'
         write(*,*)
         goto 9999
      endif
      
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*!!!  Scale of TM-score in search is based on the smaller protein --------->
      d0_min=0.5
      if(m_d0_min.eq.1)then
         d0_min=d0_min_input    !for search
      endif
      anseq_min=min(nseq1,nseq2)
      anseq=anseq_min           !length for defining TMscore in search
      d8=1.5*anseq_min**0.3+3.5 !remove pairs with dis>d8 during search & final
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      if(m_d0.eq.1)d0=d0_fix
      d00=d0                    !for quickly calculate TM-score in searching
      if(d00.gt.8)d00=8
      if(d00.lt.4.5)d00=4.5
      d002=d00**2
      nseq=max(nseq1,nseq2)
      do i=1,nseq
         n1(i)=i
         n2(i)=i
      enddo
******** Quickly calculate TM-score FOR ALL CHAINS INDIVIDUALLY *********
      ichainterm1(0)=1
      ichainterm2(0)=1
      do ii=1,itercount
         do jj=1,itercount1
            nseq11=ichainterm1(ii)-ichainterm1(ii-1)
            nseq22=ichainterm2(jj)-ichainterm2(jj-1)
            nseq11_p=ichainterm1(ii-1)
            nseq22_p=ichainterm2(jj-1)
            !print *,nseq11,nseq22
**** RUN INITIAL1 TO GET TM-SCORE
            d01=d0+1.5
            if(d01.lt.d0_min)d01=d0_min
            d02=d01*d01
            aL=min(nseq11,nseq22)
            idel=aL/2.5         !minimum size of considered fragment
            if(idel.le.5)idel=5
            n11=-nseq22+idel
            n22=nseq11-idel
            GL_max=0
            !print *,nseq22_p,nseq22+nseq22_p
            do ishift=n11,n22
               L=0
               K=0
               nseq=nseq22+nseq22_p
               do j=nseq22_p,nseq
                  k=k+1
                  i=k+ishift
                  if(i.ge.1.and.i.le.nseq11)then
                     L=L+1
                     invmap(j)=i+nseq11_p
                  else
                     invmap(j)=-1
                  endif
               enddo
               !print *,L
               if(L.ge.idel)then
                  call get_GL11(GL)
                  !print *,GL,GL_max
                  if(GL.gt.GL_max)then
                     GL_max=GL
                     do i=nseq22_p,nseq22+nseq22_p
                        invmap_i(i)=invmap(i)
                     enddo
                  endif
               endif
            enddo
            GL_max=GL_max/nseq22
            TMsin(ii,jj)=GL_max
            !stop
         enddo
         !stop
      enddo
      !stop
      
      
*********  iterating over all chain combos *****************
      if(itercount.le.itercount1)then
         ibig = itercount1
         ismall = itercount
      else
         ibig = itercount
         ismall = itercount1 
      endif
      
*************************************************************
      mtc=.false.
      icomb1=0
      
 18   call nexper(ismall,a,mtc,even)
      if(a(1).ne.0)then
         icomb1=icomb1+1
         do i=1,ismall
            ia(icomb1,i)=a(i)
         enddo
      endif
      if(icomb1.eq.2)then
c     write(*,*) " Too many Possible Combinations of Complex 1"//
c     &        " Only First 1 million are tried out"
         goto 19
      endif
      if(mtc) goto 18
      
 19   continue
cccc  the smaller of the two proteins will be kept fixed
cccc  and various chain combinations of the bigger to be tried
      
***   Generating combinations of the bigger chain: nPm ****
***   where n = no of chains in bigger proteins        ****
***         m = no of chains in smaller proteins       ****
***********************************************************
      
      mtc=.false.
      icomb=0 

 17   call nexper(ibig,a,mtc,even)
      if(a(1).ne.0)then
         icomb=icomb+1
         do i=1,ibig
            ib(icomb,i)=a(i)
         enddo
      endif
     
      if(icomb.eq.1000)then
cc         write(*,*) " Too many Possible Combinations of Complex 2"//
c     &     " Only First 1 million are tried out"
c         icomb=2
         goto 16   
      endif
      if(mtc) goto 17
 16   continue
      icnt=0
      
      do i=1,icomb
         do j=1,ismall
            ic(i,j)=ib(i,j)
         enddo
     
         do ii=0,i-1
            do jj=1,ismall
               if(ic(i,jj).eq.ic(ii,jj))then
                  icnt=icnt+1
               endif
            enddo
            
            if(icnt.eq.ismall)then
              
               do jj=1,ismall
                  ic(i,jj)=0
               enddo
               icnt=0
               goto 117
            endif
            icnt=0
         enddo
 117     continue
         icnt=0
      enddo
      
      ind=0
      do i=1,icomb
         if(ic(i,1).ne.0)then
            ind=ind+1
            do j=1,ismall
               id(ind,j)=ic(i,j)
            enddo
         endif
      enddo
               
           
***   calculating npm *******
      
      !icomb=icomb+1
      fact1=1
      fact2=1
      fact3=1
      fnpm=1
      do k=2,ibig
         fact1=fact1*k
      enddo
      do k=2,ismall
         fact2=fact2*k
      enddo
      icdiff=ibig-ismall
      if(icdiff.le.1)then
         fact3=1
      else
         do k=2,icdiff
            fact3=fact3*k
         enddo
      endif
      if(ibig.eq.ismall)then
         fnpm=fact1
        ! print *,"here"
      else
        ! print *,"no here"
         fnpm=real(fact1)/real(fact3)
      endif
      tempTM=0
      n_al=0 
      itemp=0
      rmsdtp=0
      seqtp=0
      !print *,ibig,ismall,fact1,fact3,fnpm
      !stop
      if(fnpm.gt.100)then
         fnpm=100
      endif
      iss=0
      ibb=0
      sTM=0
      sTMmax=0
      iteration=0
      
      icombo1=1
      do icombo2=1,fnpm
c         write(*,*)icombo2

         ichain=0
         sTM=0
         do i=1,ismall
            iss=ia(icombo1,i)
            ibb=id(icombo2,i)
            sTM=TMsin(iss,ibb)+sTM
         enddo
                                !print *,sTM,sTMmax
         if(sTM.lt.sTMmax*0.95)goto 9998
         if(sTM.ge.sTMmax)then
            sTMmax=sTM
         endif
         
         iteration=iteration + 1
         if(fnpm.gt.5.and.iteration.eq.5)goto 9997
         do j=1,nseq1
            xtm1(j)=0
            ytm1(j)=0
            ztm1(j)=0
         enddo
         do j=1,nseq2
            xtm2(j)=0
            ytm2(j)=0
            ztm2(j)=0
         enddo
         
***** do alignment **************************
         CALL super_align       !to find invmap(j)
      !enddo
************************************************************
***   resuperpose to find residues of dis<d8 ------------------------>
         
         in_align=0
         n_int=0
         n_int1=0
         n_int2=0
                                !itemp=0
         ix=0
         iy=0
         cov=0
         n_al=0
         
         do j=1,nseq2
            if(invmap0(j).gt.0)then
               i=invmap0(j)
c               print *,i,j
               n_al=n_al+1
               xtm1(n_al)=xa(1,i,0)
               ytm1(n_al)=xa(2,i,0)
               ztm1(n_al)=xa(3,i,0)
               xtm2(n_al)=xa(1,j,1)
               ytm2(n_al)=xa(2,j,1)
               ztm2(n_al)=xa(3,j,1)
               m1(n_al)=i       !for recording residue order
               m2(n_al)=j
            endif
         enddo
                                !stop
         d0_input=d0
         call TMscore8(d0_input,n_al,xtm1,ytm1,ztm1,n1,n_al,
     &        xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) !TM-score with dis<d8 only
         
*!!!  Output TM-score is based on the second protein------------------>
         d0_min=0.5             !for output
         anseq=nseq2            !length for defining final TMscore
         if(m_ave.eq.1)anseq=(nseq1+nseq2)/2.0 !<L>
         if(m_ave.eq.2)anseq=min(nseq1,nseq2)
         if(m_ave.eq.3)anseq=max(nseq1,nseq2)
         if(anseq.lt.anseq_min)anseq=anseq_min
         if(m_fix.eq.1)anseq=L_fix !input length
         if(anseq.gt.15)then
            d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
         else
            d0=d0_min
         endif
         if(d0.lt.d0_min)d0=d0_min
         if(m_d0.eq.1)d0=d0_fix
***   calculating interface aligned residues for final report---->
         if(intflag.eq.1)then
            iy=0
            ix=0
            i=0
            icount1=0
            icount2=0
            dis2=0
            icount1=ifcounter1*2
            icount2=ifcounter2*2
            do j=0,len_first2
               if(invmap0(j).gt.0)then
                  i=invmap0(j)
                  if(inter1(i).eq.1 .and. inter3(j).eq.1)then
                     dis2=sqrt((xa(1,i,0)-xa(1,j,1))**2+(xa(1,i,0)-
     &                    xa(2,j,1))**2+(xa(3,i,0)-xa(3,j,0))**2)
                     iy=iy+1
                     if(dis2.le.d8)then
                        ix=ix+1
                     endif
                  endif
               endif
            enddo
            dis2=0
            i=0
            do j=len_first2,nseq2
               if(invmap0(j).gt.0)then
                  i=invmap0(j)
                  if(inter2(i).eq.1 .and. inter4(j).eq.1)then
                     dis2=sqrt((xa(1,i,0)-xa(1,j,1))**2+(xa(1,i,0)-
     &                    xa(2,j,1))**2+(xa(3,i,0)-xa(3,j,0))**2)
                     iy=iy+1
                     if(dis2.le.d8)then
                        ix=ix+1
                     endif
                  endif
               endif
            enddo
            n_int=iy            !n_int1+n_int2
            cov=real(iy)/real(icount1)
            n_int5=ix
            cov5=real(ix)/real(icount1)
         endif
         
***   remove dis>d8 in normal TM-score calculation for final report----->
         j=0
         n_eq=0
         do i=1,n_al
            dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &           (ztm1(i)-ztm2(i))**2)
            if(dis2.le.d8)then
               j=j+1
               xtm1(j)=xtm1(i)
               ytm1(j)=ytm1(i)
               ztm1(j)=ztm1(i)
               xtm2(j)=xtm2(i)
               ytm2(j)=ytm2(i)
               ztm2(j)=ztm2(i)
               m1(j)=m1(i)
               m2(j)=m2(i)
               if(ss1(m1(i)).eq.ss2(m2(i)))then
                  n_eq=n_eq+1
               endif
            endif
         enddo
         
         n8_al=j
         seq_id=float(n_eq)/(n_al+0.00000001)
         d0_input=d0
         call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n1,n8_al,
     &        xtm2,ytm2,ztm2,n2,TM8,Rcomm,Lcomm) !normal TMscore
         rmsd8_al=Rcomm
         TM8temp=TM8*n8_al/anseq
         if(TM8temp.gt.tempTM)then
            tempTM=TM8temp
            itemp=n8_al
            rmsdtp=rmsd8_al
            seqtp=seq_id
            do i=1,itemp
               xtmf1(i)=xtm1(i)
               ytmf1(i)=ytm1(i)
               ztmf1(i)=ztm1(i)
               xtmf2(i)=xtm2(i)
               ytmf2(i)=ytm2(i)
               ztmf2(i)=ztm2(i)
               mf1(i)=m1(i)
               mf2(i)=m2(i)
            enddo
         endif
 9998    continue
      enddo
 9997 continue                  !enddo
      TM8=tempTM
      rmsd8_al=rmsdtp

      n8_al=itemp
      seq_id=seqtp

********* for output summary ******************************
      write(*,*)
      write(*,*)'*****************************************************',
     &     '*********************'
      write(*,*)'*                               MM-align             ',
     &     '                    *'
      write(*,*)'*       Aligning protein complex structure by Dynamic',
     &     ' Programming        *'
      write(*,*)'*       Comments on the program, please email to: z',
     &     'hng@umich.edu         *'
      write(*,*)'*****************************************************',
     &     '*********************'
      write(*,*)
      write(*,101)pdb(1),nseq1
 101  format('Protein 1:',A10,'  Size=',I4)
      write(*,102)pdb(2),nseq2,int(anseq)
 102  format('Protein 2:',A10,'  Size=',I4,
     &     ' (TM-score is normalized by ',I4,')')
      write(*,*)
      if (intflag.eq.1)then
      write(*,111)du,ifcounter1,ifcounter1
 111  format('No. Of Interface Residue:'A4,'Protein 1,Chain 1 :',I4
     &       ' | Protein 1,Chain 2 :',I4)
      write(*,112)du,ifcounter2,ifcounter2
 112  format(A29,'Protein 2,Chain 1 :',I4,' | Protein 2,Chain 2 :',I4)
      write(*,*)
      write(*,113)n8_al,rmsd8_al,TM8,seq_id
 113  format('Aligned length=',I4,
     &        ', RMSD=',f6.2,', TM-score=',f7.5,', ID=',f5.3)
      write(*,*)
      else
      write(*,103)n8_al,rmsd8_al,TM8,seq_id
 103  format('Aligned length=',I4,
     &        ', RMSD=',f6.2,', TM-score=',f7.5,', ID=',f5.3)
      write(*,*)
      endif
********* extract rotation matrix ------------>
      L=0
      do i=1,n8_al
         k=mf1(i)
         L=L+1
         r_1(1,L)=xa(1,k,0)
         r_1(2,L)=xa(2,k,0)
         r_1(3,L)=xa(3,k,0)
         r_2(1,L)=xtmf1(i)
         r_2(2,L)=ytmf1(i)
         r_2(3,L)=ztmf1(i)
       enddo
       if(L.gt.3)then
          call u3b(w,r_1,r_2,L,2,rms,u,t,ier) !u rotate r_1 to r_2
          armsd=dsqrt(rms/L)
          write(*,*)'-------- rotation matrix to rotate Chain-1 to ',
     &         'Chain-2 ------'
          write(*,*)'i          t(i)         u(i,1)         u(i,2) ',
     &         '        u(i,3)'
          do i=1,3
             write(*,204)i,t(i),u(i,1),u(i,2),u(i,3)
          enddo
          write(*,*)
      endif
 204  format(I2,f18.10,f15.10,f15.10,f15.10)
      
********* for output superposition ******************************
      if(m_out.eq.1)then
 1237    format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
 1238    format('TER')
 1239    format('CONECT',I5,I5)
 900     format(A)
 901     format('select atomno=',I4)
 104     format('REMARK Chain 1:',A10,'  Size=',I4)
 105     format('REMARK Chain 2:',A10,'  Size=',I4,
     &        ' (TM-score is normalized by ',I4,')')
 106     format('REMARK Aligned length=',I4,', RMSD=',f6.2,
     &        ', TM-score=',f7.5,', ID=',f5.3)
         OPEN(unit=7,file=outname,status='unknown') !pdb1.aln + pdb2.aln
***   script:
         write(7,900)'load inline'
         write(7,900)'select atomno<4000'
         write(7,900)'wireframe .45'
         write(7,900)'select none'
         write(7,900)'select atomno>4000'
         write(7,900)'wireframe .15'
         write(7,900)'color white'
         write(7,900)'background white'
         ikk1=0
         ikk2=0
         il1=0
         il2=0
         do i=1,n8_al
            if(mf1(i).ge.ichainterm1(ikk1))then
               ikk1=ikk1+1
            endif
            il1=mod(ikk1,2)
            if(mf2(i).ge.ichainterm2(ikk2))then
               ikk2=ikk2+1
            endif
            il2=mod(ikk2,2)
            
            dis2=sqrt((xtmf1(i)-xtmf2(i))**2+
     &           (ytmf1(i)-ytmf2(i))**2+(ztmf1(i)-ztmf2(i))**2)
            if(dis2.le.5)then
               if(mf1(i).ge.ichainterm1(ikk1-1).and.mf1(i).lt.
     &              ichainterm1(ikk1).and.il1.ne.0)then
                  write(7,901)mf1(i)
                  write(7,900)'color red'
               else
                  write(7,901)mf1(i)
                  write(7,900)'color blue'
               endif
               if(mf2(i).ge.ichainterm2(ikk2-1).and.mf2(i).lt.
     &              ichainterm2(ikk2).and.il2.ne.0)then
                  write(7,901)4000+mf2(i)
                  write(7,900)'color cyan'
               else
                  write(7,901)4000+mf2(i)
                  write(7,900)'color magenta'
               endif
            endif
         enddo
         write(7,900)'select all'
         write(7,900)'exit'
         write(7,104)pdb(1),nseq1
         write(7,105)pdb(2),nseq2,int(anseq)
         write(7,106)n8_al,rmsd8_al,TM8,seq_id
***   chain1:
         ik1=1
         ik2=1
         do i=1,n8_al
            write(7,1237)mf1(i),ss1(mf1(i)),mm1(mf1(i)),
     &           xtmf1(i),ytmf1(i),ztmf1(i)
         enddo
         write(7,1238)          !TER
         do i=2,n8_al
            if(mf1(i-1).le.ichainterm1(ik1).and.
     &           mf1(i).gt.ichainterm1(ik1))then
               ik1=ik1+1
            else
               write(7,1239)mf1(i-1),mf1(i) !connect atoms
            endif
 1414       continue
         enddo
***   chain2:
         do i=1,n8_al
            write(7,1237)4000+mf2(i),ss2(mf2(i)),mm2(mf2(i)),
     $           xtmf2(i),ytmf2(i),ztmf2(i)
         enddo
         write(7,1238)
         do i=2,n8_al
            if(mf2(i-1).le.ichainterm2(ik2).and.
     &           mf2(i).gt.ichainterm2(ik2))then
               ik2=ik2+1
            else
               write(7,1239)4000+mf2(i-1),4000+mf2(i)
            endif
 1515       continue
         enddo
         close(7)
ccc   
         k=0
         outnameall_tmp=outname//'_all'
         outnameall=''
         do i=1,200
            if(outnameall_tmp(i:i).ne.' ')then
               k=k+1
               outnameall(k:k)=outnameall_tmp(i:i)
            endif
         enddo
         OPEN(unit=8,file=outnameall,status='unknown') !pdb1.aln + pdb2.aln
***   script:
         write(8,900)'load inline'
         write(8,900)'select atomno<4000'
         write(8,900)'wireframe .45'
         write(8,900)'select none'
         write(8,900)'select atomno>4000'
         write(8,900)'wireframe .15'
         write(8,900)'background white'
         ikk1=0
         ikk2=0
         il1=0
         il2=0
         do i=1,nseq1
            if(i.ge.ichainterm1(ikk1))then
               ikk1=ikk1+1
            endif
            il1=mod(ikk1,2)
            if(i.ge.ichainterm1(ikk1-1).and.i.lt.
     &              ichainterm1(ikk1).and.il1.ne.0)then
               write(8,901)i
               write(8,900)'color red'
            else
               write(8,901)i
               write(8,900)'color blue'
            endif
         enddo
         do j=1,nseq2
            if(j.ge.ichainterm2(ikk2))then
               ikk2=ikk2+1
            endif
            il2=mod(ikk2,2)
            if(j.ge.ichainterm2(ikk2-1).and.j.lt.
     &              ichainterm2(ikk2).and.il2.ne.0)then
               write(8,901)4000+j
               write(8,900)'color cyan'
            else
               write(8,901)4000+j
               write(8,900)'color magenta'
            endif
         enddo
         write(8,900)'select all'
         write(8,900)'exit'
         write(8,104)pdb(1),nseq1
         write(8,105)pdb(2),nseq2,int(anseq)
         write(8,106)n8_al,rmsd8_al,TM8,seq_id
***   chain1:
         ik1=1
         ik2=1
         
         do i=1,nseq1
            ax=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
            ay=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
            az=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
            write(8,1237)i,ss1(i),mm1(i),ax,ay,az
         enddo
         write(8,1238)          !TER
         do i=2,nseq1
            ij=i-1
            if(ij.eq.ichainterm1(ik1)-1)then
               !print *,i,ichainterm1(ik1)
               ik1=ik1+1
            else
               write(8,1239)i-1,i
            endif
 1212       continue
         enddo
***   chain2:
         do i=1,nseq2
            write(8,1237)4000+i,ss2(i),mm2(i),
     &           xa(1,i,1),xa(2,i,1),xa(3,i,1)
         enddo

         write(8,1238)
         do i=2,nseq2
            ij=i-1
            if(ij.eq.ichainterm2(ik2)-1)then
               !print *,i,ichainterm2(ik2)
               ik2=ik2+1
            else
               write(8,1239)4000+i-1,4000+i
            endif
 1313       continue
         enddo
         close(8)
      endif
*^^^^^^^^^^^^^^^^^^ output finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

************  output aligned sequences **************************
      ii=0
      i1_old=1
      i2_old=1
      itimes=0
      iflag1=0
      iflag2=0
      iflag3=0
      iflag4=0
      iflag5=0
      iflag6=0
      l=0
      !k1=0
      !k2=0
      ll=0
      l1=0
      l2=0
      ik1=1
      ik2=1
      ll1=0
      ll2=0
      
      do i=1,n8_al
                 
         if(i1_old-1.lt.ichainterm1(ik1).and.
     &        mf1(i)-1.ge.ichainterm1(ik1)-1)then
            itill1=ichainterm1(ik1)-1
            iflag1=1
         else
            itill1=mf1(i)-1
            iflag1=0
         endif

         if(i2_old-1.lt.ichainterm2(ik2).and.
     &        mf2(i)-1.ge.ichainterm2(ik2)-1)then
            itill2=ichainterm2(ik2)-1
            iflag2=1
         else
            itill2=mf2(i)-1
            iflag2=0
         endif
        
         do j=i1_old,itill1 !mf1(i)-1
            
            if(j.ge.ichainterm1(ik1))then
               ik1=ik1+1
               goto 7777
            endif
            ii=ii+1
            l1=mod(ik1,2)
            if(j.ge.ichainterm1(ik1-1).and.j.lt.ichainterm1
     &          (ik1).and. l1.ne.0)then
               isequence1=ichar(seq1(j))
               isequence1=isequence1+32
               aseq1(ii)=char(isequence1)
            else
               aseq1(ii)=seq1(j)
            endif
            aseq2(ii)='-'
            aseq3(ii)=' '
         enddo
 7777    continue
         do j=i2_old,itill2  !mf2(i)-1
            if(j.ge.ichainterm2(ik2))then
               ik2=ik2+1
               !print *,"going out"
               !goto 7774
            endif
            ii=ii+1
            aseq1(ii)='-'
            l2=mod(ik2,2)
            
            if(j.ge.ichainterm2(ik2-1).and.j.lt.ichainterm2(ik2).and.
     &           l2.ne.0)then
               isequence=ichar(seq2(j))
               isequence=isequence+32
               aseq2(ii)=char(isequence)
            else
               aseq2(ii)=seq2(j)
            endif
            aseq3(ii)=' '
         enddo

 7774    continue
 7776    continue
         if(iflag1.eq.1)then
            do j=itill1+1,mf1(i)-1 
               iflag3=0
               if(j.ge.ichainterm1(ik1))then
                  ik1=ik1+1
                  itill1=j-1
                  iflag3=1
                  goto 7773
               endif
               ii=ii+1
               l1=mod(ik1,2)
               
               if(j.ge.ichainterm1(ik1-1).and.j.lt.ichainterm1
     &              (ik1).and. l1.ne.0)then
                  isequence1=ichar(seq1(j))
                  isequence1=isequence1+32
                  aseq1(ii)=char(isequence1)
               else
                  aseq1(ii)=seq1(j)
               endif
               aseq2(ii)='-'
               aseq3(ii)=' '
            enddo
                                !iflag1=0
         endif
         iflag1=0
 7773    continue
 7778    continue
         if(iflag2.eq.1)then
            do j=itill2+1,mf2(i)-1
               iflag4=0
               if(j.ge.ichainterm2(ik2))then
                  ik2=ik2+1
                  itill2=j-1
                  iflag4=1
                  goto 7775
               endif
               ii=ii+1
               aseq1(ii)='-'
               l2=mod(ik2,2)
               if(j.ge.ichainterm2(ik2-1).and.j.lt.ichainterm2(ik2).and.
     &              l2.ne.0)then
                  isequence=ichar(seq2(j))
                  isequence=isequence+32
                  aseq2(ii)=char(isequence)
               else
                  aseq2(ii)=seq2(j)
               endif
               aseq3(ii)=' '
            enddo
            !iflag2=0
         endif
         iflag2=0
 7775    continue
         if(iflag3.eq.1.and.iflag1.eq.1.and.itill1.lt.mf1(i)-1)goto 7776
         if(iflag4.eq.1.and.iflag2.eq.1.and.itill2.lt.mf2(i)-1)goto 7778

         ii=ii+1
         if(mf1(i).ge.ichainterm1(ik1))then
            ik1=ik1+1
         endif
         l=mod(ik1,2)
         if(mf2(i).ge.ichainterm2(ik2))then
            ik2=ik2+1
         endif
         ll=mod(ik2,2)
         
         k=0
        if(mf1(i).ge.ichainterm1(ik1-1).and.mf1(i).lt.ichainterm1(ik1)
     &        .and.l.ne.0)then
            isequence1=ichar(seq1(mf1(i)))
            isequence1=isequence1+32
            aseq1(ii)=char(isequence1)
         else
            aseq1(ii)=seq1(mf1(i))
         endif
        if(mf2(i).ge.ichainterm2(ik2-1).and.mf2(i).lt.ichainterm2(ik2)
     &        .and.ll.ne.0)then 
            isequence=ichar(seq2(mf2(i)))
            isequence=isequence+32
            aseq2(ii)=char(isequence)
         else
            aseq2(ii)=seq2(mf2(i))
         endif
         
         dis2=sqrt((xtmf1(i)-xtmf2(i))**2+
     &        (ytmf1(i)-ytmf2(i))**2+(ztmf1(i)-ztmf2(i))**2)
         if(dis2.le.5)then
            aseq3(ii)=':'
         else
            aseq3(ii)='.'
         endif
               
         i1_old=mf1(i)+1
         i2_old=mf2(i)+1
      enddo
**********************************     
      if(ichainterm1(ik1).lt.nseq1)then
         itill1=ichainterm1(ik1)
         iflag1=1
      else
         itill1=nseq1
         iflag1=0
      endif
      if(ichainterm2(ik2).lt.nseq2)then
         itill2=ichainterm2(ik2)
         iflag2=1
      else
         itill2=nseq2
         iflag2=0
      endif

      do i=i1_old,itill1
         if(i.ge.ichainterm1(ik1))then
            ik1=ik1+1
            goto 6669
         endif
         ii=ii+1
         ll1=mod(ik1,2)
         !print *,i,ii,seq1(i),ik1
         if(i.ge.ichainterm1(ik1-1).and.i.lt.ichainterm1(ik1)
     &      .and.ll1.ne.0)then
            isequence1=ichar(seq1(i))
            isequence1=isequence1+32
            aseq1(ii)=char(isequence1)
            aseq2(ii)='-'
            aseq3(ii)=' '
            !print *,"HERE",aseq1(ii)
         else
            isequence1=ichar(seq1(i))
            aseq1(ii)=char(isequence1)
            aseq2(ii)='-'
            aseq3(ii)=' '
            !print *,aseq1(ii),"HERE",i,ichainterm1(ik1),ll1
         endif
      enddo
 6669 continue
      do i=i2_old,itill2
         if(i.ge.ichainterm2(ik2))then
            ik2=ik2+1
            goto 6668
         endif
         ii=ii+1
         ll2=mod(ik2,2)
         if(i.ge.ichainterm2(ik2-1).and.i.lt.ichainterm2(ik2)
     &      .and.ll2.ne.0)then
            aseq1(ii)='-'
            isequence=ichar(seq2(i))
            isequence=isequence+32
            aseq2(ii)=char(isequence)
            aseq3(ii)=' '
         else
            aseq1(ii)='-'
            isequence=ichar(seq2(i))
            aseq2(ii)=char(isequence)
            aseq3(ii)=' '
         endif
      enddo
 6668  continue
 6667  continue
       if(iflag1.eq.1)then
          do i=itill1,nseq1
             iflag5=0
             if(i.ge.ichainterm1(ik1))then
                ik1=ik1+1
                itill1=i
                iflag5=1
               goto 6664
            endif
            ii=ii+1
            ll1=mod(ik1,2)
                                !print *,i,ii,seq1(i),ik1
            if(i.ge.ichainterm1(ik1-1).and.i.lt.ichainterm1(ik1)
     &           .and.ll1.ne.0)then
               isequence1=ichar(seq1(i))
               isequence1=isequence1+32
               aseq1(ii)=char(isequence1)
               aseq2(ii)='-'
               aseq3(ii)=' '
            else
               isequence1=ichar(seq1(i))
               aseq1(ii)=char(isequence1)
                                !print *, "No HERe" ,aseq1(ii),i,ichainterm1(ik1-1),ll1
               aseq2(ii)='-'
               aseq3(ii)=' '
            endif
         enddo
      endif
      iflag1=0
 6664 continue
 6665 continue
      if(iflag2.eq.1)then
         do i=itill2,nseq2
            iflag6=0
            if(i.ge.ichainterm2(ik2))then
               ik2=ik2+1
               itill2=i
               iflag6=1
               goto 6666
            endif
            ii=ii+1
            ll2=mod(ik2,2)
            if(i.ge.ichainterm2(ik2-1).and.i.lt.ichainterm2(ik2)
     &           .and.ll2.ne.0)then
               aseq1(ii)='-'
               isequence=ichar(seq2(i))
               isequence=isequence+32
               aseq2(ii)=char(isequence)
               aseq3(ii)=' '
            else
               aseq1(ii)='-'
               isequence=ichar(seq2(i))
               aseq2(ii)=char(isequence)
               aseq3(ii)=' '
            endif
         enddo
      endif
      iflag2=0
 6666 continue
      if(iflag5.eq.1.and.itill1.lt.nseq1-1.and.iflag1.eq.1)goto 6667
      if(iflag6.eq.1.and.itill2.lt.nseq2-1.and.iflag2.eq.1)goto 6665

      write(*,50)
 50   format('(":" denotes the residue pairs of distance < 5.0 ',
     &     'Angstrom)')
      write(*,10)(aseq1(i),i=1,ii)
      write(*,10)(aseq3(i),i=1,ii)
      write(*,10)(aseq2(i),i=1,ii)

 10   format(10000A1)
      write(*,*) 
      if(intflag.eq.1)then
         
         write(*,2222) cov
         
 2222    format('Total Aligned Interface Coverage= ',f0.3)
         write(*,*)
         write(*,*)'(Odd no chains of both proteins are in Lower '//
     &        'case and Even no chains of both proteins are in '//
     &        'UPPER case)'
      else
         write(*,*)'(Odd no chains of both proteins are in lower '//
     &        'case and Even no chains of both proteins are in '//
     &        'UPPER case)'
      endif   
      write(*,*)

 9999 END

***********************************************************************
***********************************************************************
*     Structure superposition
***********************************************************************
***********************************************************************
***********************************************************************
      SUBROUTINE super_align
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      common/init1/invmap_ii(nmax)
      dimension gapp(100)

      TMmax=0
      TMtemp=0
      n_gapp=2
      gapp(1)=-0.6
      gapp(2)=0

c      n_gapp=11
c      do i=1,n_gapp
c         gapp(i)=-(n_gapp-i)
c      enddo

*11111111111111111111111111111111111111111111111111111111
*     get initial alignment from gapless threading
**********************************************************
      call get_initial          !gapless threading
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
      DO 111 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 222 id=1,30         !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
             
            call get_score      !calculate TM-score, score(i,j)
            !print *,"In1a",TM,TMmax
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 33
            endif
            TM_old=TM
 222     continue
 33      continue
 111  continue

*222222222222222222222222222222222222222222222222222222222
*     get initial alignment from secondary structure alignment
**********************************************************
      call get_initial2         !DP for secondary structure
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, score(i,j)
      !print *,"Ini2",TM,TMmax     
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif
      if(TM.gt.TMmax*0.2)then
*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 1111 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 2222 id=1,30	!maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)

            call get_score      !calculate TM-score, score(i,j)
      !      print *,"Ini2a",TM,TMmax
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
 2222    continue
 333     continue
 1111 continue
      endif
*******************************************************************
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
      if(TM.gt.TMmax*0.2)then
*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 88 i_gapp=2,n_gapp	 !different gap panalties
         GAP_OPEN=gapp(i_gapp)   !gap panalty
         do 888 id=1,2           !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
            
            call get_score      !calculate TM-score, score(i,j)
c            print *,"Ini5a",TM,TMmax
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 880
            endif
            TM_old=TM
 888     continue
 880     continue
 88   continue
      endif
*333333333333333333333333333333333333333333333333333333333333
*     get initial alignment from invmap0+SS
*************************************************************
      call get_initial3         !invmap0+SS
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, score(i,j)
      ! print *,"Ini3",TM,TMmax     
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif
      if(TM.gt.TMmax*0.2)then
*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 1110 i_gapp=1,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 2220 id=1,30	!maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)

            call get_score      !calculate TM-score, score(i,j)
            !print *,"Ini3a",TM,TMmax
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
               if(diff.lt.0.000001)goto 330
            endif
            TM_old=TM
 2220    continue
 330     continue
 1110 continue
      endif
*444444444444444444444444444444444444444444444444444444444
*     get initial alignment of pieces from gapless threading
**********************************************************
      call get_initial4         !gapless threading
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      do i=1,nseq2
         do j=1,nseq1
            score(i,j)=0
         enddo
      enddo

      call get_score            !TM, matrix score(i,j)
       !print *,"Ini4",TM,TMmax     
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif
      if(TM.gt.TMmax*0.2)then
*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 44 i_gapp=2,n_gapp	!different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 444 id=1,2          !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
            
            call get_score      !calculate TM-score, score(i,j)
           ! print *,"Ini4a",TM,TMmax
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
 444     continue
 44   continue
      endif




c^^^^^^^^^^^^^^^ best alignment invmap0(j) found ^^^^^^^^^^^^^^^^^^
      RETURN
      END

**************************************************************
*     get initial alignment invmap0(i) from gapless threading
**************************************************************
      subroutine get_initial
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      common/init1/invmap_ii(nmax)
      common/d0/d0,anseq
      common/d0min/d0_min
      common/speed/align_len1
      common/nchains/itercount,itercount1,icombo1,icombo2
      common/chain/ibig,ismall
      common/chains/ichainterm1(0:nmax),ichainterm2(0:nmax)
      real n_cut
      integer n_all
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms
      data w /nmax*1.0/
      dimension xtmm1(nmax),ytmm1(nmax),ztmm1(nmax)
      dimension xtmm2(nmax),ytmm2(nmax),ztmm2(nmax)
      

c      d01=d0+1.5
c      if(d01.lt.d0_min)d01=d0_min
c      d02=d01*d01
c      n_jump=5
c      n_cut=0.75
c      aL=min(nseq1,nseq2)
c      idel=aL/2.3              !minimum size of considered fragment
c      if(idel.le.5)idel=5
c      n1=-nseq2+idel
c      n2=nseq1-idel
c      GL_max=0
c      GL_maxA=0
c      GL_max_cut=0.95
c      do ishift=n1,n2,n_jump
c         L=0
c         do j=1,nseq2
c            i=j+ishift
c            if(i.ge.1.and.i.le.nseq1)then
c               L=L+1
c               invmap(j)=i
c               xtmm1(L)=xa(1,i,0)
c               ytmm1(L)=xa(1,i,0)
c               ztmm1(L)=xa(1,i,0)
c               xtmm2(L)=xa(1,j,1)
c               ytmm2(L)=xa(1,j,1)
c               ztmm2(L)=xa(1,j,1)
c            else
c               invmap(j)=-1
c            endif
c         enddo
      do ii=1,itercount
         do jj=1,itercount1
            nseq11=ichainterm1(ii)-ichainterm1(ii-1)
            nseq22=ichainterm2(jj)-ichainterm2(jj-1)
            nseq11_p=ichainterm1(ii-1)
            nseq22_p=ichainterm2(jj-1)
                                !print *,nseq11,nseq22
****  RUN INITIAL1 TO GET TM-SCORE
            d01=d0+1.5
            if(d01.lt.d0_min)d01=d0_min
            d02=d01*d01
            aL=min(nseq11,nseq22)
            idel=aL/2.5         !minimum size of considered fragment
            if(idel.le.5)idel=5
            n11=-nseq22+idel
            n22=nseq11-idel
            GL_max=0
                                !print *,nseq22_p,nseq22+nseq22_p
            do ishift=n11,n22
               L=0
               K=0
               nseq=nseq22+nseq22_p
               do j=nseq22_p,nseq
                  k=k+1
                  i=k+ishift
                  if(i.ge.1.and.i.le.nseq11)then
                     L=L+1
                     invmap(j)=i+nseq11_p
                     xtmm1(L)=xa(1,i,0)
                     ytmm1(L)=xa(1,i,0)
                     ztmm1(L)=xa(1,i,0)
                     xtmm2(L)=xa(1,j,1)
                     ytmm2(L)=xa(1,j,1)
                     ztmm2(L)=xa(1,j,1)
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
         enddo
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
      call smooth               !smooth the assignment

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

*************************
           
*************************************
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
      common/initial4/mm1(nmax),mm2(nmax)
      logical contin

      dimension ifr2(2,nmax,nmax),Lfr2(2,nmax),Lfr_max2(2),i_fr2(2)
      dimension ifr(nmax)
      dimension mm(2,nmax)
      
      fra_min=4                 !>=4,minimum fragment for search
      fra_min1=fra_min-1        !cutoff for shift, save time
      dcu0=3.85
c      dcu_min=3.65

ccc   Find the smallest continuous fragments -------->
      do i=1,nseq1
         mm(1,i)=mm1(i)
      enddo
      do i=1,nseq2
         mm(2,i)=mm2(i)
      enddo
      do k=1,2
         dcu=dcu0
         if(k.eq.1)then
            nseq0=nseq1
            r_min=nseq1/3.0     !minimum fragment, in case too small protein
         else
            nseq0=nseq2
            r_min=nseq2/3.0     !minimum fragment, in case too small protein
         endif
         if(r_min.gt.fra_min)r_min=fra_min
 20      nfr=1                  !number of fragments
         j=1                    !number of residue at nf-fragment
         ifr2(k,nfr,j)=1        !what residue
         Lfr2(k,nfr)=j          !length of the fragment
         do i=2,nseq0
            dis=diszy(k-1,i-1,i)
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
            dcu=1.1*dcu
            goto 20
         endif
      enddo
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
ccc   select what piece will be used (this may araise ansysmetry, but
ccc   only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1
ccc   if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1
      mark=1
      if(Lfr2(1,i_fr2(1)).lt.Lfr2(2,i_fr2(2)))then
         mark=1
      elseif(Lfr2(1,i_fr2(1)).gt.Lfr2(2,i_fr2(2)))then
         mark=2
      else                      !Lfr1=Lfr2
         if(nseq1.le.nseq2)then
            mark=1
         else
            mark=2
         endif
      endif
ccc   
      L_fr=Lfr2(mark,i_fr2(mark))
      do i=1,L_fr
         ifr(i)=ifr2(mark,i_fr2(mark),i)
      enddo
ccc   
      if(mark.eq.1)then         !non-redundant to get_initial1
         nseq0=nseq1
      else
         nseq0=nseq2
      endif
      if(L_fr.eq.nseq0)then
         n1=int(nseq0*0.1)      !0
         n2=int(nseq0*0.89)     !2
         j=0
         do i=n1,n2
            j=j+1
            ifr(j)=ifr(n1+j)
         enddo
         L_fr=j
      endif
      
ccc   get initial ------------->
      if(mark.eq.1)then    !nseq1 as the smallest one
         nseq1_=L_fr
         aL=min(nseq1_,nseq2)
         idel=aL/2.5            !minimum size of considered fragment
         if(idel.le.fra_min1)idel=fra_min1
         n1=-nseq2+idel         !shift1
         n2=nseq1_-idel         !shift2
         GL_max=0
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
      else                      !@@@@@@@@@@@@@@@@@@@@
         nseq2_=L_fr
         aL=min(nseq1,nseq2_)
         idel=aL/2.5            !minimum size of considered fragment
         if(idel.le.fra_min1)idel=fra_min1
         n1=-nseq2_+idel
         n2=nseq1-idel
         GL_max=0
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
      endif
      
      return
      end
**************************************************************
*    fifth initial alignement. Using local structure super   *
*    position.                                               *
**************************************************************
      subroutine get_initial5
      PARAMETER (nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/d0/d0,anseq
      common/d0min/d0_min
      common/init/invmap_i(nmax)
      common/initial4/mm1(nmax),mm2(nmax)
      common/tinit/invmap_ii(nmax)
      common/local/lflag
      common/sec/isec(nmax),jsec(nmax)
      logical contin
      integer n_all,secmatch,coilcnt
      real n_max_cut,n_all_max,n_cutt
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms
      data w /nmax*1.0/
      common/inv/invmap_a(nmax)
      common/speed/align_len1
************* points needed to tuned ******************
c     1, reduce number of iteration in get_GL from 3 to 2 or 1?
c        this should reduce the time for get_initial1 as well.
c     2, put get_DP into get_initial1 as well? Then you need to 
c        reduced the number of gapless threading by jump.
c     3, in this subroutine, n_frag, n_jump, GLmaxA_cut, d02, n_frag2
c        need to be tuned
c     4, there are serveral bugs in your codes such as the edge-condition
c        Fixing these bugs will also improve the results.
c        Also, always put the common calculations outside of the circles
************************************************************

***** setting parameters ************************************
      n_frag=5
c      n_frag=15
      n_jump=n_frag
c      n_jump=n_frag/1.5
c      n_jump=1
      n_all=0
      n_all_max=0
      n_max_cut=0.7
      n_cutt=0
      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01
      GLmaxB=0
      GLmaxA=0
      GLmaxA_cut=0.12
      n_frag1=n_frag-1
      n_frag2=n_frag/2
      m1=nseq1-n_frag+1
      m2=nseq2-n_frag+1
      do ii=1,m1,n_jump
         do jj=1,m2,n_jump
            secmatch=0
            coilcnt=0
            do k=1,n_frag
               iii=ii+k-1
               jjj=jj+k-1
               if(isec(iii).eq.jsec(jjj))then
                  secmatch=secmatch+1
               endif
               if(isec(iii).eq.1)then !.or.isec(iii).ne.4)then
                  coilcnt=coilcnt+1
               endif
            enddo
            if(secmatch.gt.n_frag2.and.coilcnt.lt.n_frag2)then
               do j=1,nseq2
                  invmap(j)=-1
               enddo
               do k=0,n_frag1
                  i=ii+k
                  j=jj+k
                  invmap(j)=i
               enddo

               call get_GL(GL)
               if(GL.gt.GLmaxA*GLmaxA_cut)then
                  !GLmaxB=GL
                  do k=1,n_frag
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
                  call u3b(w,r_1,r_2,n_frag,1,rms,u,t,ier) !u rotate r_1 to r_2
                  do i=1,nseq1
              xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
              yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
              zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
                     do j=1,nseq2
                dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
                        score(i,j)=1/(1+dd/d02)
                     enddo
                  enddo
                  
*********extract alignement with score(i,j) *****************
                  align_len1=0
                  call DP(NSEQ1,NSEQ2)
                  !print *,align_len1
                  !n_all=0
                  !do j=1,nseq2
                  !   if(invmap(j).gt.0)then
                  !      n_all=n_all+1
                  !   endif
                  !enddo
                  n_cutt=n_all_max*n_max_cut
                  !print *,n_all,n_cutt,n_all_max,n_max_cut
                  if(align_len1.gt.n_cutt)then
                     if(align_len1.gt.n_all_max)then
                        n_all_max=align_len1
                     endif
                     call get_GL(GL)
                     
                     if(GL.gt.GLmaxA)then
                        GLmaxA=GL
                        do j=1,nseq2
                           invmap_i(j)=invmap(j)
                        enddo
                     endif
                  endif
               endif
            endif
         enddo
      enddo
      
      return
      end
**************************************************************
*     smooth the secondary structure assignment
**************************************************************
      subroutine smooth
      PARAMETER(nmax=5000)
      common/sec/isec(nmax),jsec(nmax)
      common/length/nseq1,nseq2

***   smooth single -------------->
***   --x-- => -----
      do i=1,nseq1
         if(isec(i).eq.2.or.isec(i).eq.4)then
            j=isec(i)
            if(isec(i-2).ne.j)then
               if(isec(i-1).ne.j)then
                  if(isec(i+1).ne.j)then
                     if(isec(i+1).ne.j)then
                        isec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).eq.2.or.jsec(i).eq.4)then
            j=jsec(i)
            if(jsec(i-2).ne.j)then
               if(jsec(i-1).ne.j)then
                  if(jsec(i+1).ne.j)then
                     if(jsec(i+1).ne.j)then
                        jsec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo

***   smooth double -------------->
***   --xx-- => ------
      do i=1,nseq1
         if(isec(i).ne.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
         if(isec(i+3).eq.2)then
         if(isec(i+4).ne.2)then
         if(isec(i+5).ne.2)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(isec(i).ne.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
         if(isec(i+3).eq.4)then
         if(isec(i+4).ne.4)then
         if(isec(i+5).ne.4)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).ne.2)then
         if(jsec(i+1).ne.2)then
         if(jsec(i+2).eq.2)then
         if(jsec(i+3).eq.2)then
         if(jsec(i+4).ne.2)then
         if(jsec(i+5).ne.2)then
            jsec(i+2)=1
            jsec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(jsec(i).ne.4)then
         if(jsec(i+1).ne.4)then
         if(jsec(i+2).eq.4)then
         if(jsec(i+3).eq.4)then
         if(jsec(i+4).ne.4)then
         if(jsec(i+5).ne.4)then
            jsec(i+2)=1
            jsec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo

***   connect -------------->
***   x-x => xxx
      do i=1,nseq1
         if(isec(i).eq.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
            isec(i+1)=2
         endif
         endif
         endif

         if(isec(i).eq.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
            isec(i+1)=4
         endif
         endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).eq.2)then
         if(jsec(i+1).ne.2)then
         if(jsec(i+2).eq.2)then
            jsec(i+1)=2
         endif
         endif
         endif

         if(jsec(i).eq.4)then
         if(jsec(i+1).ne.4)then
         if(jsec(i+2).eq.4)then
            jsec(i+1)=4
         endif
         endif
         endif
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
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/d00/d00,d002

      dimension xo1(nmax),yo1(nmax),zo1(nmax)
      dimension xo2(nmax),yo2(nmax),zo2(nmax)
      dimension dis2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         !print *,j,i
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
*     quickly calculate TM-score with given invmap(i) in 3 iterations
****************************************************************
      subroutine get_GL11(GL)
      PARAMETER(nmax=5000)
      common/length1/nseq11,nseq22,nseq11_p,nseq22_p
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/d00/d00,d002

      dimension xo1(nmax),yo1(nmax),zo1(nmax)
      dimension xo2(nmax),yo2(nmax),zo2(nmax)
      dimension dis2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=nseq22_p,NSEQ22+nseq22_p
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
      !print *,n_al
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
      !print *,n_al
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
      common/n1n2/n1(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
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
      call TMscore8_search(d0_input,n_al,xtm1,ytm1,ztm1,n1,
     &     n_al,xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) !simplified search engine
      TM=TM*n_al/anseq          !TM-score
***   calculate score matrix score(i,j)------------------>
      do i=1,n_al
         r_2(1,i)=xtm1(i)
         r_2(2,i)=ytm1(i)
         r_2(3,i)=ztm1(i)
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
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
*     with invmap(i) calculate TM-score and martix score(i,j) for rotation 
****************************************************************
      subroutine get_score11
      PARAMETER(nmax=5000)
      common/length1/nseq11,nseq22
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ22
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
      call TMscore8_search(d0_input,n_al,xtm1,ytm1,ztm1,n1,
     &     n_al,xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) !simplified search engine
      TM=TM*n_al/nseq22          !TM-score
***   calculate score matrix score(i,j)------------------>
      do i=1,n_al
         r_2(1,i)=xtm1(i)
         r_2(2,i)=ytm1(i)
         r_2(3,i)=ztm1(i)
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      do i=1,nseq11
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq22
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
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
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


*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     n1(i)--Residue sequence number of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     n2(i)--Residue sequence number of i'th residue at the second structure
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
      subroutine TMscore8_search(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=5000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)
      dimension iL0(nmax)
      common/flag/intflag
      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:arrays 
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
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
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
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
      d_output=5                !for output alignment
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
*     find the maximum scarrays ore starting from local structures superposition
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
           call u3b(w,r_1,r_2,LL,2,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun8      !init, get scores, n_cut+i_ali(i) for iteration
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
              call score_fun8   !get scores, n_cut+i_ali(i) for iteration
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
*     n1(i)--Residue sequence number of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     n2(i)--Residue sequence number of i'th residue at the second structure
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
      subroutine TMscore8(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=5000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)
      common/flag/intflag
      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
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
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
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
      d_output=5                !for output alignment
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
        do 300 iL=1,iL_max    !on aligned residues, [1,nsene.jqA]
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
           call u3b(w,r_1,r_2,LL,2,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun8      !init, get scores, n_cut+i_ali(i) for iteration
  
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
  
              call score_fun8   !get scores, n_cut+i_ali(i) for iteration

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
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score,score_max
      common/d8/d8

      d_tmp=d
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d_tmp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         if(dis.le.d8)then
            score_sum=score_sum+1/(1+(dis/d0)**2)
         endif
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
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
*     n1(i)--Residue sequence number of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     n2(i)--Residue sequence number of i'th residue at the second structure
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
      subroutine TMscore(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=5000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)
      common/flag/intflag
      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
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
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
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
      d_output=5                !for output alignment
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
           call u3b(w,r_1,r_2,LL,2,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun       !init, get scores, n_cut+i_ali(i) for iter
              
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
      subroutine score_fun
      PARAMETER(nmax=5000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score

      d_tmp=d
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d_tmp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         score_sum=score_sum+1/(1+(dis/d0)**2)
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
         goto 21
      endif
      score=score_sum/float(nseqB) !TM-score

      return
      end

**************************************************************************
*     Dynamic programming for alignment.                                 *
*     Input: score(i,j), and gap_open                                    *
*     Output: invmap(j)                                                  *
*                                                                        *
*     Please note this subroutine is not a correct implementation of     *
*     the N-W dynamic programming because the score tracks back only     *
*     one layer of the matrix. This code was exploited in TM-align       *
*     because it is about 1.5 times faster than a complete N-W code      *
*     and does not influence much the final structure alignment result.  *
*     To prevent cross-alignment the cross-aligned regions are ignored   *
*     both while filling up the scoring matrix and also during traceback.*
*     When the interface alignment option is used the DP is carried out  *
*     with a higher weight given for aligning interface residues.        *
**************************************************************************
      SUBROUTINE DP(NSEQ1,NSEQ2)
      PARAMETER(nmax=5000)
      PARAMETER(nmax2=10000) 
      LOGICAL*1 DIR
      LOGICAL*1 DIR1
      LOGICAL*1 DIR2
      LOGICAL*1 mtc 
      common/combi/npm
      common/arr/ ib(1000000,50),ic(1000000,50),id(1000000,50)  
      common/arr1/ia(1000000,50)
      common/nchains/itercount,itercount1,icombo1,icombo2
      common/chain/ibig,ismall
      common/chains/ichainterm1(0:nmax),ichainterm2(0:nmax)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/initial5/inter1(nmax),inter2(nmax),inter3(nmax)
      common/initial5/inter4(nmax)
      common/initial6/ifcounter1,ifcounter2
      common/flag/intflag,domflag
      common/choice/ichachoice,ichain,ii
      common/lens/len_first1,len_first2,len_second1,len_second2
      dimension DIR(0:nmax,0:nmax),VAL(0:nmax,0:nmax)
      dimension DIR1(0:nmax,0:nmax),VAL1(0:nmax,0:nmax)
      dimension DIR2(0:nmax,0:nmax),VAL2(0:nmax,0:nmax)
      integer fact1,fact2,fact3
      REAL H,V
      common/speed/align_len1
      Dflag=0
      len_second1=0
      len_second2=0
      len1=0
      len2=0
      align_len1=0
******  IF Domain Swapping or Gene Fusion Cross Align Allowed ******

      IF(domflag.eq.1)then

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
          align_len1=align_len1+1
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



******   If it is a dimer and interface alignment is required ******

      ELSEIF(intflag.EQ.1)then

***   define matrix for both chains of both proteins
      len_second1=nseq1-len_first1
      len_second2=nseq2-len_first2
     
***   initialize the first matrix:
      val1(0,0)=0
      do i=1,len_first1
        dir1(i,0)=.false.
        val1(i,0)=0
      enddo
      do j=1,len_first2
        dir1(0,j)=.false.
        val1(0,j)=0
        invmap(j)=-1
      enddo
***   decide first matrix and path:
      DO j=1,len_first2
         DO i=1,len_first1
            Dflag=0
            if(intflag.eq.1)then
               if(inter1(i).eq.1 .and. inter3(j).eq. 1)then
                  Dflag=2
               endif
               
               if(Dflag.eq.2)then
                  D=VAL1(i-1,j-1)+(100*SCORE(i,j))  !+0.8
                  H=VAL1(i-1,j)
                  GP=5*GAP_OPEN
                  if(DIR1(i-1,j))H=H+GP
                  V=VAL1(i,j-1)
                  if(DIR1(i,j-1))V=V+GP
               else         
                  D=VAL1(i-1,j-1)+(SCORE(i,j))
                  H=VAL1(i-1,j)
                  if(DIR1(i-1,j))H=H+GAP_OPEN
                  V=VAL1(i,j-1)
                  if(DIR1(i,j-1))V=V+GAP_OPEN
               endif
            else
               D=VAL1(i-1,j-1)+SCORE(i,j)
               H=VAL1(i-1,j)
               if(DIR1(i-1,j))H=H+GAP_OPEN
               V=VAL1(i,j-1)
               if(DIR1(i,j-1))V=V+GAP_OPEN
            endif   
         
            IF((D.GE.H).AND.(D.GE.V)) THEN
               DIR1(I,J)=.true.
               VAL1(i,j)=D
            ELSE
               DIR1(I,J)=.false.
               if(V.GE.H)then
                  val1(i,j)=v
               else
                  val1(i,j)=h
               end if
            ENDIF
         ENDDO
      ENDDO
***   Initialize second matrix 
      len1=len_first1+1
      len2=len_second2+1
      val2(len1,len2)=val1(len_first1,len_first2)      
     
      do i=len1,nseq1
         dir2(i,len_first2)=.false.
         val2(i,len_first2)=val1(len_first1,len_first2)
      enddo
      do j=len2,nseq2
         dir2(len_first1,j)=.false.
         val2(len_first1,j)=val1(len_first1,len_first2)
         invmap(j)=-1
      enddo
      if(ifcounter1.ge.ifcounter2)then
         ifcount=ifcounter1
      else
         ifcount=ifcounter2
      endif
***   Decide second matrix and path
        DO j=len2,nseq2
         DO i=len1,nseq1
            Dflag=0
            if(intflag.eq.1)then

               if(inter2(i).eq.1 .and. inter4(j).eq. 1)then
                  Dflag=2
               endif
             
               if(Dflag.eq.2)then
                  D=VAL2(i-1,j-1)+(100*SCORE(i,j))
                  H=VAL2(i-1,j)
                  GP=5*GAP_OPEN
                  if(DIR2(i-1,j))H=H+GP
                  V=VAL2(i,j-1)
                  if(DIR2(i,j-1))V=V+GP
               else
                  D=VAL2(i-1,j-1)+(SCORE(i,j))
                  H=VAL2(i-1,j)
                  if(DIR2(i-1,j))H=H+GAP_OPEN
                  V=VAL2(i,j-1)
                  if(DIR2(i,j-1))V=V+GAP_OPEN
               endif
            else
               D=VAL2(i-1,j-1)+SCORE(i,j)
               H=VAL2(i-1,j)
               if(DIR2(i-1,j))H=H+GAP_OPEN
               V=VAL2(i,j-1)
               if(DIR2(i,j-1))V=V+GAP_OPEN
            endif   
            IF((D.GE.H).AND.(D.GE.V)) THEN
               DIR2(I,J)=.true.
               VAL2(i,j)=D
            ELSE
               DIR2(I,J)=.false.
               if(V.GE.H)then
                  val2(i,j)=v
               else
                  val2(i,j)=h
               endif
            ENDIF
         ENDDO
      ENDDO
***   extract the alignment:
      i=NSEQ1
      j=NSEQ2
      len1=len_first1+1
      len2=len_first2+1
 
      DO WHILE((i.GT.0).AND.(j.GT.0))
         IF(i.eq.len_first1.AND.j.eq.len_first2)THEN
            i=i-1
            j=j-1
         ELSE
            if(i.gt.len_first1.AND.j.eq.len_first2)then
               i=i-1
            else
               if(i.eq.len_first1.AND.j.gt.len_first2)then
                  j=j-1
               else
                  if(i.gt.len_first1.AND.j.gt.len_first2)then
                     IF(DIR2(i,j))THEN
                        align_len1=align_len1+1
                        invmap(j)=i
                        i=i-1
                        j=j-1
                     ELSE
                        H=VAL2(i-1,j)
                        if(DIR2(i-1,j))H=H+GAP_OPEN
                        V=VAL2(i,j-1)
                        if(DIR2(i,j-1))V=V+GAP_OPEN
                        IF(V.GE.H) THEN
                           j=j-1
                        ELSE
                           i=i-1
                        ENDIF
                     ENDIF
                  else
                     IF(DIR1(i,j))THEN
                        invmap(j)=i
                        i=i-1
                        j=j-1
                     ELSE
                        H=VAL1(i-1,j)
                        if(DIR1(i-1,j))H=H+GAP_OPEN
                        V=VAL1(i,j-1)
                        if(DIR1(i,j-1))V=V+GAP_OPEN
                        IF(V.GE.H) THEN
                           j=j-1
                        ELSE
                           i=i-1
                        ENDIF
                     ENDIF
                  endif
               endif
            endif
         ENDIF
      ENDDO
 
******** IF IT IS A MULTIMER AND NO DOM SWAP OR NO INTERFACE ALIGNMENT *****

      ELSE

    
***   initialize the first matrix:
      val1(0,0)=0
      do i=1,nseq1
         do j=1,nseq2
            dir1(i,j)=.false.
            val1(i,j)=0
         enddo
      enddo
      do j=1,nseq2
         invmap(j)=-1
      enddo
      ichainterm2(0)=1
      ichainterm1(0)=1
      ichain=0
      j=0
      i=0
      do jj=1,ismall
         ichain=ia(icombo1,jj)
         ichachoice=id(icombo2,jj)
         !print *, ichain,ichachoice,ichainterm2(ichachoice)
         DO j=ichainterm2(ichachoice-1)-1,ichainterm2(ichachoice)
            DO i=ichainterm1(ichain-1)-1,ichainterm1(ichain)
               D=VAL1(i-1,j-1)+SCORE(i,j)
               H=VAL1(i-1,j)
               if(DIR1(i-1,j))H=H+GAP_OPEN
               V=VAL1(i,j-1)
               if(DIR1(i,j-1))V=V+GAP_OPEN
               IF((D.GE.H).AND.(D.GE.V)) THEN
                  DIR1(I,J)=.true.
                  VAL1(i,j)=D
               ELSE
                  DIR1(I,J)=.false.
                  if(V.GE.H)then
                     val1(i,j)=v
                  else
                     val1(i,j)=h
                  end if
               ENDIF
            ENDDO
 1234       continue
         ENDDO
      enddo
    

***   extract the alignment:
      i=ichainterm1(ichain)
      j=ichainterm2(ichachoice)
      ibegin1=0
      ibegin2=0
      !print *,i,j,ismall
      do jj=ismall,0,-1
         ibegin2=id(icombo2,jj)
         ibegin1=ia(icombo1,jj)
	 !print *,i,j,ichainterm1(ibegin1),ichainterm2(ibegin2) 
         IF((i.eq.ichainterm1(ibegin1)).AND.(j.eq.ichainterm2
     &        (ibegin2)))THEN
            i=i-1
            j=j-1
	    !print *,i,j,ichainterm1(ibegin1-1)-1,ichainterm2(ibegin2-1)-1   	
            !invmap(j)=i
         endif
         DO WHILE((i.GT.ichainterm1(ibegin1-1)-1).AND.(j.GT.ichainterm2
     &          (ibegin2-1)-1))
            IF((i.eq.ichainterm1(ibegin1)).AND.(j.eq.ichainterm2
     &              (ibegin2)))then
                  i=i-1
                  j=j-1
                  !invmap(j)=i
               else
               if((i.gt.ichainterm1(ibegin1-1)).AND.(j.eq.ichainterm2
     &                 (ibegin2-1)))then
                  i=i-1
               else
                if((i.eq.ichainterm1(ibegin1-1)).AND.(j.gt.ichainterm2
     &                 (ibegin2-1)))then
                     j=j-1
                  else
                     IF(DIR1(i,j))THEN
                        !print *,i,j
                        invmap(j)=i
                        align_len1=align_len1+1
                        i=i-1
                        j=j-1
                     ELSE
                        H=VAL1(i-1,j)
                        if(DIR1(i-1,j))H=H+GAP_OPEN
                        V=VAL1(i,j-1)
                        if(DIR1(i,j-1))V=V+GAP_OPEN
                        IF(V.GE.H) THEN
                           j=j-1
                        ELSE
                           i=i-1
                        ENDIF
                     ENDIF
                  endif
               endif
            ENDIF
         ENDDO
      enddo
      
      ENDIF

c^^^^^^^^^^^^^^^Dynamical programming done ^^^^^^^^^^^^^^^^^^^
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Subtourine to determine interface residues based on distance   c
c       distance cutoff of 3.4.If the distnce between the CA atoms of  c
c       any residue in different chains of the same protein are within c
c       the cutoff distance then the pair is considered part of the    c
c       interface region.                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine interface(pdb1,pdb2)
      PARAMETER (nmax=5000)
      COMMON/BACKBONE/XXA(3,nmax,0:1)   
      COMMON/BACKBONE/XXB(3,nmax,0:1)
      common/length/inseq1,inseq2
      common/initial4/imm1(nmax),imm2(nmax),imm3(nmax),imm4(nmax)
      common/initial5/inter1(nmax),inter2(nmax),inter3(nmax)
      common/initial5/inter4(nmax)
      common/initial6/ifcounter1,ifcounter2
      common/lens/len_first1,len_first2,len_second1,len_second2
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      character*100 pdb1,pdb2
      character*100 is,idu
      character*3 iaa(-1:20),iaanam,iss1(nmax),iss2(nmax),iss3(nmax)
      character*3 iss4(nmax)
      character iseq1(0:nmax),iseq2(0:nmax),iseq3(0:nmax),iseq4(0:nmax)
      dimension ixx(nmax),iyy(nmax),zz(nmax)
      dimension im1(nmax),im2(nmax),im3(nmax),im4(nmax)
      common/cutoff/int_cut
      common/initial4/mempo(nmax),mempo1(nmax)
      integer ii,jj,chk
      real*8 dist,xaa,yaa,zaa
ccc

      data iaa/ 'BCK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP','CYX'/
      character*1 islc(-1:20)
      data islc/'X','G','A','S','C',
     &     'V','T','I','P','M',
     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W','C'/

cccccc read data from first CA file :
       open(unit=10,file=pdb1,status='old')
       i=0
       iterflag=0
       ilen_first1=0
       ilen_first2=0
       !if(int_cut.eq."")then
       !   int_cut=8
       !endif
       do while (.true.)
         read(10,7001,end=7010) is
         if(i.gt.0.and.is(1:3).eq.'TER'.and.iterflag.eq.1)goto 7010
         if(i.gt.0.and.is(1:3).eq.'TER'.and.iterflag.eq.0) then
               iterflag=1      
               goto 7016
         endif
         if(is(1:3).eq.'ATO'.and.iterflag.eq.0)then
                goto 700
         else   
                goto 7015
         endif 
 700        if(is(13:16).eq.'CA  '.or.is(13:16).eq.' CA '.or
     &           .is(13:16).eq.'  CA')then
               if(is(17:17).eq.' '.or.is(17:17).eq.'A')then
                  i=i+1
                  ilen_first1=i
                  read(is,7000)idu,iaanam,idu,imm1(i),idu,
     $                 xxa(1,i,0),xxa(2,i,0),xxa(3,i,0)
                  do j=-1,20
                     if(iaanam.eq.iaa(j))then
                        iseq1(i)=islc(j)
                        goto 71
                     endif
                  enddo
                  iseq1(i)=islc(-1)
 71               continue
                  iss1(i)=iaanam
                  if(i.ge.nmax)goto 7010
               endif
            endif
         goto 7016

 7015       if(is(13:16).eq.'CA  '.or.is(13:16).eq.' CA '.or
     &           .is(13:16).eq.'  CA')then
               if(is(17:17).eq.' '.or.is(17:17).eq.'A')then
                  i=i+1
                  ilen_first2=i
                  iaanam=is(17:21)
                  read(is,7000)idu,iaanam,idu,imm2(i),idu,
     $                 xxa(1,i,1),xxa(2,i,1),xxa(3,i,1)
                  do j=-1,20
                     if(iaanam.eq.iaa(j))then
                        iseq1(i)=islc(j)
                        goto 72
                     endif
                  enddo
                  iseq2(i)=islc(-1)
 72               continue
                  iss2(i)=iaanam
                  if(i.ge.nmax)goto 7010
               endif
           endif
 7016 continue     
      enddo
 7010 continue
 7000 format(A17,A3,A3,i4,A4,3F8.3)
 7001 format(A100)
      close(10)
      inseq=ilen_first2 
      ilen_first2=ilen_first2-ilen_first1
      len_first1=ilen_first1
      ii=0
      jj=0
      kk=0
      ll=0
      distance=0
      xaa=0
      yaa=0
      zaa=0
      dist=0 
      temp=0
      icount=0
      int1=0
      int2=0
****** Initializing interface array ********      
      do i=0,inseq
         inter1(i)=0
      enddo
      do j=0,inseq
         inter2(j)=0
      enddo

***** Storing interface residues in array with bit value
      do ii=1,ilen_first1
         do jj=ilen_first2,inseq
               xaa=(xxa(1,ii,0)-xxa(1,jj,1))**2
               yaa=(xxa(2,ii,0)-xxa(2,jj,1))**2
               zaa=(xxa(3,ii,0)-xxa(3,jj,1))**2
            distance=(xaa+yaa+zaa)
            if(distance.lt.0)then
              distance=0-distance
            endif  
            dist=sqrt(distance) 
            if(dist.le.int_cut)then
                if(temp.eq.imm1(ii))goto 7017
                do chk=1,inseq
                   if (imm2(jj).eq.mempo(chk))goto 7018
                enddo   
                ifcounter1=ifcounter1+1
                int1=imm1(ii)
                int2=ilen_first1+imm2(jj)
                inter1(int1)=1  
                inter2(int2)=1  
                temp=imm1(ii)
                mempo(ii)=imm2(jj)
            endif
 7018       continue
         enddo
 7017    continue
      enddo
cccccc read second CA file  ccccccccccccccccccccccccccc     
      open(unit=10,file=pdb2,status='old')
       i=0
       iterflag=0
       ilen_first1=0
       ilen_first2=0
       do while (.true.)
         read(10,8003,end=8010) is
         if(i.gt.0.and.is(1:3).eq.'TER'.and.iterflag.eq.1)goto 8010
         if(i.gt.0.and.is(1:3).eq.'TER'.and.iterflag.eq.0) then
               iterflag=1      
               goto 8016
         endif
         if(is(1:3).eq.'ATO'.and.iterflag.eq.0)then
                goto 800
         else   
                goto 8015
         endif 
 800        if(is(13:16).eq.'CA  '.or.is(13:16).eq.' CA '.or
     &           .is(13:16).eq.'  CA')then
               if(is(17:17).eq.' '.or.is(17:17).eq.'A')then
                  i=i+1
                  ilen_first1=i
                  read(is,8002)idu,iaanam,idu,imm3(i),idu,
     $                 xxb(1,i,0),xxb(2,i,0),xxb(3,i,0)
                  do j=-1,20
                     if(iaanam.eq.iaa(j))then
                        iseq3(i)=islc(j)
                        goto 81
                     endif
                  enddo
                  iseq3(i)=islc(-1)
 81               continue
                  iss3(i)=iaanam
                  if(i.ge.nmax)goto 8010
               endif
            endif
         goto 8016

 8015       if(is(13:16).eq.'CA  '.or.is(13:16).eq.' CA '.or
     &           .is(13:16).eq.'  CA')then
               if(is(17:17).eq.' '.or.is(17:17).eq.'A')then
                  i=i+1
                  ilen_first2=i
                  read(is,8002)idu,iaanam,idu,imm4(i),idu,
     $                 xxb(1,i,1),xxb(2,i,1),xxb(3,i,1)
                  do j=-1,20
                     if(iaanam.eq.iaa(j))then
                        iseq4(i)=islc(j)
                        goto 82
                     endif
                  enddo
                  iseq4(i)=islc(-1)
 82               continue
                  iss4(i)=iaanam
                  if(i.ge.nmax)goto 8010
               endif
           endif
 8016 continue     
      enddo
 8010 continue
 8002 format(A17,A3,A2,i4,A4,3F8.3)
 8003 format(A100)
      close(10)
      inseq=ilen_first2 
      ilen_first2=ilen_first2-ilen_first1
      len_first2=ilen_first1
      ii=0
      jj=0
      kk=0
      ll=0
      distance=0
      xaa=0
      yaa=0
      zaa=0
      dist=0 
      temp=0
      icount=0
      int1=0
      int2=0
****** Initializing interface array ********      
      do i=0,ilen_first1
         inter3(i)=0
      enddo
      do j=0,ilen_first2
         inter4(j)=0
      enddo

***** Storing interface residues in array with bit value
      do ii=1,ilen_first1
         do jj=ilen_first2,inseq
               xaa=(xxb(1,ii,0)-xxb(1,jj,1))**2
               yaa=(xxb(2,ii,0)-xxb(2,jj,1))**2
               zaa=(xxb(3,ii,0)-xxb(3,jj,1))**2
            distance=(xaa+yaa+zaa)
            if(distance.lt.0)then
              distance=0-distance
            endif  
            dist=sqrt(distance)
            if(dist.le.int_cut)then
               if(temp.eq.imm3(ii))goto 8017
                do chk=1,inseq
                   if (imm4(jj).eq.mempo1(chk))goto 8018
                enddo   
                ifcounter2=ifcounter2+1
                int1=imm3(ii)
                int2=ilen_first1+imm4(jj)
                inter3(int1)=1  
                inter4(int2)=1 
                temp=imm3(ii)
                mempo1(ii)=imm4(jj)
            endif
 8018       continue
         enddo
 8017 continue
      enddo

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
      
c------------------------------------------------------------------
c     This program generates nPk.Please change n and k accordingly -
c------------------------------------------------------------------  
      subroutine nexper(ibig,a,mtc,even)
      PARAMETER (nmax2=10000)
c     next permutation of {1,...,n}. Ref NW p 59.
      integer a(nmax2),s,d
      logical mtc,even

      n=ibig
      if(mtc)goto 10
      nm3=n-3
      do 1 i=1,n
    1 a(i)=i
      mtc=.true.
    5 even=.true.
      if(n.eq.1)goto 8
    6 if(a(n).ne.1.or.a(1).ne.2+mod(n,2))return
      if(n.le.3)goto 8
      do 7 i=1,nm3
      if(a(i+1).ne.a(i)+1)return
    7 continue
    8 mtc=.false.
      return
   10 if(n.eq.1)goto 27
      if(.not.even)goto 20
      ia=a(1)
      a(1)=a(2)
      a(2)=ia
      even=.false.
      goto 6
   20 s=0
      do 26 i1=2,n
   25 ia=a(i1)
      i=i1-1
      d=0
      do 30 j=1,i
   30 if(a(j).gt.ia) d=d+1
      s=d+s
      if(d.ne.i*mod(s,2)) goto 35
   26 continue
   27 a(1)=0
      goto 8
   35 m=mod(s+1,2)*(n+1)
      do 40 j=1,i
      if(isign(1,a(j)-ia).eq.isign(1,a(j)-m))goto 40
      m=a(j)
      l=j
   40 continue
      a(l)=ia
      a(i1)=m
      even=.true.
      return
      end

