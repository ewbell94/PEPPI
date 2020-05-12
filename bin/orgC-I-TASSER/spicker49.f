ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This program, Structure-PCIKER (SPICKER), is aimed at selecting the 
c     best fold by clustering trajectors produced by CABS model. 
c     Comments and bug reports should be addressed to zhang6@buffalo.edu.
c
c     Input files includes:
c       'rmsinp'---Mandatory, length of protein & piece for RMSD calculation;
c       'seq.dat'--Mandatory, sequence file, for output of PDB models.
c       'tra.in'---Mandatory, list of trajectory names used for clustering.
c       files in 'tra.in'---Mandatory, trajectories of structures.
c       'CA'-------Optional, native structure, for comparison to native.
c       'TEMP'-----Optional, template file, for structure close to template.
c
c     Output files includes:
c       'rst.dat'-----summary of clustering results;
c       'str.txt'-----list of structure in cluster;
c       'combo*.pdb'--PDB format of cluster centroids;
c       'closc*.pdb'--PDB format of structures closest to centroids;
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        1         2         3         4         5         6         7 !
c 345678901234567890123456789012345678901234567890123456789012345678901234567
      program cluster
      parameter(ndim=2000)      !Length
      parameter(nst=20200)      !number of used structure, maximum allowed
c     parameter(nst=100)      !number of used structure, maximum allowed
      parameter(ntr=20000)      !number of trajectory files
      parameter(ncl=10)         !number of clusters
      parameter(nclc=50)        !number of closc1_i.pdb
      character filen(ntr)*72,c2*6,c3*20,c4*20,seq(ndim)*3,protein*10
      character txt1*70,txt2*70,txt3*70
ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   
      dimension xtemp(ndim),ytemp(ndim),ztemp(ndim) !structure close to templ
      dimension xt(ndim),yt(ndim),zt(ndim) !temporal coordinate
      dimension x_n(ndim),y_n(ndim),z_n(ndim) !native structure
      dimension x(ndim,nst),y(ndim,nst),z(ndim,nst) !used structures
      dimension amat(nst,nst)    !RMSD matrics
      dimension mark(nst)       !for removing used structures
      dimension n_str_near(nst) !numbef of neighboring structures
      dimension itra(nst),istr(nst),E(nst)
      dimension n_str_cl(ncl),n_str_cl_ex(ncl)
      dimension i_cl(ncl),rmsd_cl_cut(ncl)
      dimension E_combo(ncl)
      dimension i_str_cl(ncl,nst),i_str_cl_ex(ncl,nst)
      dimension rmsd_str(nst)
      dimension xs(ndim),ys(ndim),zs(ndim)
      dimension xc(ncl,ndim),yc(ncl,ndim),zc(ncl,ndim)
      dimension xcl(ncl,ndim),ycl(ncl,ndim),zcl(ncl,ndim)
      dimension xc3(ncl,ndim),yc3(ncl,ndim),zc3(ncl,ndim)
      dimension rmsd_close_min(ncl) !RMSD between 'combo.pdb' and 'closm.pdb'
      dimension i1(100),i2(100)
      integer q(ndim)
      dimension R_in(100),R_ex(100),Rc_in(100),Rc_ex(100)
      dimension iden1(100),iden2(100)
      dimension order(nst)
      dimension ires(ndim)
      double precision rmsd_a,rmsd2_a
      
      dimension x1(ndim),y1(ndim),z1(ndim),nn1(ndim)
      dimension x2(ndim),y2(ndim),z2(ndim),nn2(ndim)
      dimension r_decoy_combo(nst),key(nst)
      dimension m_closc(nclc),rmsd_cc(ncl,nclc)
      dimension x_cc(ncl,nclc,ndim)
      dimension y_cc(ncl,nclc,ndim)
      dimension z_cc(ncl,nclc,ndim)
      
**********************************************************************
****  input:
      open(1,file='rmsinp',status='old') !read Lch
      open(2,file='seq.dat',status='old') !read SEQ for model.pdb
      open(3,file='CA',status='unknown') !for comparison to native
      open(4,file='tra.in',status='old') !information for trajectories
      
****  output:
      open(20,file='rst.dat',status='unknown')
      open(21,file='str.txt',status='unknown') !structures in cluster
      open(22,file='RMSD.list',status='unknown') !structures in cluster
**********************************************************************

      read(1,*)n1,n2            !calculate RMSD between [n1,n2]
      read(1,*)Lch
      read(1,*)protein
      close(1)

******** following parameter has been optimized ********************
      nc_max=10
      RMSD_cut_initial=8
      RMSD_cut_min=3.5
      RMSD_cut_max=12
      ratio1=0.7
      ratio2=0.15
***   n_para>0: old scheme good for TASSER (2 RMSD_cut, 2 ratio_cuts)
***   n_para=-1, R_cut based on <rmsd> and delta_rmsd, ra/rb=0.5/0.8
***   n_para=-2, R_cut based on <rmsd> and delta_rmsd, ra/rb=0.3/0.15
***   n_para=-3, R_cut based on ratio0=15% and RMSD_cut0=5.25A
***   n_closc>0: closc from all decoys; n_closc<0, closc from clustered
      read(4,*)n_tra,n_para,n_closc
********************************************************************

***   find out the total number of structures & structure closest to template:
***   The total number of structure will be used for picking clustering struct.
      R_temp_min=10000
      i_str_all=0
      do i_tra=1,n_tra
         i_str_tra=0
         read(4,'(A72)')filen(i_tra)
         open(10,file=filen(i_tra),status='old')
 10      read(10,*,end=19,err=19)Lch1,energ
         do i=1,Lch1
            read(10,*,end=19,err=19)xt(i),yt(i),zt(i)
         enddo
         i_str_tra=i_str_tra+1
         i_str_all=i_str_all+1
         goto 10
 19      close(10)
      enddo
      n_str_all=i_str_all
c     write(*,*)'Total number structures=',n_str_all
      close(4)
c^^^^^^^^^^^^^^^^^^^^^^^^ n_str_all done ^^^^^^^^^^^^^^^^^^^^^^^^
      
cccccccccccc read native structure cccccccccccccccccccc
      i=0
 11   read(3,1239,end=5)tx,ii,tx,tx,ii,tx,a1,a2,a3
      if(ii.ge.1.and.ii.le.Lch)then
         i=i+1
         ires(i)=ii
         x_n(i)=a1
         y_n(i)=a2
         z_n(i)=a3
         r_1(1,i)=a1
         r_1(2,i)=a2
         r_1(3,i)=a3
      endif
      goto 11
 5    continue
      Lch_n=i                   !length of native
 1239 format(A4,I7,A4,A7,I4,A4,3F8.3)
c^^^^^^^^^^^^^^^ read native done ^^^^^^^^^^^^^^^^^^^^^^
      
ccc read structures from trajectories and find best structure in 20001 cccccccc
      i_str_tra_a=0             !<i_str_tra>
      rmsd_min_all=100
      rmsd_min_use=100
      n_max=nst                 !maximum structures to handle
      delta=float(n_str_all)/float(n_max) !take a structure in each delta
      if(delta.lt.1)delta=1
      i_str=1			!number of structure used for clustering
      i_str_all=0
      E_min=100000
      E_min_all=100000
      write(22,*)'RMSD   Energy   #str   trajectory'
      do 100 i_tra=1,n_tra
         open(10,file=filen(i_tra),status='old')
         i_str_tra=0            !order number of the structure in this file
 20      read(10,*,end=29,err=29)Lch1,energ
         do i=1,Lch1
            read(10,*,end=29,err=29)xt(i),yt(i),zt(i)
         enddo
         i_str_tra=i_str_tra+1
         i_str_all=i_str_all+1
***   calculate RMSD of structure to native----------->
         if(Lch_n.gt.1)then
            do i=1,Lch_n
               r_2(1,i)=xt(ires(i))
               r_2(2,i)=yt(ires(i))
               r_2(3,i)=zt(ires(i))
            enddo
            call u3b(w,r_1,r_2,Lch_n,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch_n)
            if(rmsd_min_all.gt.armsd)then
               rmsd_min_all=armsd
               E_rma=energ
               itra_rmsd_min_all=i_tra
               istr_rmsd_min_all=i_str_tra
            endif
            if(E_min_all.gt.energ)then
               rmsd_E_min_all=armsd
               E_min_all=energ
               itra_E_min_all=i_tra
               istr_E_min_all=i_str_tra
            endif
         endif
***   select structure for clustering-------->
         if(i_str_all.ge.i_str*delta)then
            if(i_str.gt.n_max) goto 29
            order(i_str_tra)=order(i_str_tra)+1
            i_str_tra_a=i_str_tra_a+i_str_tra
            itra(i_str)=i_tra
            istr(i_str)=i_str_tra
            E(i_str)=energ
            n_str=i_str
            do i=1,Lch
               x(i,i_str)=xt(i)
               y(i,i_str)=yt(i)
               z(i,i_str)=zt(i)
            enddo
            if(Lch_n.gt.1)then
               if(rmsd_min_use.gt.armsd)then
                  rmsd_min_use=armsd
                  E_rmu=energ
                  itra_rmsd_min_use=i_tra
                  istr_rmsd_min_use=i_str_tra
               endif
               if(E_min.gt.energ)then
                  E_min=energ
                  rmsd_E_min=armsd
                  itra_E_min=i_tra
                  istr_E_min=i_str_tra
               endif
               rmsd_str(i_str)=armsd
               write(22,61)armsd,energ,i_str_tra,filen(i_tra)
 61            format(f8.3,1x,f11.2,1x,i10,1x,a18)
            else
               rmsd_str(i_str)=-1
            endif
            i_str=i_str+1
         endif
***   
         goto 20
 29      close(10)
 100  continue
      close(22)
c^^^^^^^^^^^^^^^^^read structures finished ^^^^^^^^^^^^^^^^^^^^
c     write(*,*)i_str,n_str,itra(n_str),istr(n_str),n_max
      
ccc   output distribution of i_str_tra ------->
      i_str_tra_a=float(i_str_tra_a)/float(n_str)
      write(21,*)'delta=',delta
      write(21,*)'<i_str_tra>=',i_str_tra_a
      write(21,*)'distribution------->'
      do i=1,nst
         if(order(i).gt.0.5)then
            write(21,*)i,order(i)
         endif
      enddo
      
cccccccccccccc calculate RMSD matrics ccccccccccccccccccccccccccc
c     write(*,*)'number of used structures=',n_str
      rmsd_a=0
      rmsd2_a=0
      n_rmsd=0
      do i=1,n_str
         mark(i)=1
         do j=i,n_str
            do k=1,Lch
               r_1(1,k)=x(k,i)
               r_1(2,k)=y(k,i)
               r_1(3,k)=z(k,i)
               r_2(1,k)=x(k,j)
               r_2(2,k)=y(k,j)
               r_2(3,k)=z(k,j)
            enddo
            call u3b(w,r_1,r_2,Lch,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch) !RMSD12
            amat(i,j)=armsd
            amat(j,i)=armsd

            n_rmsd=n_rmsd+1
            rmsd_a=rmsd_a+armsd
            rmsd2_a=rmsd2_a+armsd*armsd
         enddo
      enddo
      rmsd_a=rmsd_a/n_rmsd
      rmsd2_a=rmsd2_a/n_rmsd
      rmsd_delta=sqrt(rmsd2_a-rmsd_a**2) ! 68.2% is in [-d,+d]
c^^^^^^^^^^^^RMSD matrics finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
******** reset RMSD_cutoff parameters:
      if(n_para.eq.-1)then      !parameters for quark
         ra=0.5
         rb=0.8
         RMSD_cut_initial=rmsd_a*ra
         if(RMSD_cut_initial.lt.2)RMSD_cut_initial=2
         RMSD_cut_min=RMSD_cut_initial-rb*rmsd_delta
         RMSD_cut_max=RMSD_cut_initial+rb*rmsd_delta
         if(RMSD_cut_min.lt.1)RMSD_cut_min=1
      elseif(n_para.eq.-2)then  !parameters for rosetta n_tra>3000
         ra=0.3
         rb=0.15
         RMSD_cut_initial=rmsd_a*ra
         if(RMSD_cut_initial.lt.2)RMSD_cut_initial=2
         RMSD_cut_min=RMSD_cut_initial-rb*rmsd_delta
         RMSD_cut_max=RMSD_cut_initial+rb*rmsd_delta
         if(RMSD_cut_min.lt.1)RMSD_cut_min=1
      endif
      write(*,*)RMSD_cut_initial,RMSD_cut_min,RMSD_cut_max
      
**************************************************************
*     find out the biggest cluster and decide the RMSD_cut_min
**************************************************************
      if(n_para.eq.-3)then
         ratio0=0.15
         RMSD_cut0=5.25
         RMSD_cut=1
 141     do j=1,n_str
            n_str_near(j)=0     !number of structure in cluster
            do k=1,n_str
               if(amat(j,k).lt.RMSD_cut)then
                  n_str_near(j)=n_str_near(j)+1
               endif
            enddo
         enddo
***   find out the biggest cluster ----------------->
         n_str_cl_max=0         !maximum number of structures in cluster
         do j=1,n_str
            if(n_str_near(j).gt.n_str_cl_max)then
               n_str_cl_max=n_str_near(j)
            endif
         enddo
***   check the size of the cluster ------------->
         ratio=float(n_str_cl_max)/float(n_str)
         if(ratio.lt.ratio0.and.RMSD_cut.lt.RMSD_cut0)then
            RMSD_cut=RMSD_cut+0.1
            goto 141
         endif
      else
         RMSD_cut=RMSD_cut_initial
 140     do j=1,n_str
            n_str_near(j)=0     !number of structure in cluster
            do k=1,n_str
               if(amat(j,k).lt.RMSD_cut)then
                  n_str_near(j)=n_str_near(j)+1
               endif
            enddo
         enddo
***   find out the biggest cluster ----------------->
         n_str_cl_max=0         !maximum number of structures in cluster
         do j=1,n_str
            if(n_str_near(j).gt.n_str_cl_max)then
               n_str_cl_max=n_str_near(j)
            endif
         enddo
***   check the size of the cluster ------------->
         ratio=float(n_str_cl_max)/float(n_str)
         if(ratio.gt.ratio1.and.RMSD_cut.gt.RMSD_cut_min)then
            RMSD_cut=RMSD_cut-0.1
            goto 140
         endif
         if(ratio.lt.ratio2.and.RMSD_cut.lt.RMSD_cut_max)then
            RMSD_cut=RMSD_cut+0.2
            goto 140
         endif
         if(RMSD_cut.lt.RMSD_cut_min)RMSD_cut=RMSD_cut_min
         if(RMSD_cut.gt.RMSD_cut_max)RMSD_cut=RMSD_cut_max
      endif
***** RMSD_cut of first cluster, i.e. RMSD_cut_min, is decided ***
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccc find out nc_max clusters cccccccccccccccccccccccc      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nc=0                      !number of clusters
      do i=1,nc_max
***   calculate nc(j), number of neightboring structures ----------------->
         do j=1,n_str
            n_str_near(j)=0     !number of structure in cluster
            do k=1,n_str
               if(mark(j).eq.1.and.mark(k).eq.1)then
                  if(amat(j,k).lt.RMSD_cut)then
                     n_str_near(j)=n_str_near(j)+1
                  endif
               endif
            enddo
         enddo
***   find out the biggest cluster ----------------->
         n_str_cl_max=0         !maximum number of structures in cluster
         do j=1,n_str
            if(n_str_near(j).gt.n_str_cl_max)then
               n_str_cl_max=n_str_near(j)
               i_cl(i)=j        !structure center of i'th cluster
            endif
         enddo
***   check the size of the cluster ------------->
         if(n_str_cl_max.lt.1) goto 41 !no remaining clusters.
***   remove the structures in the cluster-------->
         n_str_cl(i)=0          !# of structure including some used struct
         n_str_cl_ex(i)=0       !# of structure excluding used structure
         R_in(i)=0              !average pairwise RMSD including
         R_ex(i)=0              !average pairwise RMSD excluding
         do j=1,n_str
c     if(mark(j).eq.1)then
            if(amat(j,i_cl(i)).lt.RMSD_cut)then
               if(mark(j).ne.0)then
                  n_str_cl_ex(i)=n_str_cl_ex(i)+1
                  R_ex(i)=R_ex(i)+amat(j,i_cl(i))
                  i_str_cl_ex(i,n_str_cl_ex(i))=j
               endif
               mark(j)=0
               n_str_cl(i)=n_str_cl(i)+1
               i_str_cl(i,n_str_cl(i))=j
               R_in(i)=R_in(i)+amat(j,i_cl(i))
            endif
c     endif
         enddo
         R_in(i)=R_in(i)/n_str_cl(i)
         R_ex(i)=R_ex(i)/n_str_cl_ex(i)
         rmsd_cl_cut(i)=rmsd_cut
         nc=i
c     write(*,*)i,n_str_cl_ex(i),n_str_cl_max
      enddo
c^^^^^^^^^^^^^^^^5 cluster find out ^^^^^^^^^^^^^^^^^^^^^^
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 41   continue

*****************************************************************
c     calculate 'combo.pdb' and E_combo
COMBO************************************************************
      do i=1,nc
         ns=0
         E_combo(i)=0
         rmsd_nat_a=0
         rmsd_cen_a=0
         do ii=1,Lch
            xs(ii)=0
            ys(ii)=0
            zs(ii)=0
         enddo
***
         k=i_cl(i)              !center structure of the cluster
         do ii=1,Lch
            r_2(1,ii)=x(ii,k)   !center structure
            r_2(2,ii)=y(ii,k)
            r_2(3,ii)=z(ii,k)
         enddo
         nn=n_str_cl(i)         !number of structures in cluster (including)
***
         write(21,1200)i,rmsd_str(k),istr(k),filen(itra(k))
 1200    format('#Cluster',i3,f8.3,i8,2x,a15)
         write(21,*)'i_cl   i_str  R_nat   R_cen  E    #str     traj'
         write(21,*)'-----------------------------------------------'
         write(21,*)'Nstr=',nn
         do l=1,nn
            m=i_str_cl(i,l)     !l'th structure in i'th cluster
            write(21,1201)i,l,rmsd_str(m),amat(k,m),E(m),
     $           istr(m),filen(itra(m))
 1201       format(i3,i8,2f7.3,f9.1,i8,2x,a15)
***   E_combo:
            E_combo(i)=E_combo(i)+E(m)
            ns=ns+1
            rmsd_nat_a=rmsd_nat_a+rmsd_str(m)
            rmsd_cen_a=rmsd_cen_a+amat(k,m)
***   rotate nearby structure into the center structure---->
            do n=1,Lch
               r_1(1,n)=x(n,m)
               r_1(2,n)=y(n,m)
               r_1(3,n)=z(n,m)
            enddo
            call u3b(w,r_1,r_2,Lch,1,rms,u,t,ier) !u rotate r_1 to r_2
            do j=1,Lch
               xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
               yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
               zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
               xs(j)=xs(j)+xx
               ys(j)=ys(j)+yy
               zs(j)=zs(j)+zz
            enddo
         enddo
***   averaging:
         do j=1,Lch
            xc(i,j)=xs(j)/float(ns)
            yc(i,j)=ys(j)/float(ns)
            zc(i,j)=zs(j)/float(ns)
         enddo
         E_combo(i)=E_combo(i)/float(ns)
         rmsd_nat_a=rmsd_nat_a/ns
         rmsd_cen_a=rmsd_cen_a/ns
         write(21,*)'-----------------------------------------'
         write(21,1222)rmsd_nat_a,rmsd_cen_a,E_combo(i)
 1222    format('   average=',2f7.3,f9.1)
      enddo
c^^^^^^^^^^^^^^^^^^^^ combo finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*****************************************************************
c     calculate <RMSD to centroid> & find closc: key_closc(i,k)
*****************************************************************
      do i=1,nc
         m_closc(i)=nclc        !top 10 closc recorded
         if(m_closc(i).gt.n_str_cl(i)) m_closc(i)=n_str_cl(i)
         do k=1,Lch
            r_1(1,k)=xc(i,k)
            r_1(2,k)=yc(i,k)
            r_1(3,k)=zc(i,k)
         enddo
         nn=0
         Rc_in(i)=0             !RMSD to centroid including
         do j=1,n_str_cl(i)
            m=i_str_cl(i,j)
            do k=1,Lch
               r_2(1,k)=x(k,m)
               r_2(2,k)=y(k,m)
               r_2(3,k)=z(k,m)
            enddo
            call u3b(w,r_1,r_2,Lch,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch) !RMSD12
            r_decoy_combo(j)=armsd
            Rc_in(i)=Rc_in(i)+armsd
            nn=nn+1
         enddo
         Rc_in(i)=Rc_in(i)/nn
         call SORT(n_str_cl(i),r_decoy_combo,key)
ccc   
         do j=1,m_closc(i)
            rmsd_cc(i,j)=r_decoy_combo(key(j))
            do k=1,Lch
               x_cc(i,j,k)=x(k,i_str_cl(i,key(j)))
               y_cc(i,j,k)=y(k,i_str_cl(i,key(j)))
               z_cc(i,j,k)=z(k,i_str_cl(i,key(j)))
            enddo
         enddo
ccc   
      enddo
      do i=1,nc
         do k=1,Lch
            r_1(1,k)=xc(i,k)
            r_1(2,k)=yc(i,k)
            r_1(3,k)=zc(i,k)
         enddo
         nn=0
         Rc_ex(i)=0             !RMSD to centroid including
         do j=1,n_str_cl_ex(i)
            m=i_str_cl_ex(i,j)
            do k=1,Lch
               r_2(1,k)=x(k,m)
               r_2(2,k)=y(k,m)
               r_2(3,k)=z(k,m)
            enddo
            call u3b(w,r_1,r_2,Lch,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch) !RMSD12
            Rc_ex(i)=Rc_ex(i)+armsd
            nn=nn+1
         enddo
         Rc_ex(i)=Rc_ex(i)/nn
      enddo
      do i=1,nc
         iden1(i)=0             ! which decoy at the trajectory file
         iden2(i)=0             ! which trajectory file
      enddo
      
*****************************************************************
c     calculate 'closc.pdb', closest structure to combo
CLOSC************************************************************
      if(n_closc.lt.0)goto 130  !use decoys in the same cluster
      do i=1,nc
         rmsd_close_min(i)=100.
      enddo
      do i_tra=1,n_tra
         i12=0
         open(10,file=filen(i_tra),status='old')
 110     read(10,*,end=119,err=119)Lch1,energ
         i12=i12+1
         do i=1,Lch1
            read(10,*,end=119,err=119)xt(i),yt(i),zt(i)
         enddo
         do i=1,Lch
            r_1(1,i)=xt(i)
            r_1(2,i)=yt(i)
            r_1(3,i)=zt(i)
         enddo
***   check whether this structure close to 'combo.pdb':
         do j=1,nc
            do i=1,Lch
               r_2(1,i)=xc(j,i)
               r_2(2,i)=yc(j,i)
               r_2(3,i)=zc(j,i)
            enddo
            call u3b(w,r_1,r_2,Lch,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch)
            if(armsd.lt.rmsd_close_min(j))then
               iden1(j)=i12     ! which decoy at the trajectory file
               iden2(j)=i_tra   ! which trajectory file
               rmsd_close_min(j)=armsd
               do i=1,Lch
                  xcl(j,i)=r_1(1,i)
                  ycl(j,i)=r_1(2,i)
                  zcl(j,i)=r_1(3,i)
               enddo
            endif
         enddo
***   
         goto 110
 119     close(10)
      enddo
ccc   check and replace the closest closc:
      do i=1,nc
         if(rmsd_close_min(i).lt.rmsd_cc(i,1))then
            do j=m_closc(i),2,-1
               rmsd_cc(i,j)=rmsd_cc(i,j-1)
               do k=1,Lch
                  x_cc(i,j,k)=x_cc(i,j-1,k)
                  y_cc(i,j,k)=y_cc(i,j-1,k)
                  z_cc(i,j,k)=z_cc(i,j-1,k)
               enddo
            enddo
            j=1
            rmsd_cc(i,j)=rmsd_close_min(i)
            do k=1,Lch
               x_cc(i,j,k)=xcl(i,k)
               y_cc(i,j,k)=ycl(i,k)
               z_cc(i,j,k)=zcl(i,k)
            enddo
         endif
      enddo
 130  continue
c^^^^^^^^^^^^^^^^^^^^closc finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
ccccccccccccccoutput the 5 clusters in PDB format ccccccccccccccc
      do i=1,Lch
         read(2,*)ii,seq(i)
      enddo
      do i=1,nc
         if(i.lt.10)then
            c2=char(48+i)//'.pdb'
         else
            c2=char(48+i/10)//char(48+i-10*(i/10))//'.pdb'
         endif
***   'combo.pdb'->average of structures of cluster-------->
         open(10,file='combo'//c2,status='unknown')
         do j=1,Lch
            xx=xc(i,j)
            yy=yc(i,j)
            zz=zc(i,j)
            write(10,1237)'ATOM',j,'CA',seq(j),j,'',xx,yy,zz
         enddo
         write(10,322)
 322     format('TER')
         do j=2,Lch
            write(10,1238)'CONECT',j-1,j
         enddo
         close(10)
***   'closc.pdb'->structure closest to combo.pdb-------->
         open(10,file='closc'//c2,status='unknown')
         do j=1,Lch
            xx=x_cc(i,1,j)
            yy=y_cc(i,1,j)
            zz=z_cc(i,1,j)
            write(10,1237)'ATOM',j,'CA',seq(j),j,'',xx,yy,zz
         enddo
         write(10,322)
         do j=2,Lch
            write(10,1238)'CONECT',j-1,j
         enddo
         close(10)
***   'closci_j.pdb'->structure closest to combo.pdb-------->
         do j=1,m_closc(i)
            if(i.lt.10)then
               c3=char(48+i)
            else
               c3=char(48+i/10)//char(48+i-10*(i/10))
            endif
            if(j.lt.10)then
               c4=c3(1:len_trim(c3))//'_'//char(48+j)//'.pdb'
            else
               c4=c3(1:len_trim(c3))//'_'//
     &              char(48+j/10)//char(48+j-10*(j/10))//'.pdb'
            endif
            open(10,file='closc'//c4,status='unknown')
            do k=1,Lch
               xx=x_cc(i,j,k)
               yy=y_cc(i,j,k)
               zz=z_cc(i,j,k)
               write(10,1237)'ATOM',k,'CA',seq(k),k,'',xx,yy,zz
            enddo
            write(10,322)
            do k=2,Lch
               write(10,1238)'CONECT',k-1,k
            enddo
            close(10)
         enddo
      enddo
 1237 format(A4,I7,A4,A5,I6,A4,3F8.3)
 1238 format(A6,i5,i5)
c^^^^^^^^^^^^^^output PDB structures finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      write(20,*)'Target: ',protein
      write(20,*)'Modeling Length: ',Lch
      write(20,*)'Native Length: ',Lch_n
      write(20,*)'Number of clusters: ',nc
      write(20,*)'Total number of structures=',n_str_all
      write(20,*)'Number of structure in use=',n_str
      write(20,*)
cccccccRMSD and TM-score of cluster to native cccccccccccccccccccccccccc
      if(Lch_n.gt.1)then
         write(20,*)'--- comparison to native structure-------'
         write(20,*)'  i  R_combo  TM_combo  R_closc  TM_closc'
         write(20,*)'A----------------------------------------'
         do i=1,nc
            do j=1,Lch_n
               r_1(1,j)=x_n(j)  !CA
               r_1(2,j)=y_n(j)
               r_1(3,j)=z_n(j)
               nn2(j)=j
               x2(j)=x_n(j)
               y2(j)=y_n(j)
               z2(j)=z_n(j)
            enddo
***   combo->
            do j=1,Lch_n
               r_2(1,j)=xc(i,ires(j))
               r_2(2,j)=yc(i,ires(j))
               r_2(3,j)=zc(i,ires(j))
               nn1(j)=j
               x1(j)=xc(i,ires(j))
               y1(j)=yc(i,ires(j))
               z1(j)=zc(i,ires(j))
            enddo
            call u3b(w,r_1,r_2,Lch_n,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch_n) !RMSD12
            rmsd_combo=armsd
            call TMscore(Lch_n,x1,y1,z1,nn1,Lch_n,x2,y2,z2,nn2,
     &           TM1,Rcomm,Lcomm)
***   closc->
            do j=1,Lch_n
               r_2(1,j)=xcl(i,ires(j))
               r_2(2,j)=ycl(i,ires(j))
               r_2(3,j)=zcl(i,ires(j))
               nn1(j)=j
               x1(j)=xcl(i,ires(j))
               y1(j)=ycl(i,ires(j))
               z1(j)=zcl(i,ires(j))
            enddo
            call u3b(w,r_1,r_2,Lch_n,0,rms,u,t,ier)
            armsd=dsqrt(rms/Lch_n) !RMSD12
            rmsd_closc=armsd
            call TMscore(Lch_n,x1,y1,z1,nn1,Lch_n,x2,y2,z2,nn2,
     &           TM2,Rcomm,Lcomm)
***   
            write(20,50)i,rmsd_combo,TM1,rmsd_closc,TM2
 50         format(i5,f9.3,f9.4,f9.3,f9.4)
         enddo
         write(20,*)
         txt1='                         '
         txt2='RMSD   Energy   #str    Trajectory'
         write(20,*)txt1(1:23)//txt2(1:35)
         write(20,*)'-------------------------------------------'
         i=itra_rmsd_min_all
         j=itra_rmsd_min_use
         k=itra_E_min
         l=itra_E_min_all
         write(20,124)rmsd_min_all,E_rma,istr_rmsd_min_all,filen(i)
         write(20,125)rmsd_min_use,E_rmu,istr_rmsd_min_use,filen(j)
         write(20,126)rmsd_E_min,E_min,istr_E_min,filen(k)
         write(20,127)rmsd_E_min_all,E_min_all,istr_E_min_all,filen(l)
 124     format('Minimum RMSD in all=',f8.3,f9.1,i8,2x,a20)
 125     format('Minimum RMSD in use=',f8.3,f9.1,i8,2x,a20)
 126     format('Minimum E    in all=',f8.3,f9.1,i8,2x,a20)
 127     format('Minimum E    in use=',f8.3,f9.1,i8,2x,a20)
      endif
      write(20,*)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
ccccccccccccccccccoutput cluster analysis cccccccccccccccccccc
      write(20,*)'------------ summary of clusers -------------------'
      txt1='i Size R_cut density R_cl_co '
      txt2='<E>    E_model  #str(mod) Traj(mod)'
      txt3='  #str(clo)  Traj(clo)'
      write(20,*)txt1(1:31)//txt2(1:35)//txt3(1:35)
      write(20,*)'B--------------------------------------------------'
      do i=1,nc
         k=i_cl(i)
         density=n_str_cl(i)/rmsd_cl_cut(i)
         if(iden1(i)>0)then
            write(20,123)i,n_str_cl(i),rmsd_cl_cut(i),density,
     $           rmsd_close_min(i),E_combo(i),E(k),istr(k),
     $           filen(itra(k)),iden1(i),filen(iden2(i))
         else
            write(20,123)i,n_str_cl(i),rmsd_cl_cut(i),density,
     $           rmsd_close_min(i),E_combo(i),E(k),istr(k),
     $           filen(itra(k))
         endif
      enddo
 123  format(i2,i6,f6.2,f7.0,f6.2,2f9.1,i6,1x,a13,i6,1x,a13)
      
ccccccccccccccccccoutput cluster analysis cccccccccccccccccccc
      write(20,*)
      write(20,*)' include used structure  exclude used structre'
      write(20,*)'  ---------------------  ---------------------'
      write(20,*)'i  N_in  <R_in> <Rc_in>   N_ex  <R_ex> <Rc_ex>'
      write(20,*)'C---------------------------------------------'
      do i=1,nc
         write(20,133)i,n_str_cl(i),R_in(i),Rc_in(i),
     $        n_str_cl_ex(i),R_ex(i),Rc_ex(i)
      enddo
 133  format(i2,i6,f7.2,f7.2,i8,f7.2,f7.2)
      
      write(20,*)
      write(20,*)'RMSD_cut_initial=',RMSD_cut_initial
      write(20,*)'RMSD_cut_min=',RMSD_cut_min
      write(20,*)'RMSD_cut_max=',RMSD_cut_max
      write(20,*)'ratio_min=',ratio1
      write(20,*)'ratio_max=',ratio2
      write(20,*)
      write(20,*)'trajectories used:',n_tra
      do i=1,n_tra
         write(20,*)filen(i)
      enddo
c^^^^^^^^^^^^^^^^^^^^output done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      stop
      end
      
cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m           (given)
c  x    - x(i,m) are coordinates of atom c m in set x       (given)
c  y    - y(i,m) are coordinates of atom c m in set y       (given)
c  n    - n is number of atom pairs                         (given)
c  mode  - 0:calculate rms only                             (given)
c          1:calculate rms,u,t                              (takes longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
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
      
      double precision sqrt3, tol, zero
      
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
      
      wc = zero
      rms = zero
      e0 = zero
      
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
         if( w(m) .lt. 0.0 ) return
         wc = wc + w(m)
         do i=1, 3
            xc(i) = xc(i) + w(m)*x(i,m)
            yc(i) = yc(i) + w(m)*y(i,m)
         end do
      end do
      if( wc .le. zero ) return
      do i=1, 3
         xc(i) = xc(i) / wc
         yc(i) = yc(i) / wc
      end do
      
      do m=1, n
         do i=1, 3
            e0=e0+w(m)*((x(i,m)-xc(i))**2+(y(i,m)-yc(i))**2)
            d = w(m) * ( y(i,m) - yc(i) )
            do j=1, 3
               r(i,j) = r(i,j) + d*( x(j,m) - xc(j) )
            end do
         end do
      end do
      
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
      
      rms = (e0 - d) - d
      if( rms .lt. 0.0 ) rms = 0.0
      
      return
      end

*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*     Reference: Yang Zhang, Jeffrey Skolnick, Proteins 2004 57:702-10.
*
*     Explanations:
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
      subroutine TMscore(L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,TM,Rcomm,Lcomm)
      PARAMETER(nmax=2000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ****************
      nseqA=L1
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
      enddo
      nseqB=L2
      do i=1,L2
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
      enddo

******************************************************************
*     pickup the aligned residues:
******************************************************************
      k=0
      do i=1,nseqA
         do j=1,nseqB
            if(nresA(i).eq.nresB(j))then
               k=k+1
               iA(k)=i
               iB(k)=j
               goto 205
            endif
         enddo
 205     continue
      enddo
      n_ali=k                   !number of aligned residues
      Lcomm=n_ali
      if(n_ali.lt.1)then
c         write(*,*)'There is no common residues in the input structures'
         TM=0
         Rcomm=0
         return
      endif

************/////
*     parameters:
*****************
***   d0------------->
      d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      if(d0.lt.0.5)d0=0.5
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
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
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
      PARAMETER(nmax=2000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score,score_max

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

*********************************************************
*     sort A(N) and give KEY, coded by Y Zhang 2013/12/7
*     input: N, A(N); output: KEY(N), smallest is first
*     A(N) is unchanged before and after sorting
*********************************************************
      SUBROUTINE SORT(N,A,KEY)
      dimension A(N),KEY(N)
      dimension B(N),mk(N)
      
ccc   initialization------>
      do i=1,N
         B(i)=A(i)
         mk(i)=-1
      enddo
      
ccc   sort B(N) --------->
      DO 30 I=2,N
         X=B(I)
         J=I
 10      J=J-1
         IF(J.EQ.0 .OR. B(J).LE.X) GO TO 20
         B(J+1)=B(J)
         GO TO 10
 20      B(J+1)=X
 30   CONTINUE
      
ccc   find KEY(N) --------->
      do i=1,N
         do j=1,N
            if(mk(j).eq.-1)then
               if(abs(A(j)-B(i)).lt.0.00000001)then
                  KEY(i)=j
                  mk(j)=1
                  goto 40
               endif
            endif
         enddo
 40      continue
      enddo
      
      END
