
*****************************************************************************
*     This program is to generate protein structural decoys by Monte Carlo  *
*     simulations under an on-and-off-lattice CAS model. It is illigal to   *
*     distribute any part of this code without writing permission from the  *
*     author. Please address comments/bug-reports to: zhng@umich.edu        *
*****************************************************************************
*
*     This program should be compiled by 
*     'gfortran -static -O3 -ffast-math -lm -o cas cas.f'
*
*     Last update by yzhang on March 18 2013
*

c        1         2         3         4         5         6         7 !
c 3456789012345678901234567890123456789012345678901234567890123456789012345678
      program TASSER
      implicit integer(i-z)
      parameter(ndim=1999)      !maximum length of chain-length
      parameter(nrep=100)       !maximum number of replicas
      parameter(nvec=416)       !number of vectors on lattice
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/three/angle(nvec,nvec)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b
      common/hopp/eonehw(0:19)
      common/looks/exc,exc1,exc2
      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc

      common/arandom/  aarand,abrand,acrand,adrand
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12) !safe when vr^2<30
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/fr/frga(ndim),frgb(ndim)
      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/maxi/maxin,vect1,vect2
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/forpreparemove4/ asrr(0:19,0:19,-12:12)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/distres/er4,es3c
      common/rmsdrange/nca1,nca2
      common/CA/dx(ndim),dy(ndim),dz(ndim)
      common/msichores/msicho
      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,EHB5,EHB6
      common/ehbenergy1/EHB5a,EHB5b
      common/eshortenergy1/ESHORT1,ESHORT2,ESHORT3,ESHORT4,ESHORT11
      common/eshortenergy2/ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c
      common/eshortenergy3/ESHORT6,ESHORT7,ESHORT8,ESHORT9,ESHORT10
      common/eshortenergy4/ESHORT12
      common/otherenergy/E_cord,E_cnum
      common/resnumber/Ncom,Ndis,accur
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/thrtem/aT_ann(100)
      common/outputxyz/fxyz(3,ndim)

      common/temperature/itemp,atemp
      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)
      common/pair1/eh2,eh1b,eh1c
      dimension E_s(nrep),E_ss(nrep)
      common/paircut/ash
      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
      common/countres1/t_combCA8,t_distL,s_combCA8,s_distL,N_resc

      common/nswap/bNSa(100),bNSt(100)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/pairmap2/apabla(0:19,0:19),apablp(0:19,0:19)
      common/pairmap3/apablm(0:19,0:19),amap(ndim,ndim),nmap

      common/moveretio1/h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,h7,hend
      common/moveretio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
      common/moveretio3/bh7a,bh7b,bhendn,bhendc
      common/icgg/ icg(ndim), EH6  
      common/rs/i_thr0
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      common/nrepfile/n_repf

      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/mng/m_g(100)
      common/acct/accept0
      character*6 mname
      character fn
      common/movename/mname(100)
      common/readinitial/m_initial
      common/eall/N_sum(100),energ_sum(100),energ_sum2(100),E_min

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/chainm/mv(ndim)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep) !Ca
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) !SG
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) !cc
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) !hb
      common/weight/chuan
      common/bigbond/i_bigbond,teco
      common/ssp/ssp
      common/defoangle/defo_angle
      common/fractpair/fract_pair1,fract_pair3
      common/zscore/izscore
      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs
      common/ichos/ichos
      common/nana1/nana
      common/ranzy/nozy
      common/hours/hour_max
      common/stick1/nstick,astick,nrmsd,ermsd
      common/stick2/iq(ndim,nrep)
      common/stick3/ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
      common/stick6/bx00(ndim),by00(ndim),bz00(ndim),armsd_min0
      common/trackn/n_tem(100)
      dimension eshort12_a(nrep)
      common/lattice/m_latt,latt1,latt2
      common/nwmax1/nvecr
      common/iter/n_run
      common/aminoacid/sequ(ndim)
      character*3 sequ
      common/stick7/itemp0,icycle,icycle0
      common/mloopf/mloop
      common/res2/er14,er15,er16,er17
      common/svm1/Mcon(3,ndim),Kcon(3,ndim,300),awei(3,ndim,ndim),fw
      common/svm2/acc_cut,dist_svm2(15),nk_svm(3),ik_svm(3,15)
      common/svm3/er21,er22,er23
      common/svm4/acc_cutB,acc_cutG
      common/svm5/acc_cut0,aw_cont,aw_cont0,aLo1,aLo2,aLo3,anf,itype(10)
      common/svm6/npot,dist_svm(15),dist_svmU(15),dist_svmU2(15)
      common/concuttu/npot1,concutU2(0:19,0:19),concutU(0:19,0:19),aw1
      common/concuttu4/npot4,aw4,Cr20
      common/CAc1/npot2,d_CA_cut,d_CA_cut2,d_CA_cutU,d_CA_cutU2,aw2
      common/CA8a/npot3,e_CA_cut,e_CA_cut2,e_CA_cutU,e_CA_cutU2,aw3
      common/CA8b/eq2a,eq2b,eq2c,eq2d
      common/CA8c/dcut3
      common/CB6a/npot5,aw5,f_CB_cut,f_CB_cut2,f_CB_cutU,f_CB_cutU2
      common/CB6b/fq2a,fq2b,fq2c,fq2d
      common/CB6c/McomCB6(ndim),KcomCB6(ndim,100),aweigCB6(ndim,ndim)
      common/CB8a/npot6,aw6,g_CB_cut,g_CB_cut2,g_CB_cutU,g_CB_cutU2
      common/CB8b/gq2a,gq2b,gq2c,gq2d
      common/CB8c/McomCB8(ndim),KcomCB8(ndim,100),aweigCB8(ndim,ndim)
      common/dwell/dwell,ha(15),hb(15),hc(15),hd(15)
      
      character*20 cfile(100)
      common/cfile/cfile,n_conf
      dimension awei0(15,ndim,ndim)
      
cccc  RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   

***************************************************
*     three types of contacts:
*
*     1, 'quasi3.comm', generic, dmin(A,A)<d<dmax(A,A) from 'quasi3.comm'
*        'profile3.comm', generic environment, d<dmax(A,A) from 'quasi3.comm'
*     
*     2, 'pair3.dat', prot-spec, merge to 'quasi3', dmin(A,A)<d<dmax(A,A) from 'quasi3.comm'
*        'pair1.dat', prot-spec, to pair3 but orien-indp, d<dmax(A,A) from 'contact.comm'
*     
*     3, 'comb.dat', SG from 'init.dat', d<dmax(A,A) from 'contact.comm'
*        'combCA.dat', CA from 'init.dat', d<6A
*        'comb8CA.dat', CA from 'init.dat', d<8A
*        'par.dat', SG from 'init.dat', merge to 'pair1.dat', d<dmax(A,A) from 'contact.comm'
***************************************************

ccccccccccccccccccccccccc common input files cccccccccccccccccccccccccccccc
      open(unit=1,file='contact.comm', status='old') !cutoff for 'comb.dat' and 'pair1.dat'
      open(unit=2,file='profile3.comm',status='old') !contact_envir, d_max from 'quasi3.comm'
      open(unit=3,file='quasi3.comm',  status='old') !SG-potential, [d_min,d_max]
      open(unit=4,file='sidechain.comm',status='old') !for Sc position
      open(unit=5,file='r13.comm',    status='old') !E_short of (i,i+2)
      open(unit=12,file='r14.comm',   status='old') !E_short of (i,i+3)
      open(unit=7,file='r14h.comm',   status='old') !stress helical-stru.
      open(unit=8,file='r14e.comm',   status='old') !stress extended-stru.
      open(unit=9,file='r15.comm',   status='old') !E_short(i,i+4)
      open(unit=10,file='r15h.comm',   status='old')
      open(unit=11,file='r15e.comm',   status='old')

ccccccccccccccc sequence specified input files cccccccccccccccc///////////
      open(unit=14,file='seq.dat',     status='old') !note format
      open(unit=15,file='par.dat',   status='old') !SG-pair potential
      open(unit=16,file='comb.dat',   status='old') !SG-contact, d_cut from contact.comm
      open(unit=17,file='dist.dat',    status='old') !CA-distant restraints
      open(unit=18,file='combCA.dat',status='old') !CA-contact restraints
c      open(unit=28,file='comb8CA.dat',status='old') !CA-contact at 8A <-not used
c      open(unit=29,file='comb6CB.dat',status='old') !CB-contact at 6A
c      open(unit=30,file='comb8CB.dat',status='old') !CB-contact at 8A
      open(unit=21,file='rmsinp',status='old') !chain length for Lch
      open(unit=22,file='distL.dat',status='old') !long range CA dist restraint
      open(unit=25,file='pair3.dat',status='unknown') !contact_ori,[d_min,d_max] from quasi3.comm
      open(unit=26,file='pair1.dat',status='unknown') !no_ori, d_cut from contact.comm
      open(unit=27,file='exp.dat',status='unknown')
      open(unit=24,file='init.dat',status='unknown')
      open(unit=50,file='contact.map',status='unknown') !unified contact map
      
      open(unit=19, file='in.dd',     status='old')
      open(unit=20, file='out.d',     status='unknown')
c      open(unit=93, file='potential.d', status='unknown')
ccccccccccccccccfor E-t cccccccccccccccccccccccccccccccccccccccccccccccccc
c     open(unit=91, file='swepa.d',    status='unknown') !$$
c     open(unit=92, file='swepb.d',    status='unknown') !$$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      read(19,*) hour_max
      read(19,*) random,ncycle,phot,N_rep,n_run
      read(19,*) h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend
      read(19,*) atemp2,atemp1,exc,exc1,exc2,Mend,defo_angle
      read(19,*) d_xyz0,angle0,L_cut,teco
      read(19,*) switch,i_thr0,ssp !i_thr0: 0->consensus; 1->top-1; 2->2th
      read(19,*) eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh !alpha-type HB
      read(19,*) eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh !alpha-type HB
      read(19,*) eh1,eh2,eh3,eh4 !for EHBs (6)
      read(19,*) eh1a,eh1b
      read(19,*) es2,es3,es4,es5,es6     !for ESHORTs (6)
      read(19,*) es3a,es3b,es3c
      read(19,*) en1,en2,en3    !for ensemble (3)
      
      read(19,*) chuan,nana,m_latt,ifa !decide if we should increase/decreas weight
      read(19,*) aTs2,aTs1    !SQ2
      read(19,*) aTTs2,aTTs1  !ARCSH
      read(19,*) nstick,astick,nrmsd,ermsd !whether stick to templates
c      read(19,*) npot1,aw1,npot2,aw2,npot4,aw4,Cr20 !comb,combCA,par/pair1
      read(19,*) aw1,aw2,aw4,Cr20 !comb,combCA,par/pair1
      read(19,*) er6,dwell,fw   !weight for distL
      
***   
      write(20,*)'starting time: ',fdate() !pgf77 has problem on fdate()
      write(*,*)'starting time: ',fdate()
***   
      
      call set_common           !set common parameters
      call read_seq             !read seq(i),sec(i) from 'seq.dat'
      call read_centro          !read eoinp(i,dis) from 'centro.comm'
      call read_profile         !read envir(ia,im,ip,i,j) from 'profile3.comm'
      call read_E13             !read 1-3 short-range E from 'r13.comm'
      call read_E14             !read 1-4 potential from 'r14*.dat'
      call read_E15             !read 1-5 potential from 'r15*.dat'
      call read_quarsi3         !read 'quarsi3.comm'
      call read_par             !read 'par.dat', pair-wise potential
      call read_concut          !read cut-off for contact prediction
      call read_contactrestrain !read contact restrains from 'comb.dat'
      call read_distantrestrain !read distant restrains from 'dist.dat'
      call read_longdistantrestrain !read long distant restrain from 'distL.dat'
      call read_exp             !read slovent expose prediction
      call read_combCA          !read CAcontact restrains from 'combCA.dat'
c      call read_comb8CA         !read CAcontact restrains from 'comb8CA.dat'
c      call read_comb6CB         !read CBcontact restrains from 'comb6CB.dat'
c      call read_comb8CB         !read CBcontact restrains from 'comb8CB.dat'
      call read_seqcontact      !read seq_based_contact 'svmseqca6.dat'
      call reset_temperature    !reset temperature according N_rest
      call set_temperature      !set temperature for different replic
      call set_EHB              !set structure-specitic H-bond, EHBIJ(i,j)
      
      call prepare_vectors      !prepare all possible bond-vectors
      call prepare_neighbors    !define goodc(i,j), angle(i,j), prod(i,j)
      call prepare_beta         !define C_beta, C_group, and hydrogen-bond
      call prepare_frg          !compute the secondary fragment biases

      call get_acorder          !calculate contact order
      call write_parameter      !print out initial parameters

      call set_move_retio       !set movement percentage
      call prepare_move2        !2-bond move, num2=26784
      
      n_repf=5                  !number of output replicas
      if(izscore.eq.1)then      !easy target
         n_repf=8               !number of output replicas
      endif
      if(switch.gt.1)then
         n_repf=3               !number of output replicas
         call template_initial  !initial model from templates
      endif
      
ccccccccccccccc trajectory files cccccccccccccccccccccccccccccccccccccc
      if(n_repf.gt.N_rep)n_repf=N_rep
      do i=1,n_repf
         if(i.lt.10)then
            fn=char(48+i)
            open(unit=30+i,file='rep'//fn//'.tra',status='unknown')
         else
            fn=char(48+(i-10))
            open(unit=30+i,file='rep1'//fn//'.tra',status='unknown')
         endif
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(switch.gt.1)then
         call template_simulation !simulate the templates
         stop
      endif

      do i=1,100
         bNSa(i)=0              !aceptance for swep
         bNSt(i)=0

         bNa(i)=0               !acceptance for move2,3,4,5,6,7
         bNt(i)=0

         bNNa(i)=0              !acceptance for different temperature.
         bNNt(i)=0

         N_sum(i)=0
         energ_sum(i)=0         !<E_tot>
         energ_sum2(i)=0        !<E_tot^2>
         eshort12_a(i)=0
      enddo
      i_tr=0                    !order number of output trajectory
      E_min=10000
      mcycle=0
      
      i_run=0
 15   i_run=i_run+1
c     call random_initial       !produce initial structure randomly.
      call read_initial         !read initial (x,y,z) from 'init.dat'

      do i=1,n_rep
         n_tem(i)=i             !i'th replica from n_tem(i)'th template
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc       The main cycle start from here !                         ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do k=1,3                  !CA,CB,SG
         do i=1,Lch
            do j1=1,Mcon(k,i)   !number of contacts with 'i'
               j=Kcon(k,i,j1)   !residue-ID contact to 'i'
               awei0(k,i,j)=awei(k,i,j)
            enddo
         enddo
      enddo
      i_an=-1
      n_an=20                   !each annealing contains 20 moves
      do 1111 icycle=1,ncycle
cccc  
         if(ifa.eq.1)then       !<--increase weight gradually
            i_an=i_an+1
            afa=float(i_an)/float(n_an)
            if(i_an.ge.n_an)then
               i_an=-1
            endif
            do k=1,3
               do i=1,Lch
                  do j1=1,Mcon(k,i) !number of contacts with 'i'
                     j=Kcon(k,i,j1) !residue-ID contact to 'i'
                     awei(k,i,j)=awei0(k,i,j)*afa
                  enddo
               enddo
            enddo
c            write(*,*)'afa=',icycle,afa
         endif
         if(ifa.eq.-1)then      !<--decrease weight gradually
            i_an=i_an+1
            afa=1-float(i_an)/float(n_an)
            if(i_an.ge.n_an)then
               i_an=-1
            endif
            do k=1,3
               do i=1,Lch
                  do j1=1,Mcon(k,i) !number of contacts with 'i'
                     j=Kcon(k,i,j1) !residue-ID contact to 'i'
                     awei(k,i,j)=awei0(k,i,j)*afa
                  enddo
               enddo
            enddo
c            write(*,*)'afa=',icycle,afa
         endif
ccccc 
         do 2222 itemp=1,N_rep  !iterate for all the replicas
            atemp=aT_rep(itemp)	!current temperature
            aTs=aTs_rep(itemp)
            aTTs=aTTs_rep(itemp)
            call set_current    !get current (x,y,z,ica)
            call initial_move   !update center, axis, energy
ccc
            do 3333 iphot=1,phot !N_swap, iterate at fixed temperature
               do 4444 i_lch=1,Lch
                  fff=aranzy(nozy)
                  if(fff.le.bh2)then
                     call move2
                  elseif(fff.le.bh3s)then
                     call move3s
                  elseif(fff.le.bh3d)then
                     call move3d
                  elseif(fff.le.bh4s)then
                     call move4s
                  elseif(fff.le.bh4d)then
                     call move4d
                  elseif(fff.le.bh5s)then
                     call move5s
                  elseif(fff.le.bh5d)then
                     call move5d
                  elseif(fff.le.bh6)then
                     call move6
c                     call move7a
c                     call move7b
c                     call move9
                  elseif(fff.le.bhendn)then
                     call move_n_end
                  else
                     call move_c_end
                  endif
                  atime=second()/3600.0 !pgf77:etime(tarray),real tarray(2)
                  if(atime.gt.hour_max)goto 901
 4444          continue
 3333       continue
ccc   
ccccccccccrecord energy and (x,y,z) cccccccccccccccccccc
            E_rep(itemp)=energy_tot() !whole energy
c            energy_total=
c     $            eh1*EHB1       !general soft-core energy for Ca-SC
c     $           +eh1a*EHB1a    !general soft-core energy for CA-CA 
c     $           +eh1b*EHB1b    !general soft-core energy for SC-SC 
c     $           +eh1c*EHB1c    !pair-wise potential for SC-SC 
c     $           +eh2*EHB2      !soft-core energy for SC-SC (quarsi3)
c     $           +eh3*EHB3      !coupling of secondary structure and pairwise
c     $           +eh4*EHB4      !enhanced quarsi3
c     $           +eh5a*EHB5a    !H-bond (alpha)
c     $           +eh5b*EHB5b    !H-bond (beta)
c     $           +es2*ESHORT2   !bury potential for SC
c     $           +er1*ESHORT3   !distmap
c     $           +er3*ESHORT4   !contact restrain
c     $           +er4*ESHORT4a  !deviation of contact restrain
c     $           +er5*ESHORT9   !CAcontact restrain
c     $           +er6*ESHORT10   !longrange CA-dist restrain
c     $           +er7*ESHORT11   !derivation to template
c     $           +astick*ESHORT12   !derivation to ax00
c     $           +es3*ESHORT5   !general bias to protein-like structure
c     $           +es3a*ESHORT5a  !panality on crumpling
c     $           +es3b*ESHORT5b  !bias to predicted alpha/beta structures
c     $           +es3c*ESHORT5c  !bias to possible alpha/beta structures
c     $           +es4*ESHORT6   !R13
c     $           +es5*ESHORT7   !R14
c     $           +es6*ESHORT8   !R15
c     $           +en1*eprofo    !E_environment
c     $           +en2*E_cord    !Contact order
c     $           +en3*E_cnum    !Contact number
c     parameters with input: eh1,eh1a,eh1b,eh2,eh3,eh4,eh5a,eh5b,
c                            es2,es3,es3a,es3b,es3c,es4,es5,
c                            en1,en2,en3
c     parameters mandtoried: eh1c,er1,er2,er3,er4,er5,er6,er7
c            write(*,*)E_rep(itemp),energy_total
c            write(*,*)icycle,itemp,E_rep(itemp),eshort12,'=='
            if(E_rep(itemp).lt.E_min) E_min=E_rep(itemp)
            do i=1,Lch
               xrep(i,itemp)=x(i)
               yrep(i,itemp)=y(i)
               zrep(i,itemp)=z(i)
            enddo
            eshort12_a(itemp)=eshort12_a(itemp)+eshort12
 2222    continue
         
cccccccccccccccccc print out 'swep.d' cccccccccccccccccccccccccccc
c         write(91,91)icycle,(E_rep(i),i=1,20)       !$$
c         write(92,91)icycle,(E_rep(i),i=21,N_rep)   !$$
c 91      format(i10,21f9.1)

ccccccccccccccccc<RMSD>, <E> cccccccccccccccccccccccccc
         do i=1,N_rep
            energ_sum(i)=energ_sum(i)+E_rep(i)
            energ_sum2(i)=energ_sum2(i)+E_rep(i)*E_rep(i)
            N_sum(i)=N_sum(i)+1
         enddo
         
ccccccccccccccccccccc snapshots of E(1), E(N_rep) ccccccccccccc
         if(icycle.eq.icycle/1*1)then
            i_tr=i_tr+1
            do k=1,n_repf
               write(30+k,401)Lch,E_rep(k),i_tr,icycle
               do i=1,Lch
                  abx=xrep(i,k)*0.87
                  aby=yrep(i,k)*0.87
                  abz=zrep(i,k)*0.87
                  write(30+k,402)abx,aby,abz
               enddo
            enddo
         endif
 401     format(i8,1x,f10.1,2i8)
 402     format(f10.3,1x,f10.3,1x,f10.3)

         call count_restrains   !count number of satisfied restraints
         
ccccccccccccccccc swap replicas cccccccccccccccccccccccccccccccc
         if(icycle.eq.icycle/2*2)then
            do i=1,N_rep-1,2    !swap odd replicas
               call swap(i,i+1)
            enddo
         else
            do i=2,N_rep-1,2
               call swap(i,i+1) !swap even replicas
            enddo
          endif
          mcycle=mcycle+1
 1111 continue
      if(i_run.lt.n_run)goto 15
 901  continue
c--------------------------Main cycle ended here!!---------------

      call test_neighbor        !check all the neighboring residues
      call test_overlap         !test the overlap of C_a and C_b

      energy_tot_tmp=energy_tot()
      write(20,*)'E_final=',energy_tot_tmp

      write(20,*)
      write(20,*)'<s_comb>=',s_comb/float(N_resc)
      write(20,*)'<t_comb>=',t_comb/float(N_resc)
      write(20,*)'s_comb/t_comb=',float(s_comb)/(t_comb+0.001)
      write(20,*)
      write(20,*)'<s_dist>=',s_dist/float(N_resc)
      write(20,*)'<t_dist>=',t_dist/float(N_resc)
      write(20,*)'s_dist/t_dist=',float(s_dist)/(t_dist+0.001)
      write(20,*)
      write(20,*)'<s_distL>=',s_distL/float(N_resc)
      write(20,*)'<t_distL>=',t_distL/float(N_resc)
      write(20,*)'s_distL/t_distL=',float(s_distL)/(t_distL+0.001)
      write(20,*)
      write(20,*)'<s_combCA>=',s_combCA/float(N_resc)
      write(20,*)'<t_combCA>=',t_combCA/float(N_resc)
      write(20,*)'s_combCA/t_combCA=',float(s_combCA)/(t_combCA+0.001)
      write(20,*)
      write(20,*)'<s_combCA8>=',s_combCA8/float(N_resc)
      write(20,*)'<t_combCA8>=',t_combCA8/float(N_resc)
      write(20,*)'s_combCA8/t_combCA8=',
     &     float(s_combCA8)/(t_combCA8+0.001)
      write(20,*)
      
cccccc output 'stick.pdb' ccccccccccccccccccccccccccccccc
c     sticki.pdb is the conformation closest to i'th structure in init.pdb
      if(nrmsd.eq.1)then
        write(20,*)' ------- RMSD to templates ------------'
        write(20,*)'rmsd_min0=',armsd_min0*0.87
        write(20,*)'icycle0=',icycle0
        write(20,*)'itemp0=',itemp0
        
        open(unit=70,file='stick.pdb',status='unknown')
        do i=1,Lch
          write(70,1037)i,sequ(i),i,bx00(i)*0.87,by00(i)*0.87,
     &      bz00(i)*0.87
        enddo
        write(70,*)'TER'
        close(70)
      endif
 1037 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
      
cccccccccccccccccccccccc Na/Nt cccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*) 'i_move    move   Na(i)  Nt(i)   Na(i)/Nt(i)'
      do i=2,15
         if(bNt(i).gt.1)then
            write(20,5004) i,mname(i),bNa(i),bNt(i),bNa(i)/bNt(i)
         else
            write(20,5004) i,mname(i),bNa(i),bNt(i)
         endif
      enddo
 5004 format(I4,A9,2f15.1,f11.6)
      
ccccccccccccccccccccccccccE_final, NSa/NSt ccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'-------------- E_final, Na_swap/Nt_swap ---------'
      WRITE(20,*) 'i  T(i) final_E(i)  Nsa(i)  Nst(i)  Nsa(i)/Nst(i)'
      do i=1, n_rep
         if(bNSt(i).gt.1)then
            write(20,5005) i,aT_rep(i),E_rep(i),
     $           bNSa(i),bNSt(i),bNSa(i)/bNSt(i)
         else
            write(20,5005) i,aT_rep(i),E_rep(i),bNSa(i),bNSt(i)
         endif
      enddo
 5005 format(I4,f7.2,f8.1,2f15.1,f11.6)

ccccccccccccccccccccccc <E>, NNa/NNt ccccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'------------ <energy>, Na(i)/Nt(i) ----------------'
      write(20,*)'i_rep  T   <E>  c  NNa(i_temp)  NNt(i_temp)  Na/Nt
     &     <ESHORT12>'
      do i=1,N_rep
         energ_sum(i)=energ_sum(i)/N_sum(i)
         energ_sum2(i)=energ_sum2(i)/N_sum(i)
         if(bNNt(i).gt.1)then
            cheat=(energ_sum2(i)-energ_sum(i)**2)/(aT_rep(i)**2)
            write(20,5006)i,aT_rep(i),energ_sum(i),cheat,bNNa(i)
     $           ,bNNt(i),bNNa(i)/bNNt(i),eshort12_a(i)/N_sum(i)
         else
            write(20,5006)i,aT_rep(i),energ_sum(i),cheat,bNNa(i)
     $           ,bNNt(i)
         endif
      enddo
 5006 format(I4,f7.2,f8.1,f15.3,2f12.1,f11.6,6f8.1)
      write(20,*)'E_min=',E_min

      write(20,*)
      write(20,*)'ncycle_max=',ncycle*n_run
      write(20,*)'ncycle_real=',mcycle
      write(20,*)
      write(20,*)'hour_max=',hour_max
      write(20,*)'hour_real=',atime
      write(20,*)
      write(20,*)'ending time: ',fdate()
      
      write(*,*)
      write(*,*)'ncycle_max=',ncycle*n_run
      write(*,*)'ncycle_real=',mcycle
      write(*,*)
      write(*,*)'hour_max=',hour_max
      write(*,*)'hour_real=',atime
      write(*,*)
      write(*,*)'ending time: ',fdate()
      
      STOP
      END
ccccccc=======================================================cccccccccc
cc          The main program ended!
ccccccc=======================================================cccccccccc




































cccccccccccccccccc set common used parammeters cccccccccccc
      subroutine set_common
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      character protein*10
      common/lengths/Lch,Lch1,Lch2
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/maxdis2/maxdis2(ndim)
      common/arandom/  aarand,abrand,acrand,adrand
      common/nswap/bNSa(100),bNSt(100)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/sw1/aT_rep(nrep),E_rep(nrep)
      character*6 mname
      character*80 line
      common/movename/mname(100)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/distres/er4,es3c
      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      common/zscore/izscore
      common/excluded/vvv(ndim,ndim)
      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
      common/countres1/t_combCA8,t_distL,s_combCA8,s_distL,N_resc
      common/pair1/eh2,eh1b,eh1c
      common/weight/chuan
      character type*10
      common/bigbond/i_bigbond,teco
      common/ssp/ssp
      common/fractpair/fract_pair1,fract_pair3
      common/pair33/i_pair3,mk_pair3
      common/ichos/ichos
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2
      common/nwmax1/nvecr
      common/stick1/nstick,astick,nrmsd,ermsd
      common/stick2/iq(ndim,nrep)
      common/stick3/ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
      common/stick6/bx00(ndim),by00(ndim),bz00(ndim),armsd_min0
      common/res2/er14,er15,er16,er17
      common/svm2/acc_cut,dist_svm2(15),nk_svm(3),ik_svm(3,15)
      common/svm3/er21,er22,er23
      common/svm4/acc_cutB,acc_cutG
      common/svm5/acc_cut0,aw_cont,aw_cont0,aLo1,aLo2,aLo3,anf,itype(10)
      common/svm6/npot,dist_svm(15),dist_svmU(15),dist_svmU2(15)
      
cccccccccccccc set the random generator cccccccccccccccccccccccc
      nozy=random
      if(nozy.gt.0)nozy=-nozy
      firstrandom=aranzy(nozy)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      read(21,*)
      read(21,*)Lch             !length of chain
      read(21,*)protein
      Lch1=Lch-1
      Lch2=Lch-2
      anvec=nvec-0.00001
      contt=1.5*float(Lch)      !TARGET NUMBER OF CONTACTS   1.5*N
      do i=1,Lch
         maxdis2(i)=20*i*i      !maximum distance of walk in i steps
      enddo
      write(20,*)'Target:',protein
      write(20,*)'Length:',Lch
      ichos=1

      if(m_latt.eq.1)then
         latt1=14
         latt2=25
         nvecr=312
      else
         latt1=12
         latt2=26
         nvecr=416
      endif
      anvec=nvecr-0.00001
      
      do i=1,Lch
         do j=1,Lch
            vvv(i,j)=1          !every pair should be checked
         enddo
      enddo

      armsd_min0=1000

c      if(Lch.ge.300)ncycle=int(ncycle*0.8)
c      if(Lch.ge.400)ncycle=int(ncycle*0.8)
c      if(Lch.ge.500)ncycle=int(ncycle*0.8)
c      if(Lch.ge.600)ncycle=int(ncycle*0.8)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      N_resc=0
      t_comb=0
      t_dist=0
      t_distL=0
      t_combCA=0
      t_combCA8=0
      s_comb=0
      s_dist=0
      s_distL=0
      s_combCA=0
      s_combCA8=0

ccccccccccccccccccccccc Temperature ccccccccccccccccccccccccccc
c     [80,130] is the standard:
      if(Lch.lt.80)then
         atemp1=atemp1*0.97
         atemp2=atemp2*0.91
         if(Lch.lt.55)then
            atemp1=atemp1*0.97
            atemp2=atemp2*0.91
         endif
      endif
      if(Lch.gt.130)then
         atemp1=atemp1*1.05
         atemp2=atemp2*1.1
         if(Lch.gt.165)then
            atemp1=atemp1*1.05
            atemp2=atemp2*1.1
            if(Lch.gt.200)then
               atemp1=atemp1*1.05
               atemp2=atemp2*1.1
            endif
         endif
      endif

ccccccccccccc Number of replicas #################
      if(Lch.gt.165)then        !50
         N_rep=N_rep+10
      endif
      if(Lch.gt.240)then        !60
         N_rep=N_rep+10
      endif
      if(Lch.gt.300)then        !70
         N_rep=N_rep+10
      endif
      if(Lch.gt.400)then      !80
         N_rep=N_rep+10
      endif
      if(N_rep.gt.80)N_rep=80

ccccccccccccccc movement name cccccccccccccccccccccccc
      mname(2)='move2a'
      mname(3)='move3s'
      mname(4)='move3d'
      mname(5)='move4s'
      mname(6)='move4d'
      mname(7)='move8'
      mname(8)='move5s'
      mname(9)='move5d'
      mname(10)='move6'
      mname(11)='move_n'
      mname(12)='move_c'
      mname(13)='move7a'
      mname(14)='move7b'
      mname(15)='move9'
      mname(16)='tran_N'
      mname(17)='tran_M'
      mname(18)='tran_C'
      mname(19)='rot_N' !no
      mname(20)='rot_M'
      mname(21)='rot_C' !no
      mname(22)='trot_N'
      mname(23)='trot_M'
      mname(24)='trot_C'
      mname(25)='defo_N'
      mname(26)='defo_M'
      mname(27)='defo_C'

ccccccccccreset restraints weights according to zscore cccccccccccc
      rewind(24)
      read(24,*)n_thr,type
      if(type.eq.'easy')then    !--------------->easy
         izscore=1
         er1=chuan*3.6          !for dist.dat
         er3=chuan*0.765        !for comb.dat
         er4=chuan*0.45         !for comb.dat of deviation
         er5=chuan*2.7          !for combCA.dat
         eh1c=chuan*1.8         !for par.dat
c         er6=chuan*0.45         !for distL.dat
         er7=chuan*500          !for RMSD
         if(ssp.eq.1)er7=chuan*100 !for RMSD
         i_bigbond=3            !decide fragment base on distance
         fract_pair1=0.4
         fract_pair3=0.3
         
         acc_cut=0.4            !control what predicted_accuracy should be used
         acc_cutB=1.1           !ampliphier of CB
         acc_cutG=1.1           !ampliphier of GB
         aw_cont=0.2            !weight of svmseq
         er21=2.5
         er22=2.5
         er23=2.5
      elseif(type.eq.'medm')then !--------------->medium
         izscore=2
         er1=chuan*4.05         !for dist.dat
         er3=chuan*0.81         !for comb.dat
         er4=chuan*0.405        !for comb.dat of deviation
         er5=chuan*1.08         !for combCA.dat
         eh1c=chuan*1.0         !for par.dat
c         er6=chuan*1.0          !for distL.dat
         er7=chuan*100          !for RMSD
         if(ssp.eq.1)er7=chuan*10 !for RMSD
         i_bigbond=1            !decide fragment base on +-2
         fract_pair1=0.7
         fract_pair3=0.3
         
         acc_cut=0.4
         acc_cutB=1.1
         acc_cutG=1.1
         aw_cont=1.0            !weight of svmseq
         er21=2.5
         er22=2.5
         er23=2.5
      else                      !--------------->hard
         izscore=3
         er1=chuan*2.7          !for dist.dat
         er3=chuan*0.4          !for comb.dat
         er4=chuan*0.27         !for comb.dat of deviation
         er5=chuan*0.4          !for combCA.dat
         eh1c=chuan*1.5         !for par.dat
c         er6=chuan*0.5          !for distL.dat
         er7=chuan*5            !for RMSD
         i_bigbond=3            !decide fragment base on distance
         fract_pair1=0.3
         fract_pair3=0.7
         
         acc_cut=0.375
         acc_cutB=1.1
         acc_cutG=1.1
         aw_cont=1.0            !weight of svmseq
         er21=2.0
         er22=2.0
         er23=2.0
      endif
      er21=1.0
      er22=1.0
      er23=1.0
      
      if(teco.eq.1) i_bigbond=1 !teco=1, template from teco
      if(chuan.le.0.0001)then   !without using restraints
         eh1c=1                 !for par.dat
         fract_pair1=1          !do not use par.dat
      endif
c      er14=er5
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c^^^^^^^^^^^^^^^^ common parameters finished ^^^^^^^^^^^^^^^^^^^^^^^      
      return
      end

ccccccccccccccccccc read sequence ccccccccccccccccccccccc
c     Only seq(i) is useful
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_seq
      implicit integer(i-z)
      parameter(ndim=1999)
                parameter(nvec=416)
      character*3 aa(-1:20), NAME,sequ
      common/seqe/seq(ndim),sec(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/icgg/ icg(ndim), EH6  
      common/aminoacid/sequ(ndim)

      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
c                -1    0     1     2     3     4     5     6
     &     'PRO','MET','ASP','ASN','LEU',
c            7     8     9    10    11
     &     'LYS','GLU','GLN','ARG',
c           12    13    14    15
     &     'HIS','PHE','TYR','TRP','CYX'/
c           16    17    18    19    20

      do 121 i=1,Lch
         read(14,707) k,NAME,SEC(I),tmp
         do j=0,19
            if(NAME.eq.aa(j)) then
               SEQ(i)=j
               icg(i)=0
               sequ(i)=name
               if(NAME.eq.'ASP'.or.NAME.eq.'GLU') icg(i)=-1
               if(NAME.eq.'LYS'.or.NAME.eq.'ARG') icg(i)= 1	
               go to 121
            endif
         enddo
         SEQ(i)=0
         icg(i)=0
         sequ(i)='GLY'
 121  continue
 707  format(i5,3x,a3,2i5)
      close(14)

c^^^^^^^^^^^^^^^^^^^^^ read sequence finished ^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read centrosymmetric potential cccccccccccccccccccc
c     eonekd(A) controls centrosymmetric potential of C_g.
c
      subroutine read_centro
      implicit integer(i-z)
      parameter(nvec=416)
      character*3 NAME
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/lengths/Lch,Lch1,Lch2
      common/hopp/eonehw(0:19)

c     read hydrophobic potential for Sg, positive for hydrophobic residue--->
      data eonekd /-0.4, 1.8, -0.8, 2.5, 4.2, -0.7, 4.5,
     &     -1.6, 1.9,  -3.5, -3.5, 3.8,
     &     -3.9, -3.5, -3.5, -4.5,
     &     -3.2, 2.8, -1.3, -0.9/
c            ^          ^     ^     !contradict with 'centro.comm'
c     read hydrophilic potential for Sg, positive for hydrophilic residue--->
      data eonehw /0.0, -0.5, 0.3, -1.0, -1.5, -0.4, -1.8,
     &     0.0, -1.3, 3.0, 0.2, -1.8,
     &     3.0, 3.0, 0.2, 3.0,
     &     -0.5, -2.5, -2.3, -3.4/

c     expected gyration radius:
      acrit=2.2*exp(0.38*alog(float(Lch)))/0.87 !gyrat-radius~2.2*l^0.38
*     Defination of gyration-radius: acrit=sqrt(<(r-r0)^2>)

c^^^^^^^^^^^^^^^^^ read centersymmetric potential finished ^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read environment cccccccccccccccccccccccccccc
c         G  A  V  L  I  S  T  C  M  P  D  N  E  Q  K  R  H  F  Y  W
c (0,0,0) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c (0,0,1) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c (0,0,2) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c   ...                    ...
c (4,4,3) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c (4,4,4) *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_profile
      implicit integer(i-z)
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1

      do kkk=1,4
      do i=0,19
      do im=0,15
      do ip=0,15
      do ia=0,15
         envir(im,ip,ia,i,kkk)=2.0
      end do
      end do
      end do
      end do
      enddo

c     PROFILE3 potential =envir(#of antiparalel contacts,
c     #of orthogonal, # of parallel, aminoacid's type)
c     ia,im, ip - taken modulo 2, i.e. 0-1 contact, 2-3,...
c     profile3.comm is a score table, i.e. envir(envir_class,A)
c     here environment class is number of contacts on residue A.

      do i=0,19                 !from column
         read(2,*)
         do im=0,8
         do ia=0,8
            read(2,*)(envir(ia,im,ip,i,3),ip=0,8) !this is used
         end do
         read(2,*)
         end do
         read(2,*)
      enddo

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,1),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,2),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c      do i=0,19                 !from column
c         read(2,*)
c         do im=0,8
c         do ia=0,8
c            read(2,*)(envir(ia,im,ip,i,4),ip=0,8)
c         end do
c         read(2,*)
c         end do
c         read(2,*)
c      end do

c^^^^^^^^^^^^^^^^^^^^^^^^ read profile finished ^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read 1-3 potential cccccccccccccccccccc
c     interaction between 1'th CA and 3'th CA
      subroutine read_E13
      implicit integer(i-z)
      parameter(ndim=1999)
                parameter(nvec=416)
      common/short2/codevsum,didevsum,csr(ndim,2)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      dimension csre(0:19,0:19,2)

c     R13 potential - two bins only (helical and expanded)
c     r2<48, E=csre(i,j,1); r2>48, E=csre(i,j,2)
      do i=0,19
         do j=0,19
            read(5,*)
            read(5,*) (csre(i,j,k),k=1,2)
         enddo
      enddo

      do i=1,Lch2
         do k=1,2
            csr(i,k)=2.0*csre(seq(i),seq(i+2),k)
         enddo
      enddo

c^^^^^^^^^^^^^^^^^ read E13 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read 1-4 potential cccccccccccccccccccc
c     interaction between 1'th CA and 4'th CA
c     the aim is to obtain IBIN(r14), asr(i,IBIN)
      subroutine read_E14
      implicit integer(i-z)
      parameter(ndim=1999)
                parameter(nvec=416)
      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      common/forpreparemove4/asrr(0:19,0:19,-12:12)
      DIMENSION asrh(0:19,0:19,-12:12),asre(0:19,0:19,-12:12)
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

ccccccccc read asrr,asrh,asre------------------>
      do i=0,19                 !asrr(Ai,Bi,dis) from 'r14.comm'
         do j=0,19
            read(12,*)
            read(12,*) (asrr(i,j,k),k=-12,-5)
            read(12,*) (asrr(i,j,k),k=-4,3) !without k=4
            read(12,*) (asrr(i,j,k),k=5,12)
            do k=4,1,-1
               asrr(i,j,k)=asrr(i,j,k-1) !without k=0
            enddo
         enddo
      enddo
      do i=0,19                 !asrh(Ai,Bi,dis) from 'r14h.comm'
         do j=0,19
            read(7,*)
            read(7,*) (asrh(i,j,k),k=-12,-5)
            read(7,*) (asrh(i,j,k),k=-4,3)
            read(7,*) (asrh(i,j,k),k=5,12)
            do k=4,1,-1
               asrh(i,j,k)=asrh(i,j,k-1)
            enddo
         enddo
      enddo
      do i=0,19                 !asre(Ai,Bi,dis) from 'r14e.comm'
         do j=0,19
            read(8,*)
            read(8,*) (asre(i,j,k),k=-12,-5)
            read(8,*) (asre(i,j,k),k=-4,3)
            read(8,*) (asre(i,j,k),k=5,12)
            do k=4,1,-1
               asre(i,j,k)=asre(i,j,k-1)
            enddo
         enddo
      enddo
c^^^^^^^^^ read asrr,asrh,asre finished ^^^^^^^^^^^^^^^^^

      do i=1,Lch-3
         do k=-12,12
            asr(i,k)=asrr(seq(i+1),seq(i+2),k) !general
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
               if(sec(i).eq.2) then !helix
c                  asr(i,k)=(asr(i,k)+asrh(seq(i+1),seq(i+2),k))/2.0
                  asr(i,k)=asrh(seq(i+1),seq(i+2),k)
               endif
            endif
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
               if(sec(i).eq.4) then !sheet
c                  asr(i,k)=(asr(i,k)+asre(seq(i+1),seq(i+2),k))/2.0
                  asr(i,k)=asre(seq(i+1),seq(i+2),k)*1.5
               endif
            endif
         enddo
      enddo
c^^^^^^^^^^^^ asr(i,ibin(r14)) finished ^^^^^^^^^^^^^^^^^^^^^
c     r(i,i+3)=k, E=asr(i,k), 12 bins (24 bins when considering chiral)
      do i=1,300
         kk=int((sqrt(float(i))*0.87))+1
         if(kk.gt.12) kk=12
         IBIN(I) = kk           !convert lattice r^2 into real r
         IBIN(-I)=-kk
      ENDDO
      IBIN(0)=IBIN(1)

c^^^^^^^^^^^^^^^^^ read E14 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read 1-5 potential cccccccccccccccccccc
c     interaction between 1'th CA and 5'th CA
c     the aim is to obtain JBIN(r15), bsr(i,JBIN)
      subroutine read_E15
      implicit integer(i-z)
      parameter(ndim=1999)
                parameter(nvec=416)
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      DIMENSION bsrh(0:19,0:19,16)
      dimension bsre(0:19,0:19,16)
      dimension bsrr(0:19,0:19,16)
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c

cccccccc read bsrr,bsrh,bsre ----------------->
      do i=0,19                 !read bsrr from 'r15.dat'
         do j=0,19
            read(9,*)
            read(9,*) (bsrr(i,j,k),k=1,8)
            read(9,*) (bsrr(i,j,k),k=9,16)
         enddo
      enddo
      do i=0,19                 !read bsrh from 'r15h.dat'
         do j=0,19
            read(10,*)
            read(10,*) (bsrh(i,j,k),k=1,8)
            read(10,*) (bsrh(i,j,k),k=9,16)
         enddo
      enddo	
      do i=0,19                 !read bsre from 'r15e.dat'
         do j=0,19
            read(11,*)
            read(11,*) (bsre(i,j,k),k=1,8)
            read(11,*) (bsre(i,j,k),k=9,16)
         enddo
      enddo	

      do i=1,Lch-4
         do k=1,16
            bsr(i,k)=bsrr(seq(i+1),seq(i+3),k)
            if(sec(i+1).eq.2.AND.sec(i+2).eq.2.AND.sec(i+3).eq.2)then
c               bsr(i,k)=(bsr(i,k)+bsrh(seq(i+1),seq(i+3),k))/2.0 !helix
               bsr(i,k)=bsrh(seq(i+1),seq(i+3),k)
            endif
            if(sec(i+1).eq.4.AND.sec(i+2).eq.4.AND.sec(i+3).eq.4)then
c               bsr(i,k)=(bsr(i,k)+bsre(seq(i+1),seq(i+3),k))/2.0 !sheet
               bsr(i,k)=bsre(seq(i+1),seq(i+3),k)*1.5
            endif
         enddo
      enddo
c^^^^^^^^^^^^^^^^^^^^ E_15(Ai,Aj,dis) prepared ^^^^^^^^^^^^^^^^^^^

c     prepare distance bin-------------------->
      do i=0,500
         kk=int((sqrt(float(i))*0.87))+1 !i, lattice-dist; kk, real distance
         if(kk.gt.16) kk=16
         JBIN(I) = kk           !jbin: real distance
      ENDDO

ccccc acops(i,jbin) to enhance the contacts between gragments cccc
      do i=1,Lch-4
         acops(i,1)=(min(bsr(i,1),0.0))/2.0 !acpos<0
         do k=2,15
            acops(i,k)=min(0.0,bsr(i,k-1)+2.0*bsr(i,k)+bsr(i,k+1)) !<0
         enddo
         acops(i,16)=(min(bsr(i,16),0.0))/2.0
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c^^^^^^^^^^^^^^^^^ read E15 finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read contact-pair potential cccccccccccccccccccc
c     potential: app(Ai,Aj), apa(Ai,Aj), apm(Ai,Aj)
c     distance range: [arlp,alp], [arla,ala], [arlm,alm]
c     all the general data in 'quarsi3.comm'.
      subroutine read_quarsi3
      implicit integer(i-z)
      parameter(ndim=1999)
                parameter(nvec=416)
      character*3 NAME
      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)
      common/pair1/eh2,eh1b,eh1c
      common/pairmap2/apabla(0:19,0:19),apablp(0:19,0:19)
      common/pairmap3/apablm(0:19,0:19),amap(ndim,ndim),nmap
      common/paircut/ash
      common/zscore/izscore
      common/fractpair/fract_pair1,fract_pair3
      common/pair33/i_pair3,mk_pair3
      common/weight/chuan

      dimension apba(ndim,ndim),apbp(ndim,ndim),apbm(ndim,ndim)

c     Pairwise interactions apablp ... and cut-off parmeters
c     arlp, orientation dependent, pairwise specific, sequence
c     independent

c     read contact-pair potential from 'quarsi3.comm' ------->
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (apablp(i,j),j=0,19) !for app
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (apablm(i,j),j=0,19) !for apm
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (apabla(i,j),j=0,19) !for apa
      enddo
c     read distance-range from 'quarsi3.comm' ------->
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (arlp(i,j),j=0,19) !max distance for parallel
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (arlm(i,j),j=0,19) !for perpendicular contact
      enddo
      read(3,*)
      read(3,*)
      do i=0,19
         read(3,725) NAME, (arla(i,j),j=0,19) !for antiparellel pair
      enddo
 725  format(a3,1x,20f5.1)

c     width of energy-well: [arla,ala]
      ash=0.17
      ash_min=1-ash
      ash_max=1+ash
      do i=0,19
         do j=0,19
            ala(i,j)=(arla(i,j)*ash_max/0.87)**2
            alm(i,j)=(arlm(i,j)*ash_max/0.87)**2
            alp(i,j)=(arlp(i,j)*ash_max/0.87)**2
            arla(i,j)=(arla(i,j)*ash_min/0.87)**2
            arlm(i,j)=(arlm(i,j)*ash_min/0.87)**2
            arlp(i,j)=(arlp(i,j)*ash_min/0.87)**2
         enddo
      enddo
c     E=EH1/2, for r in [0,arlp];
c     E=app-es*fs,  for [0,alp];
c     E=0,     for r in [alp,00].
c^^^^^^^^^^^^^contact interaction range finished ^^^^^^^^^^^^^^^^^

*>>>>>>>>>>>>>
c     read from 'pair3.dat'-------------->
      i_pair3=-1                !without pair3.dat exist
      read(25,*,end=1000)
      rewind(25)
      i_pair3=1
      Nline=100                 !maximum Lch is 2500 !!!
      do i=1,Lch
         read(25,*)
         do i_line=1,Nline
            line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
            read(25,*)(apba(i,j),j=(i_line-1)*25+1,line_end)
            if(line_end.eq.Lch) go to 4074
         enddo
 4074    continue
         read(25,*)
      enddo
      do i=1,Lch
         read(25,*)
         do i_line=1,Nline
            line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
            read(25,*)(apbm(i,j),j=(i_line-1)*25+1,line_end)
            if(line_end.eq.Lch) go to 4075
         enddo
 4075    continue
         read(25,*)
      enddo
      do i=1,Lch
         read(25,*)
         do i_line=1,Nline
            line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
            read(25,*)(apbp(i,j),j=(i_line-1)*25+1,line_end)
            if(line_end.eq.Lch) go to 4076
         enddo
 4076    continue
         read(25,*)
      enddo
      close(25)
 1000 continue
*<<<<<<<<<<<<<<<<

c     combine the data in 'quarsi3.comm' and 'pair3.dat' to get
c     contact potential-------------------->
      do i=1,Lch
         ii=SEQ(i)
         do j=1,Lch
            jj=SEQ(j)
            if(iabs(i-j).lt.5) then
               dd=0.0
            else
               dd=0.25          !encourage contact of distant residues.
            endif

            apa(i,j)=apabla(ii,jj)-dd !quasi3.comm only
            apm(i,j)=apablm(ii,jj)-dd !quasi3.comm only
            app(i,j)=apablp(ii,jj)-dd !quasi3.comm only
            if(i_pair3.gt.0)then
               apa(i,j)=fract_pair3*apba(i,j)+(1-fract_pair3)* !apba: pair3; apabla:quasi3
     &              apabla(ii,jj)-dd
               apm(i,j)=fract_pair3*apbm(i,j)+(1-fract_pair3)*
     &              apablm(ii,jj)-dd
               app(i,j)=fract_pair3*apbp(i,j)+(1-fract_pair3)*
     &              apablp(ii,jj)-dd
            endif
         enddo
      enddo

c^^^^^^^^^^^^^^^ pair-potential is obtained ^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc read contact-pair potential cccccccccccccccccccc
c     the sequence-dependent contact-pair data in 'par.dat'.
      subroutine read_par
      implicit integer(i-z)
      parameter(ndim=1999)
      common/lengths/Lch,Lch1,Lch2
      common/par/apar(ndim,ndim)
      common/fractpair/fract_pair1,fract_pair3
      common/pair33/i_pair3,mk_pair3
      common/weight/chuan
      common/pair1/eh2,eh1b,eh1c
      common/concuttu4/npot4,aw4,Cr20

      dimension apar1(ndim,ndim) !from par.dat
      dimension apar2(ndim,ndim) !from pair1.dat

c     read from 'par.dat'-------------->
      Nline=1000
      do i=1,Lch
         read(15,*)
         do i_line=1,Nline
            line_end=min(10*i_line,Lch) !25,50,75,100,...., ending point
            read(15,*)(apar1(i,j),j=(i_line-1)*10+1,line_end)
            if(line_end.ge.Lch) go to 1
         enddo
 1       continue
      enddo

c     read from 'pair1.dat'-------------->
      if(i_pair3.gt.0)then      ! 'pair1.dat' exist
         Nline=100              !maximum Lch is 2500 !!!
         do i=1,Lch
            read(26,*)
            do i_line=1,Nline
               line_end=min(25*i_line,Lch) !25,50,75,100,...., ending point
               read(26,*)(apar2(i,j),j=(i_line-1)*25+1,line_end)
               if(line_end.eq.Lch) go to 2
            enddo
 2          continue
            read(26,*)
         enddo
      endif
      
***   combine and rescale the pair interaction ---------------->
      apar1_max=-100
      apar2_max=-100
      do i=1,Lch
         do j=1,Lch
            if(abs(apar1(i,j)).gt.apar1_max)apar1_max=abs(apar1(i,j))
            if(abs(apar2(i,j)).gt.apar2_max)apar2_max=abs(apar2(i,j))
         enddo
      enddo
      do i=1,Lch
         do j=1,Lch
            apar(i,j)=0
            if(i_pair3.gt.0)then !'pair3.dat' and 'pair1.dat' exist
               if(apar1_max.lt.0.00001)then !par.dat is wrong
                  apar(i,j)=apar2(i,j) !'pair1.dat' only
               else
                  apar(i,j)=apar1(i,j)*apar2_max/apar1_max*
     &                 (1-fract_pair1)+apar2(i,j)*fract_pair1
               endif
            else
               apar(i,j)=apar1(i,j) !'pair1.dat' not exist, 'par.dat' only
            endif
            apar(i,j)=apar(i,j)*aw4
         enddo
      enddo
      
c      write(*,*)apar1_max,apar2_max
c      do i=1,Lch
c         do j=1,Lch
c            write(*,*)i,j,apar1(i,j),apar2(i,j),apar(i,j)
c         enddo
c      enddo
c      stop

c^^^^^^^^^^^^^^^^^ pair-wise potential finished ^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read cut-off for contact predictions
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_concut
      implicit integer(i-z)
      character*3 NAME
      dimension cut(0:19,0:19),cut_dev(0:19,0:19)
      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
      common/zscore/izscore
      common/concuttu/npot1,concutU2(0:19,0:19),concutU(0:19,0:19),aw1
      common/concuttu1/dq1a(0:19,0:19),dq1b(0:19,0:19)
      common/concuttu2/dq1c(0:19,0:19),dq1d(0:19,0:19)

      rewind(1)
      read(1,*)
      do i=0,19
         read(1,*)NAME,(cut(i,j),j=0,19)
      enddo
      read(1,*)
      do i=0,19
         read(1,*)NAME,(cut_dev(i,j),j=0,19)
      enddo

      do i=0,19
         do j=0,19
            if(izscore.eq.1)then !easy
               dist=7
            elseif(izscore.eq.2)then !medium
               dist=7.5
            else                !hard
               dist=cut(i,j)+cut_dev(i,j)*2.5 !real cut-off
            endif
c            dist=8
            concut(i,j)=dist/0.87 !cut-off on lattice
            concut2(i,j)=concut(i,j)**2 !cut-off squared on lattice
            concutU(i,j)=(dist+2)/0.87 !upper cut-off on lattice
            concutU2(i,j)=concutU(i,j)**2 !upper cut-off squared
            
            dq1a(i,j)=(concutU(i,j)+concut(i,j))/2 !for QUARK potential
            dq1b(i,j)=concutU(i,j)-concut(i,j)
            dq1c(i,j)=(80/0.87+concutU(i,j))/2
            dq1d(i,j)=80/0.87-concutU(i,j)
         enddo
      enddo
      
      return
      end

cccccccccccccccc read contact restrains from 'comb.dat' cccccccccc
c     aim: Mcom(i)--->number of contacts on i-residue
c          Kcom(i,ii): ii's contact of i-residue is j=kres
c     contact restrain is only on C_g.
      subroutine read_contactrestrain
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/lengths/Lch,Lch1,Lch2
      common/distres/er4,es3c
      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      common/resnumber/Ncom,Ndis,accur
      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b
      common/freg/aweig(ndim,ndim)
      common/zscore/izscore
      common/concuttu/npot1,concutU2(0:19,0:19),concutU(0:19,0:19),aw1
      
      DIMENSION r1(8000),r2(8000)

c     READS Side group - side group contacts 
c     (from NMR or therading predictions or clusters)
      
      if(izscore.eq.1)then      !easy target
         cut_min=0.2
         cut0=0.6
      elseif(izscore.eq.2)then  !medium target
         cut_min=0.1
         cut0=0.5
      else                      !hard target
         cut_min=0.1
         cut0=0.4
      endif
ccc   pool all the restraints into r1(i),r2(i)--------------->
 11   continue
      read(16,*)ntmp
      i_c=0
      do i=1,ntmp
         read(16,*)i1,i2,conf
         if(conf.gt.1)conf=1
         if(conf.gt.cut_min)then
            i_c=i_c+1
            if(i_c.gt.8000)then
               write(*,*)'Too many COMB restraints!!!!!!!!'
               cut_min=cut_min*1.1
               rewind(16)
               goto 11
            endif
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweig(i1,i2)=(1+abs(conf-cut0)*4)*aw1
            else
               aweig(i1,i2)=(1-abs(conf-cut0)*2)*aw1
            endif
c            write(*,*)i1,i2,aw1,aweig(i1,i2)
            aweig(i2,i1)=aweig(i1,i2)
c     write(*,*)i1,i2,conf,aweig(i1,i2)
         endif
      enddo
      Ncom=i_c

ccc   map r1,2(i) into Mcom(i),Kcom(i,Mcom(i))------------>
      do i=1,Lch
         Mcom(i)=0              !number of contacts with 'i'
         do j=1,Ncom
            if(r1(j).eq.i)then
               Mcom(i)=Mcom(i)+1
               Kcom(i,Mcom(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               Mcom(i)=Mcom(i)+1
               Kcom(i,Mcom(i))=r1(j)
            endif
         enddo
      enddo

      colim=1.5*Ncom           !background number for derviation
c     the larger 'colim' is, the weaker the contact restrain is.
      
ccc   output restraints------------->
      write(20,*)'Number of restraints:',Ncom
      write(20,*)'----------- contact restraints ---------------'
      nnc=0
      do i=1,Lch
         nnc=nnc+Mcom(i)
         write(20,12)i,Mcom(i),(Kcom(i,j),j=1,Mcom(i))
 12      format(i4,'(',i2,'):',20i4)
      enddo
      write(20,*)'Number of contact=',nnc,' Lch=',Lch
      write(20,*)'fc=',float(nnc)/Lch
      write(20,*)

      if(aw1.le.0.000001)then   !donot waste time for 'comb.dat'
         do i=1,Lch
            Mcom(i)=0           !number of contacts with 'i'
         enddo
      endif

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc read distance restrains from 'dist.dat' cccccccccc
      subroutine read_distantrestrain
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/lengths/Lch,Lch1,Lch2
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/distres/er4,es3c
      common/resnumber/Ncom,Ndis,accur
      common/cutoff/cut1a,cut2a,cut3a,cut4a,cut1b,cut2b,cut3b,cut4b

      dimension r1(10000),r2(10000),dis(10000),deviation(10000)
      
      read(17,*)Ndis
      if(Ndis.gt.10000)then
         write(*,*)'warning: too many short range distance restraints!'
         Ndis=10000
      endif
      do i=1,Ndis
         read(17,*)r1(i),r2(i),nothing,dis(i),deviation(i)
         dis(i)=dis(i)/0.87
         deviation(i)=deviation(i)/0.87
         if(deviation(i).lt.0.5)deviation(i)=0.5
      enddo

      do i=1,Lch
         Mdis(i)=0              !number of prediction for 'i'
         do j=1,Ndis
            if(r1(j).eq.i)then
               Mdis(i)=Mdis(i)+1
               kdis(i,Mdis(i))=r2(j) !r2(j) with 'i'
               dist(i,Mdis(i))=dis(j) !predicted distance for i<->r2(j)
               dev(i,Mdis(i))=deviation(j)
            endif
            if(r2(j).eq.i)then
               Mdis(i)=Mdis(i)+1
               kdis(i,Mdis(i))=r1(j)
               dist(i,Mdis(i))=dis(j) !predicted distance for i<->r2(j)
               dev(i,Mdis(i))=deviation(j)
            endif
         enddo
      enddo
c      dilim=1.5*Ndis          !background number for derviation
      dilim=0.5*Ndis          !background number for derviation

c      write(*,*)
c      write(*,*)'----------- distant map ---------------'
c      nnc=0
c      do i=1,Lch
c         nnc=nnc+Mdis(i)
c         write(*,12)i,Mdis(i),(Kdis(i,j),dist(i,j)*0.87,j=1,Mdis(i))
c 12      format(i4,':',i3,20(i4,'-'f5.2))
c      enddo
c      write(*,*)'Number of distmap=',nnc,' Lch=',Lch

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc read distance restrains from 'dist.dat' cccccccccc
c   note: kdisL(ndim,i), i is the order number not the residue
c   each pair (i,j) can have many different distances predictions
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_longdistantrestrain
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      COMMON/longdist/MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
c     maximum N=4*L/10, if L=100, maximum N=400
      dimension r1(400000),r2(400000),dis(400000)

      read(22,*)NdisL
      if(NdisL.gt.400000)NdisL=400000
      do i=1,NdisL
         read(22,*)r1(i),r2(i),dis(i)
         dis(i)=dis(i)/0.87
      enddo
      
      do i=1,Lch
         MdisL(i)=0             !number of prediction for 'i'
         do j=1,NdisL
            if(r1(j).eq.i)then
               MdisL(i)=MdisL(i)+1
               kdisL(i,MdisL(i))=r2(j) !r2(j) with 'i'
               distL(i,MdisL(i))=dis(j) !predicted distance for i<->r2(j)
               if(MdisL(i).ge.500)goto 101
            endif
            if(r2(j).eq.i)then
               MdisL(i)=MdisL(i)+1
               kdisL(i,MdisL(i))=r1(j)
               distL(i,MdisL(i))=dis(j) !predicted distance for i<->r2(j)
               if(MdisL(i).ge.500)goto 101
            endif
         enddo
 101     continue
         if(MdisL(i).gt.500)then !size>5000/4
            write(*,*)i,MdisL(i),'N_res>500 for one residue, exit'
            stop
         endif
      enddo

c      write(*,*)
c      write(*,*)'----------- distant map ---------------'
c      nnc=0
c      do i=1,Lch
c        nnc=nnc+MdisL(i)
c        write(*,12)i,MdisL(i),(KdisL(i,j),distL(i,j)*0.87,j=1,MdisL(i))
c 12     format(i4,':',i3,20(i4,'-'f5.2))
c      enddo
c      write(*,*)'Number of distmap=',nnc,' Lch=',Lch
c      stop
      
c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc read distance restrains from 'dist.dat' cccccccccc
c   note: kdisL(ndim,i), i is the order number not the residue
c   each pair (i,j) can have many different distances predictions
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_exp
      implicit integer(i-z)
      parameter(ndim=1999)
      common/lengths/Lch,Lch1,Lch2
      common/expose/mp(20,ndim),area(ndim)
      common/nana1/nana

      na=nana
      read(27,*)
      do i=1,Lch
         read(27,*)itmp,(mp(j,i),j=1,na)
         area(i)=0
         do j=1,na
            if(mp(j,i).eq.0)mp(j,i)=-1 !-1, bury; 1, expose
            area(i)=area(i)+mp(j,i)
         enddo
c         write(*,*)i,mp(1,i),mp(2,i),mp(12,i)
      enddo
      close(27)
c      stop

c^^^^^^^^^^^^^^^^^^ read contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc read contact restrains from 'combCA.dat' cccccccccc
c     aim: Mcom(i)--->number of contacts on i-residue
c          Kcom(i,ii): ii's contact of i-residue is j=kres
c     CAcontact restrain is only on CA.
      subroutine read_combCA
      implicit integer(i-z)
      parameter(ndim=1999)
      common/lengths/Lch,Lch1,Lch2
      common/zscore/izscore
      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(ndim,ndim)
      common/CAc1/npot2,d_CA_cut,d_CA_cut2,d_CA_cutU,d_CA_cutU2,aw2
      common/CAcontact1a/dq2a,dq2b,dq2c,dq2d
      DIMENSION r1(8000),r2(8000)
      
      d_tmp=6.5                 !
      d_CA_cut=d_tmp/0.87
      d_CA_cut2=d_CA_cut**2
      d_CA_cutU=(d_tmp+2)/0.87
      d_CA_cutU2=d_CA_cutU**2
      
      dq2a=(d_CA_cutU+d_CA_cut)/2 !for QUARK potential
      dq2b=d_CA_cutU-d_CA_cut
      dq2c=(80/0.87+d_CA_cutU)/2
      dq2d=80/0.87-d_CA_cutU
      
      if(izscore.eq.1)then      !easy target
         cut_min=0.2
         cut0=0.5
      elseif(izscore.eq.2)then  !medium target
         cut_min=0.1
         cut0=0.4
      else                      !hard target
         cut_min=0.1
         cut0=0.3
      endif
ccc   pool all the restraints into r1(i),r2(i)--------------->
 11   continue
      read(18,*)ntmp
      i_c=0
      do i=1,ntmp
         read(18,*)i1,i2,conf
c         conf=float(i3)/float(i4)
         if(conf.gt.1)conf=1
         if(conf.gt.cut_min)then
c            write(*,*)'===----',i1,i2,conf,cut_min
            i_c=i_c+1
            if(i_c.gt.8000)then
               write(*,*)'Too many COMBCA restraints!!!!!'
               cut_min=cut_min*1.1
               rewind(18)
               goto 11
            endif
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweigCA(i1,i2)=(1+abs(conf-cut0)*4)*aw2
            else
               aweigCA(i1,i2)=(1-abs(conf-cut0)*2)*aw2
            endif
            aweigCA(i2,i1)=aweigCA(i1,i2)
c            write(*,*)i1,i2,conf,aweigCA(i1,i2)
         endif
      enddo
      NcomCA=i_c

ccc   map r1,2(i) into McomCA(i),KcomCA(i,McomCA(i))------------>
      do i=1,Lch
         McomCA(i)=0              !number of contacts with 'i'
         do j=1,NcomCA
            if(r1(j).eq.i)then
               McomCA(i)=McomCA(i)+1
               KcomCA(i,McomCA(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               McomCA(i)=McomCA(i)+1
               KcomCA(i,McomCA(i))=r1(j)
            endif
         enddo
      enddo

ccc   output restraints------------->
c      write(*,*)'Number of restraints:',NcomCA
c      write(*,*)'----------- CAcontact restraints ---------------'
c      nnc=0
c      do i=1,Lch
c         nnc=nnc+McomCA(i)
c         write(*,12)i,McomCA(i),(KcomCA(i,j),j=1,McomCA(i))
c 12      format(i4,'(',i2,'):',20i4)
c         write(*,13)i,McomCA(i),(aweigCA(i,KcomCA(i,j)),j=1,McomCA(i))
c 13      format(i4,'(',i2,'):',20f8.5)
c      enddo
c      write(*,*)'Number of CAcontact=',nnc,' Lch=',Lch
c      write(*,*)
c      stop
      
      if(aw2.le.0.000001)then   !donot waste time for 'combCA.dat'
         do i=1,Lch
            McomCA(i)=0         !number of contacts with 'i'
         enddo
      endif
      
c^^^^^^^^^^^^^^^^^^ read CAcontact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc read contact restrains from 'combCA.dat' cccccccccc
c     aim: McomCB6(i)--->number of contacts on i-residue
c          KcomCB6(i,ii): ii's contact of i-residue is j=kres
c     CAcontact restrain is only on CB.
      subroutine read_comb6CB
      implicit integer(i-z)
      parameter(ndim=1999)
      common/lengths/Lch,Lch1,Lch2
      common/zscore/izscore
      
      common/CB6a/npot5,aw5,f_CB_cut,f_CB_cut2,f_CB_cutU,f_CB_cutU2
      common/CB6b/fq2a,fq2b,fq2c,fq2d
      common/CB6c/McomCB6(ndim),KcomCB6(ndim,100),aweigCB6(ndim,ndim)
      
      DIMENSION r1(8000),r2(8000)
      
      f_tmp=6.0                 !
      f_CB_cut=f_tmp/0.87
      f_CB_cut2=f_CB_cut**2
      f_CB_cutU=(f_tmp+2)/0.87
      f_CB_cutU2=f_CB_cutU**2
      
      fq2a=(f_CB_cutU+f_CB_cut)/2 !for QUARK potential
      fq2b=f_CB_cutU-f_CB_cut
      fq2c=(80/0.87+f_CB_cutU)/2
      fq2d=80/0.87-f_CB_cutU
      
      if(izscore.eq.1)then      !easy target
         cut_min=0.2
         cut0=0.5
      elseif(izscore.eq.2)then  !medium target
         cut_min=0.1
         cut0=0.4
      else                      !hard target
         cut_min=0.1
         cut0=0.3
      endif

ccc   pool all the restraints into r1(i),r2(i)--------------->
 11   continue
      read(29,*)ntmp
      i_c=0
      do i=1,ntmp
         read(29,*)i1,i2,conf
         if(conf.gt.1)conf=1
         if(conf.gt.cut_min)then
            i_c=i_c+1
            if(i_c.gt.8000)then
               write(*,*)'Too many comb6CB restraints!!!!!'
               cut_min=cut_min*1.1
               rewind(29)
               goto 11
            endif
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweigCB6(i1,i2)=(1+abs(conf-cut0)*4)*aw5
            else
               aweigCB6(i1,i2)=(1-abs(conf-cut0)*2)*aw5
            endif
            aweigCB6(i2,i1)=aweigCB6(i1,i2)
         endif
      enddo
      NcomCB6=i_c
      
ccc   map r1,2(i) into McomCB6(i),KcomCB6(i,McomCB6(i))------------>
      do i=1,Lch
         McomCB6(i)=0           !number of contacts with 'i'
         do j=1,NcomCB6
            if(r1(j).eq.i)then
               McomCB6(i)=McomCB6(i)+1
               KcomCB6(i,McomCB6(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               McomCB6(i)=McomCB6(i)+1
               KcomCB6(i,McomCB6(i))=r1(j)
            endif
         enddo
      enddo
      
c      write(*,*)'aw5=',aw5
      if(aw5.le.0.000001)then   !donot waste time for 'combCB6.dat'
         do i=1,Lch
            McomCB6(i)=0        !number of contacts with 'i'
         enddo
      endif
      
      goto 32
cc   output restraints------------->
      write(*,*)'Number of restraints:',NcomCB6
      write(*,*)'----------- CB6contact restraints ---------------'
      nnc=0
      do i=1,Lch
         nnc=nnc+McomCB6(i)
         write(*,12)i,McomCB6(i),(KcomCB6(i,j),j=1,McomCB6(i))
 12      format(i4,'(',i2,'):',20i4)
         write(*,13)i,McomCB6(i),
     &        (aweigCB6(i,KcomCB6(i,j)),j=1,McomCB6(i))
 13      format(i4,'(',i2,'):',20f8.5)
      enddo
      write(*,*)'Number of CB6contact=',nnc,' Lch=',Lch
      write(*,*)
 32   continue
      
c      stop
      
c^^^^^^^^^^^^^^^^^^ read CB6contact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc read contact restrains from 'comb.dat' cccccccccc
c     aim: Mcom(i)--->number of contacts on i-residue
c          Kcom(i,ii): ii's contact of i-residue is j=kres
c     CAcontact restrain is only on CA.
      subroutine read_comb8CA
      implicit integer(i-z)
      parameter(ndim=1999)
      common/lengths/Lch,Lch1,Lch2
      common/zscore/izscore
      common/CA8/McomCA8(ndim),KcomCA8(ndim,100),aweigCA8(ndim,ndim)
      common/CA8a/npot3,e_CA_cut,e_CA_cut2,e_CA_cutU,e_CA_cutU2,aw3
      common/CA8b/eq2a,eq2b,eq2c,eq2d
      common/CA8c/dcut3
      
      DIMENSION r1(8000),r2(8000)
      
c      e_tmp=8.0                 !cutoff for comb8CA.dat
      e_tmp=dcut3               !cutoff for comb8CA.dat
      e_CA_cut=e_tmp/0.87
      e_CA_cut2=e_CA_cut**2
      e_CA_cutU=(e_tmp+2)/0.87
      e_CA_cutU2=e_CA_cutU**2
      
      eq2a=(e_CA_cutU+e_CA_cut)/2 !for QUARK potential
      eq2b=e_CA_cutU-e_CA_cut
      eq2c=(80/0.87+e_CA_cutU)/2
      eq2d=80/0.87-e_CA_cutU
      
      if(izscore.eq.1)then      !easy target
         cut_min=0.2
         cut0=0.5
      elseif(izscore.eq.2)then  !medium target
         cut_min=0.1
         cut0=0.4
      else                      !hard target
         cut_min=0.1
         cut0=0.3
      endif
      
ccc   pool all the restraints into r1(i),r2(i)--------------->
 11   continue
      read(28,*)ntmp
      i_c=0
      do i=1,ntmp
         read(28,*)i1,i2,conf
         if(conf.gt.1)conf=1
         if(conf.gt.cut_min)then
            i_c=i_c+1
            if(i_c.gt.8000)then
               write(*,*)'Too many COMB8CA restraints!!!!!'
               cut_min=cut_min*1.1
               rewind(28)
               goto 11
            endif
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweigCA8(i1,i2)=(1+abs(conf-cut0)*4)*aw3
            else
               aweigCA8(i1,i2)=(1-abs(conf-cut0)*2)*aw3
            endif
            aweigCA8(i2,i1)=aweigCA8(i1,i2)
         endif
      enddo
      NcomCA=i_c

ccc   map r1,2(i) into McomCA(i),KcomCA(i,McomCA(i))------------>
      do i=1,Lch
         McomCA8(i)=0           !number of contacts with 'i'
         do j=1,NcomCA
            if(r1(j).eq.i)then
               McomCA8(i)=McomCA8(i)+1
               KcomCA8(i,McomCA8(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               McomCA8(i)=McomCA8(i)+1
               KcomCA8(i,McomCA8(i))=r1(j)
            endif
         enddo
      enddo
      
ccc   output restraints------------->
c      write(*,*)'Number of restraints:',NcomCA
c      write(*,*)'----------- CAcontact restraints ---------------'
c      nnc=0
c      do i=1,Lch
c         nnc=nnc+McomCA(i)
c         write(*,12)i,McomCA(i),(KcomCA(i,j),j=1,McomCA(i))
c 12      format(i4,'(',i2,'):',20i4)
c         write(*,13)i,McomCA(i),(aweigCA(i,KcomCA(i,j)),j=1,McomCA(i))
c 13      format(i4,'(',i2,'):',20f8.5)
c      enddo
c      write(*,*)'Number of CAcontact=',nnc,' Lch=',Lch
c      write(*,*)
c      stop

      if(aw3.le.0.000001)then   !donot waste time for 'combCA8.dat'
         do i=1,Lch
            McomCA8(i)=0        !number of contacts with 'i'
         enddo
      endif
      
c^^^^^^^^^^^^^^^^^^ read CAcontact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc read contact restrains from 'comb.dat' cccccccccc
c     aim: McomCB8(i)--->number of contacts on i-residue
c          KcomCB8(i,ii): ii's contact of i-residue is j=kres
c     CAcontact restrain is only on CA.
      subroutine read_comb8CB
      implicit integer(i-z)
      parameter(ndim=1999)
      common/lengths/Lch,Lch1,Lch2
      common/zscore/izscore

      common/CB8a/npot6,aw6,g_CB_cut,g_CB_cut2,g_CB_cutU,g_CB_cutU2
      common/CB8b/gq2a,gq2b,gq2c,gq2d
      common/CB8c/McomCB8(ndim),KcomCB8(ndim,100),aweigCB8(ndim,ndim)
      
      DIMENSION r1(8000),r2(8000)
      
      g_tmp=8.0                 !cutoff for comb8CB.dat
      g_CB_cut=g_tmp/0.87
      g_CB_cut2=g_CB_cut**2
      g_CB_cutU=(g_tmp+2)/0.87
      g_CB_cutU2=g_CB_cutU**2
      
      gq2a=(g_CB_cutU+g_CB_cut)/2 !for QUARK potential
      gq2b=g_CB_cutU-g_CB_cut
      gq2c=(80/0.87+g_CB_cutU)/2
      gq2d=80/0.87-g_CB_cutU
      
      if(izscore.eq.1)then      !easy target
         cut_min=0.2
         cut0=0.5
      elseif(izscore.eq.2)then  !medium target
         cut_min=0.1
         cut0=0.4
      else                      !hard target
         cut_min=0.1
         cut0=0.3
      endif
      
ccc   pool all the restraints into r1(i),r2(i)--------------->
 11   continue
      read(30,*)ntmp
      i_c=0
      do i=1,ntmp
         read(30,*)i1,i2,conf
         if(conf.gt.1)conf=1
         if(conf.gt.cut_min)then
c            write(*,*)'----',i1,i2,conf,cut_min
            i_c=i_c+1
            if(i_c.gt.8000)then
               write(*,*)'Too many comb8CB restraints!!!!!'
               cut_min=cut_min*1.1
               rewind(30)
               goto 11
            endif
            r1(i_c)=i1
            r2(i_c)=i2
            if(conf.gt.cut0)then
               aweigCB8(i1,i2)=(1+abs(conf-cut0)*4)*aw6
            else
               aweigCB8(i1,i2)=(1-abs(conf-cut0)*2)*aw6
            endif
            aweigCB8(i2,i1)=aweigCB8(i1,i2)
         endif
      enddo
      NcomCB8=i_c

ccc   map r1,2(i) into McomCB8(i),KcomCB8(i,McomCB8(i))------------>
      do i=1,Lch
         McomCB8(i)=0           !number of contacts with 'i'
         do j=1,NcomCB8
            if(r1(j).eq.i)then
               McomCB8(i)=McomCB8(i)+1
               KcomCB8(i,McomCB8(i))=r2(j)
            endif
            if(r2(j).eq.i)then
               McomCB8(i)=McomCB8(i)+1
               KcomCB8(i,McomCB8(i))=r1(j)
            endif
         enddo
      enddo
      
      if(aw6.le.0.000001)then   !donot waste time for 'combCB8.dat'
         do i=1,Lch
            McomCB8(i)=0        !number of contacts with 'i'
         enddo
      endif
      
      goto 32
cc   output restraints------------->
      write(*,*)'Number of restraints:',NcomCB8
      write(*,*)'----------- CB8contact restraints ---------------'
      nnc=0
      do i=1,Lch
         nnc=nnc+McomCB8(i)
         write(*,12)i,McomCB8(i),(KcomCB8(i,j),j=1,McomCB8(i))
 12      format(i4,'(',i2,'):',20i4)
         write(*,13)i,McomCB8(i),
     &        (aweigCB8(i,KcomCB8(i,j)),j=1,McomCB8(i))
 13      format(i4,'(',i2,'):',20f8.5)
      enddo
      write(*,*)'Number of CB8contact=',nnc,' Lch=',Lch
      write(*,*)
 32   continue
      
c      stop
      
c^^^^^^^^^^^^^^^^^^ read comb8CB restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc read contact restrains by from 'svmseqca8.dat' cccccccccc
c     aim: Mcom(i)--->number of contacts on i-residue
c          Kcom(i,ii): ii's contact of i-residue is j=kres
c     CAcontact restrain is only on CA.
      subroutine read_seqcontact
      implicit integer(i-z)
      parameter(ndim=1999)
      common/lengths/Lch,Lch1,Lch2
      common/svm5/acc_cut0,aw_cont,aw_cont0,aLo1,aLo2,aLo3,anf,itype(10)
      
      common/svm1/Mcon(3,ndim),Kcon(3,ndim,300),awei(3,ndim,ndim),fw
      common/svm2/acc_cut,dist_svm2(15),nk_svm(3),ik_svm(3,15)
      common/svm4/acc_cutB,acc_cutG
      common/svm6/npot,dist_svm(15),dist_svmU(15),dist_svmU2(15)
      common/dwell/dwell,ha(15),hb(15),hc(15),hd(15)
      
      DIMENSION r1(8000),r2(8000),conf0(4,15),idist(15)
      dimension i_atom_type(8000)
      
ccc initialize contact, energy will neglect contact when 'contact.map' not exist
      do m=1,3                  !m=1: CA; m=2: CB; m=3: SG
         do i=1,Lch
            Mcon(m,i)=0         !number of m-th contacts with 'i'
         enddo
         
         dist=8                 !suppose all distance cutoff is 8A
         dist_svm(m)=dist/0.87
         dist_svm2(m)=(dist/0.87)**2
         dist_svmU(m)=(dist+dwell)/0.87
         dist_svmU2(m)=((dist+dwell)/0.87)**2
         
         ha(m)=(dist_svmU(m)+dist_svm(m))/2 !for QUARK potential
         hb(m)=dist_svmU(m)-dist_svm(m)
         hc(m)=(80/0.87+dist_svmU(m))/2
         hd(m)=80/0.87-dist_svmU(m)
      enddo
      
ccc   read contact predictions ------->
      nc_tot=0                  !total number of contacts
      read(50,*,end=11)         !read head of contact.map
 10   read(50,*,end=11)i1,i2,conf,ia
      nc_tot=nc_tot+1
      r1(nc_tot)=i1
      r2(nc_tot)=i2
      awei(ia,i1,i2)=conf       !weight factor of contact (i,j)
      awei(ia,i2,i1)=conf
      i_atom_type(nc_tot)=ia
      if(nc_tot.lt.8000)then
         goto 10                !maximum 8000 contacts are allowed
      endif
 11   continue
      
ccc   map r1,2(i) into Mcon(k,i),Kcon(k,i,Mcon(i))------------>
      n_max_per_residue=300
      do 72 m=1,2               !only CA,CB
         nc_tot_m=0
         do i=1,Lch
            Mcon(m,i)=0         !number of contacts with 'i'
            do j=1,nc_tot
               if(i_atom_type(j).eq.m)then
                  if(r1(j).eq.i)then
                     if(Mcon(m,i).lt.n_max_per_residue)then !maximum numb of contacts/residue <300
                        nc_tot_m=nc_tot_m+1
                        Mcon(m,i)=Mcon(m,i)+1 !number of contact on residue-i
                        Kcon(m,i,Mcon(m,i))=r2(j) !what residue contact with residue-i
                     else
                        write(*,*)m,i,j,Mcon(m,i),'warn: Mcon(m,i)>300'
                        write(20,*)m,i,j,Mcon(m,i),'warn: Mcon(m,i)>300'
                     endif
                  endif
                  if(r2(j).eq.i)then
                     if(Mcon(m,i).lt.n_max_per_residue)then !maximum numb of contacts/residue <300
                        nc_tot_m=nc_tot_m+1
                        Mcon(m,i)=Mcon(m,i)+1
                        Kcon(m,i,Mcon(m,i))=r1(j) !each contact counted twice
                     else
                        write(*,*)m,i,j,Mcon(m,i),'warn: Mcon(m,i)>300'
                        write(20,*)m,i,j,Mcon(m,i),'warn: Mcon(m,i)>300'
                     endif
                  endif
               endif
            enddo
         enddo
         
         goto 130
c        output restraints------------->
         if(m.eq.1) write(*,*)'Number of CA-contacts:',nc_tot_m
         if(m.eq.2) write(*,*)'Number of CB-contacts:',nc_tot_m
         if(m.eq.1) write(20,*)'Number of CA-contacts:',nc_tot_m
         if(m.eq.2) write(20,*)'Number of CB-contacts:',nc_tot_m
         nnc=0
         do i=1,Lch
            nnc=nnc+Mcon(m,i)
            write(*,12)i,Mcon(m,i),(Kcon(m,i,j),j=1,Mcon(m,i))
            write(*,13)i,Mcon(m,i),
     &           (awei(m,i,Kcon(m,i,j)),j=1,Mcon(m,i))
            write(20,12)i,Mcon(m,i),(Kcon(m,i,j),j=1,Mcon(m,i))
            write(20,13)i,Mcon(m,i),
     &           (awei(m,i,Kcon(m,i,j)),j=1,Mcon(m,i))
 12         format(i4,'(',i3,'):',300i5)
 13         format(i4,'(',i3,'):',300f5.1)
         enddo
 130     continue
 72   continue
      
c      stop

c^^^^^^^^^^^^^^^^^^ read CAcontact restrains finished ^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccset reset temperature ccccccccccccccccccc
      subroutine reset_temperature
      implicit integer(i-z)
      parameter(ndim=1999)      !number of residues
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)       !number of vectors
      common/lengths/Lch,Lch1,Lch2
      common/resnumber/Ncom,Ndis,accur
      common/commonuse2/atemp1,atemp2,N_rep,phot
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/distres/er4,es3c
      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      real r_dis,r_con,r_dev,T10,T20,T1a,T2a

      if(er1+er3+er4.lt.0.1)return !without restrains
      
      a_rest=float(Ncom)/(1.3*Lch)
      a_rest=sqrt(a_rest)
      if(a_rest.lt.0.875)a_rest=0.875 !because of 80/70 -> rest/without_rest
      atemp2=atemp2*a_rest
      if(atemp2.gt.115) atemp2=115.0

      return
      end

cccccccccccccccccccset temperature for different replicas ccccc
      subroutine set_temperature
      implicit integer(i-z)
      parameter(nrep=100)
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/thrtem/aT_ann(100)
      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs

***********for normal run ***********************************
      do i=1,N_REP
         aT_rep(i)=atemp1*(atemp2/atemp1)**(float(i-1)/(N_rep-1))
         aTs_rep(i)=aTs1*(aTs2/aTs1)**(float(i-1)/(N_rep-1))
         aTTs_rep(i)=aTTs1*(aTTs2/aTTs1)**(float(i-1)/(N_rep-1))
c     write(*,*)i,aT_rep(i)
      enddo

c      stop
c^^^^^^^^^^^^^^^^^ set aT_rep(i) finished ^^^^^^^^^^^^^^^
      return
      end

cccccccccccccc set EHBIJ(i,j) ccccccccccccccccccc
      subroutine set_EHB
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim)

c
c     EHBIJ - set-up secondary structure dependent
c     strength of the hyrogen bond network - stronger for helices
c     and beta-sheets
c

      do i=1,Lch
         is=sec(i)
         do j=1,Lch
            js=sec(j)
            EHBIJ(i,j)=1
            if(iabs(i-j).eq.3.and.is.eq.2.and.js.eq.2)then
               EHBIJ(i,j)=EHBIJ(i,j)+0.5 !helix, enhanced
            endif
            if(is.eq.4.or.js.eq.4) then
               if(is*js.ne.8.and.iabs(i-j).gt.4)then
                  EHBIJ(i,j)=EHBIJ(i,j)+0.5 !beta-beta, enhanced
               endif
            endif
         enddo
      enddo

c^^^^^^^^^^^ set H_bond finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     prepare possible vector v=(vx,vy,vz) satisfy |v*v|=14-25 (3.26A-4.35A)
c     and vector(-5:5,-5:5,-5:5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_vectors
      implicit integer(i-z)
      parameter(nvec=416)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension n(100)
      common/lattice/m_latt,latt1,latt2

      do i=1,100
         n(i)=0
      enddo

      nwmax=0
      aaa=0
      nn=5
      do x=-nn,nn
         do y=-nn,nn
            do z=-nn,nn
               vector(x,y,z)=0
               r=x*x+y*y+z*z
               if(r.ge.latt1.and.r.le.latt2) then
                  nwmax=nwmax+1
                  vx(nwmax)=x
                  vy(nwmax)=y
                  vz(nwmax)=z
                  vector(x,y,z)=nwmax
c                  write(*,*)nwmax,vx(nwmax),vy(nwmax),vz(nwmax),r,
c     $                 sqrt(float(r)*0.87*0.87)
                  n(r)=n(r)+1
                  aaa=aaa+sqrt(float(r)*0.87*0.87)
               endif
            enddo
         enddo
      enddo
      aaa=aaa/float(nwmax)
      write(20,*)'n1_all=',nwmax,'  <vr>=',aaa

c      do i=1,nwmax
c         write(*,*)i,vx(i),vy(i),vz(i)
c      enddo
c      stop

c      do i=10,30
c         write(*,*)i,n(i),sqrt(float(i)*0.87*0.87)
c      enddo
c      stop
      
c     i=1,5
c           10          24   2.751182    
c           11          24   2.885463    
c           12           8   3.013768    +
c           13          24   3.136830    +
c           14          48   3.255242    x
c           15           0   3.369496    
c           16           6   3.480000    x
c           17          48   3.587102    x
c           18          36   3.691097    x
c           19          24   3.792242    x
c           20          24   3.890758    x
c           21          48   3.986841    x
c           22          24   4.080662    x
c           23           0   4.172373    
c           24          24   4.262112    x
c           25          30   4.350000    x
c           26          72   4.436147    +
c           27          32   4.520653    
c           28           0   4.603607    
c           29          72   4.685093    
c           30          48   4.765186    
c nwmax=         616  <vr>=   3.982909    
c nwmax=         312  <vr>=   3.809868    

c      stop
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c define goodc(i,j), angle(i,j), prod(i,j)
c     prod(i,j): v(i)*v(j).
c     angle(i,j): angle of neighbor bonds, i,j--->(1:nvec)
c     goodc(i,j): ture, when angle(i,j) in [60,160]; false, otherwise.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_neighbors
      implicit integer(i-z)
                parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/three/angle(nvec,nvec)
      common/lengths/Lch,Lch1,Lch2

      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)
      common/nwmax1/nvecr

      max_m12=0
      do i=1,nvecr
         m12(i)=0
      enddo

      mmm=0
      nnn=0
      kkk=0
      do i=1,nvecr
      do j=1,nvecr
         u21(i,j)=0
         a2=vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)
         b2=vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)
         c2=(vx(i)+vx(j))**2+(vy(i)+vy(j))**2+(vz(i)+vz(j))**2
         cosangle=(a2+b2-c2)/(2*sqrt(a2*b2))
         angle(i,j)=acos(cosangle)*180/3.1415926
c     in database, angle is in [65,165];
         if(angle(i,j).gt.65.and.angle(i,j).lt.165)then
            goodc(i,j)=.true.
            mmm=mmm+1
            ijx=vx(i)+vx(j)
            ijy=vy(i)+vy(j)
            ijz=vz(i)+vz(j)
            do k=1,nvecr
               if(vx(k).eq.ijx.and.vy(k).eq.ijy.and.vz(k).eq.ijz)then
                  kkk=kkk+1
                  u21(i,j)=k    !vi+vj=vk
                  m12(k)=m12(k)+1
                  u1(k,m12(k))=i
                  u2(k,m12(k))=j
                  if(max_m12.lt.m12(k))max_m12=m12(k)
                  goto 10
               endif
            enddo
 10         continue
         else
            goodc(i,j)=.false.
         endif
         nnn=nnn+1
c         write(*,*)i,j,mmm,nnn,angle(i,j),goodc(i,j)
      enddo
      enddo

      n=0
      do i=1,nvecr
         r=vx(i)**2+vy(i)**2+vz(i)**2
         if(m12(i).gt.0)n=n+1
c     if(r.gt.17)write(*,*)i,r,m12(i)
      enddo
c      write(*,*)'n2_all=',nnn
c      write(*,*)'n2good=',mmm,'  n21=',kkk
c      write(*,*)'n1_all=',nvecr,'  n12=',n

c      stop
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c define hydrogen-bond, C_beta, C_group for all possible neighbors
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_beta
      implicit integer(i-z)
                parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/three/angle(nvec,nvec)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/cb/hx(nvec,nvec,0:19),hy(nvec,nvec,0:19),hz(nvec,nvec,0:19)
      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)
      common/beta3/bxalf(0:19),byalf(0:19),bzalf(0:19)
      common/beta4/bxbet(0:19),bybet(0:19),bzbet(0:19)
      common/nwmax1/nvecr


      do k=0,19
         read(4,*)
         read(4,*)axalf(k),ayalf(k),azalf(k),axbet(k),aybet(k),azbet(k)
      enddo
      do k=0,19
         read(4,*)
         read(4,*)bxalf(k),byalf(k),bzalf(k),bxbet(k),bybet(k),bzbet(k)
      enddo
      CLOSE(4)                  !sidecent.comm

ccccccccccc define hydrogen-bond, C_beta, C_group for good (i,j)----->
      esp=0.000001
      do 101 i=1,nvecr
      do 102 j=1,nvecr
         avi=sqrt(float(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i)))
         avxi=vx(i)/avi
         avyi=vy(i)/avi
         avzi=vz(i)/avi
         avj=sqrt(float(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j)))
         avxj=vx(j)/avj
         avyj=vy(j)/avj
         avzj=vz(j)/avj
         
ccc   if vi and vj is parallel, a is ok but b=c=0 ------->
         if(abs(avxi-avxj).lt.esp)then
            if(abs(avyi-avyj).lt.esp)then
               if(abs(avzi-avzj).lt.esp)then
                  ax=avxi+avxj
                  ay=avyi+avyj
                  az=avzi+avzj
                  aaa=sqrt(ax*ax+ay*ay+az*az)
                  ax=ax/aaa     !(vi+vj)/|vi+vj|
                  ay=ay/aaa     !(vi+vj)/|vi+vj|
                  az=az/aaa     !(vi+vj)/|vi+vj|
                  
c     calculate b=a(x)u, u=(1,1,1):
                  bx=ay*1-az*1
                  by=az*1-ax*1
                  bz=ax*1-ay*1
                  bbb=sqrt(bx*bx+by*by+bz*bz)
                  bx=bx/bbb     ! a(x)1/|a(x)1|
                  by=by/bbb     ! a(x)1/|a(x)1|
                  bz=bz/bbb     ! a(x)1/|a(x)1|
                  
c     calculate c=a(x)b:
                  cx=ay*bz-az*by
                  cy=az*bx-ax*bz
                  cz=ax*by-ay*bx
                  ccc=sqrt(cx*cx+cy*cy+cz*cz)
                  cx=cx/ccc     ! a(x)b/|a(x)b|
                  cy=cy/ccc     ! a(x)b/|a(x)b|
                  cz=cz/ccc     ! a(x)b/|a(x)b|
                  cax(i,j)=cx
                  cay(i,j)=cy
                  caz(i,j)=cz
                  goto 29
               endif
            endif
         endif
ccc   check if vi and vj is anti-parallel, c is ok, b=a=0 ------->
         if(abs(avxi+avxj).lt.esp)then
            if(abs(avyi+avyj).lt.esp)then
               if(abs(avzi+avzj).lt.esp)then
                  cx=avxi-avxj
                  cy=avyi-avyj
                  cz=avzi-avzj
                  ccc=sqrt(cx*cx+cy*cy+cz*cz)
                  cx=cx/ccc     !(vi-vj)/|vi-vj|
                  cy=cy/ccc     !(vi-vj)/|vi-vj|
                  cz=cz/ccc     !(vi-vj)/|vi-vj|
                  cax(i,j)=cx   !(vi-vj)/|vi-vj|
                  cay(i,j)=cy   !(vi-vj)/|vi-vj|
                  caz(i,j)=cz   !(vi-vj)/|vi-vj|

c     calculate a=c(x)u, u=(1,1,1):
                  ax=cy*1-cz*1
                  ay=cz*1-cx*1
                  az=cx*1-cy*1
                  aaa=sqrt(ax*ax+ay*ay+az*az)
                  ax=ax/aaa     ! c(x)1/|c(x)1|
                  ay=ay/aaa     ! c(x)1/|c(x)1|
                  az=az/aaa     ! c(x)1/|c(x)1|
                  
c     calculate b=c(x)a:
                  bx=cy*az-cz*ay
                  by=cz*ax-cx*az
                  bz=cx*ay-cy*ax
                  bbb=sqrt(bx*bx+by*by+bz*bz)
                  bx=bx/bbb     ! c(x)a/|c(x)a|
                  by=by/bbb     ! c(x)a/|c(x)a|
                  bz=bz/bbb     ! c(x)a/|c(x)a|
                  goto 29
               endif
            endif
         endif
         
         ax=avxi+avxj
         ay=avyi+avyj
         az=avzi+avzj
         aaa=sqrt(ax*ax+ay*ay+az*az)
         ax=ax/aaa              !(vi+vj)/|vi+vj|
         ay=ay/aaa              !(vi+vj)/|vi+vj|
         az=az/aaa              !(vi+vj)/|vi+vj|
         
         bx=avyi*avzj-avzi*avyj
         by=avzi*avxj-avxi*avzj
         bz=avxi*avyj-avyi*avxj
         bbb=sqrt(bx*bx+by*by+bz*bz)
         bx=bx/bbb             ! vi(x)vj/|vi(x)vj|
         by=by/bbb             ! vi(x)vj/|vi(x)vj|
         bz=bz/bbb             ! vi(x)vj/|vi(x)vj|

         cx=avxi-avxj
         cy=avyi-avyj
         cz=avzi-avzj
         ccc=sqrt(cx*cx+cy*cy+cz*cz)
         cx=cx/ccc              !(vi-vj)/|vi-vj|
         cy=cy/ccc              !(vi-vj)/|vi-vj|
         cz=cz/ccc              !(vi-vj)/|vi-vj|
         cax(i,j)=cx            !(vi-vj)/|vi-vj|
         cay(i,j)=cy            !(vi-vj)/|vi-vj|
         caz(i,j)=cz            !(vi-vj)/|vi-vj|
         
 29      continue

         goto 39
c     check if aaa/bbb/ccc=0 ------->
         if(aaa.lt.esp.or.bbb.lt.esp.or.ccc.lt.esp)then
            write(*,19)i,j,avxi,avyi,avzi,aaa,bbb,ccc
         endif
c     check if a.b=b.c=a.c=0 ------->
         ab=ax*bx+ay*by+az*bz
         bc=bx*cx+by*cy+bz*cz
         ac=ax*cx+ay*cy+az*cz
         if(ab.gt.esp.or.ab.gt.esp.or.ac.gt.esp)then
            write(*,19)i,j,avxi,avyi,avzi,ab,bc,ac
         endif
c         write(*,19)i,j,avxi,avyi,avzi,ab,bc,ac
 19      format(2i5,10f8.3)
 39       continue

c     H-bond (unit vector):
         hbx(i,j)=bx
         hby(i,j)=by
         hbz(i,j)=bz

c     side-chain coordinate from C_a to side-chain ---------------->
         do k=0,19
            if(angle(i,j).lt.105) then ! alpha-helix or turn like
               gx(i,j,k)=(axalf(k)*ax+ayalf(k)*bx+azalf(k)*cx)/0.87
               gy(i,j,k)=(axalf(k)*ay+ayalf(k)*by+azalf(k)*cy)/0.87
               gz(i,j,k)=(axalf(k)*az+ayalf(k)*bz+azalf(k)*cz)/0.87

               hx(i,j,k)=(bxalf(k)*ax+byalf(k)*bx+bzalf(k)*cx)/0.87
               hy(i,j,k)=(bxalf(k)*ay+byalf(k)*by+bzalf(k)*cy)/0.87
               hz(i,j,k)=(bxalf(k)*az+byalf(k)*bz+bzalf(k)*cz)/0.87
            else                ! beta-sheet
               gx(i,j,k)=(axbet(k)*ax+aybet(k)*bx+azbet(k)*cx)/0.87
               gy(i,j,k)=(axbet(k)*ay+aybet(k)*by+azbet(k)*cy)/0.87
               gz(i,j,k)=(axbet(k)*az+aybet(k)*bz+azbet(k)*cz)/0.87

               hx(i,j,k)=(bxbet(k)*ax+bybet(k)*bx+bzbet(k)*cx)/0.87
               hy(i,j,k)=(bxbet(k)*ay+bybet(k)*by+bzbet(k)*cy)/0.87
               hz(i,j,k)=(bxbet(k)*az+bybet(k)*bz+bzbet(k)*cz)/0.87
            endif
         enddo
 102  continue
 101  continue

c      stop
      return
      end

cccccccccc Compute the secondary fragment biases cccccccccccccccc
c     check local secondary structure from 'seq.dat'
c     if it is beta-structure in [i,i+6],  frga(i)=19.1/0.87;
c     if it is alpha-structure in [i,i+7], frgb(i)=10.5/0.87;
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_frg
      implicit integer(i-z)
      parameter(ndim=1999)
                parameter(nvec=416)
      common/fr/frga(ndim),frgb(ndim)
      common/seqe/seq(ndim),sec(ndim)
      common/lengths/Lch,Lch1,Lch2

      do i=1,Lch
         frga(i)=0.0
         frgb(i)=0.0
      enddo

      do i=1,Lch-7
         q=0
         do j=i,i+7
            if(sec(j).eq.2) q=q+2 !helix structure.
         enddo
         if(q.eq.16)then        !8 continue alpha-residues
            frga(i)=10.5/0.87   !distance for 7 alpha-bonds
         endif
      enddo

      do i=1,Lch-6
         q=0
         do j=i+1,i+5
            if(sec(j).eq.4) q=q+4 !beta structure
         enddo
         if(q.eq.20)then        !5 continue beta-residues
            if(sec(i).ne.2.and.sec(i+6).ne.2)then
               frgb(i)=19.1/0.87 !distance for 6 beta-bonds
            endif
         endif
      enddo

c      do i=1,Lch
c         write(*,*)i,seq(i),sec(i),frga(i),frgb(i)
c      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     count the number of restraints satisfied
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine count_restrains
      implicit integer (i-z)
      parameter(nrep=100)
      parameter(ndim=1999)	!maximum length of chain-length
      parameter(nvec=416)

      common/lengths/Lch,Lch1,Lch2
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/seqe/seq(ndim),sec(ndim)
      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(ndim,ndim)
      common/CAc1/npot2,d_CA_cut,d_CA_cut2,d_CA_cutU,d_CA_cutU2,aw2
      COMMON/longdist/MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      common/CA8/McomCA8(ndim),KcomCA8(ndim,100),aweigCA8(ndim,ndim)
      common/CA8a/npot3,e_CA_cut,e_CA_cut2,e_CA_cutU,e_CA_cutU2,aw3
      common/CA8b/eq2a,eq2b,eq2c,eq2d
 
      dimension ax(ndim),ay(ndim),az(ndim) !SG
      dimension ica(0:ndim),x(ndim),y(ndim),z(ndim) !CA
      
      common/countres/t_comb,t_dist,t_combCA,s_comb,s_dist,s_combCA
      common/countres1/t_combCA8,t_distL,s_combCA8,s_distL,N_resc
      
      do i=1,3
         N_resc=N_resc+1
         do j=1,Lch
            x(j)=xrep(j,i)
            y(j)=yrep(j,i)
            z(j)=zrep(j,i)
         enddo
         do j=1,Lch1
            wx=x(j+1)-x(j)
            wy=y(j+1)-y(j)
            wz=z(j+1)-z(j)
            ica(j)=vector(wx,wy,wz) !identify order number of each bond-vector
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         do j=1,Lch
            ax(j)=x(j)+gx(ica(j-1),ica(j),seq(j))
            ay(j)=y(j)+gy(ica(j-1),ica(j),seq(j))
            az(j)=z(j)+gz(ica(j-1),ica(j),seq(j))
         enddo

         do j=1,Lch
c     comb.dat--------->
            do k=1,Mcom(j)
               m=Kcom(j,k)      !k'th contact with j
               if(m.gt.j)then
                  t_comb=t_comb+1
                  dis=(ax(j)-ax(m))**2+(ay(j)-ay(m))**2
     $                 +(az(j)-az(m))**2
                  if(dis.le.concut2(seq(j),seq(m)))then
                     s_comb=s_comb+1
                  endif
               endif
            enddo
c     dist.dat--------->
            do k=1,Mdis(j)
               m=kdis(j,k)      !j-m restraints
               if(m.gt.j) then  !to avoid repeat
                  t_dist=t_dist+1
                  dis2=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  dis=sqrt(dis2)
                  err=abs(dis-dist(j,k)) !dist: predicted dis
                  if(err.lt.dev(j,k)) then !dev: deviation for arca
                     s_dist=s_dist+1
                  endif
               endif
            enddo
c     distL.dat--------->
            do k=1,MdisL(j)
               m=kdisL(j,k)     !j-m restraints
               if(m.gt.j) then  !to avoid repeat
                  t_distL=t_distL+1
                  dis2=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  dis=sqrt(dis2)
                  err=abs(dis-distL(j,k)) !dist: predicted dis
                  if(err.lt.1.5) then !dev: deviation for arca
                     s_distL=s_distL+1
                  endif
               endif
            enddo
c     combCA.dat--------->
            do k=1,McomCA(j)
               m=KcomCA(j,k)      !k'th contact with j
               if(m.gt.j)then
                  t_combCA=t_combCA+1
                  dis=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  if(dis.le.d_CA_cut2)then
                     s_combCA=s_combCA+1
                  endif
               endif
            enddo
c     combCA8.dat--------->
            do k=1,McomCA8(j)
               m=KcomCA8(j,k)   !k'th contact with j
               if(m.gt.j)then
                  t_combCA8=t_combCA8+1
                  dis=(x(j)-x(m))**2+(y(j)-y(m))**2+(z(j)-z(m))**2
                  if(dis.le.e_CA_cut2)then
                     s_combCA8=s_combCA8+1
                  endif
               endif
            enddo
ccc   
         enddo
      enddo

c^^^^^^^^^^^^^^^^^ restraints count finished ^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prepare 2-bond move, i.e. calculate v21(tx,ty,tz,i), v22(tx,ty,tz,i).
c v21(tx,ty,tz,i) --- 1th vector of i'th path from (0,0,0) to (tx,ty,tz)
c v22(tx,ty,tz,i) --- 2th vector of i'th path from (0,0,0) to (tx,ty,tz)
c it will be dangerous if number of used variable is more than 3,000,000
c i.e. the usage of memory of CUP can not beyond 90%.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_move2
      implicit integer(i-z)
      parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)
      common/nwmax1/nvecr
      common/lattice/m_latt,latt1,latt2

      do i=-10,10
         do j=-10,10
            do k=-10,10
               Nw(i,j,k)=0
            enddo
         enddo
      enddo

      max=0
      maxa=0
      nnn=0
      mmm=0

      do 101 i=1,nvecr
         do 102 j=1,nvecr
            if(goodc(i,j))then
               ijx=vx(i)+vx(j)  !vx in (-5,5)
               ijy=vy(i)+vy(j)
               ijz=vz(i)+vz(j)
***   based on ijx:
               Nw(ijx,ijy,ijz)=Nw(ijx,ijy,ijz)+1 !based on r
               w21(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=i
               w22(ijx,ijy,ijz,Nw(ijx,ijy,ijz))=j
               if(maxa.lt.Nw(ijx,ijy,ijz))maxa=Nw(ijx,ijy,ijz)
***   based on i,j:
               Np2(i,j)=0
               do 103 ii=1,nvecr
                  jx=ijx-vx(ii)
                  jy=ijy-vy(ii)
                  jz=ijz-vz(ii)
                  jr=jx*jx+jy*jy+jz*jz
                  if(jr.ge.latt1.and.jr.le.latt2)then
                     jj=vector(jx,jy,jz)
                     if(goodc(ii,jj))then
                        Np2(i,j)=Np2(i,j)+1 !based on i,j
                        v21(i,j,Np2(i,j))=ii
                        v22(i,j,Np2(i,j))=jj
                        mmm=mmm+1 !number of total memory occupied by vx21.
                     endif
                  endif
 103           continue
               if(max.lt.Np2(i,j))max=Np2(i,j)
c               write(*,*)i,j,Np2(i,j),max,maxa
               nnn=nnn+1        !number of possible pairs
            endif
 102     continue
 101  continue

ccc among all 312*312=97344 pairs, 67272 pairs are legal (~70%).
ccc all Np2(i,j)>=2, i.e. there are at least one other path for any pair.
ccc <Np2(i,j)>=27.
      
      write(20,*)'maximum of Np2(i,j)=',max,maxa
      write(20,*)'number of possible pair=',nnn
      write(20,*)'sum of Np2(i,j), total memory=',mmm
cccc  the following is the cases without goodc on prefabracated pairs:
c              N_v  N_pare N_pp   mmm
c     [14,25], 312, 97344, 67272, 1830000 (25% memory, within memory)

cccc the following is cases without goodc limitation:
c     [14,25], 306, 93636, 3872214 (90% memory, beyond memory)
c     [14,24], 282, 79524, 2923758 (80% memory, beyond memory)
c      stop

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c prepare 3-bond move, i.e. calculate v31(tx,ty,tz,i), v32(tx,ty,tz,i).
c v31(tx,ty,tz,i) --- 1th vector of i'th path from (0,0,0) to (tx,ty,tz)
c v32(tx,ty,tz,i) --- 2th vector of i'th path from (0,0,0) to (tx,ty,tz)
c v33(tx,ty,tz,i) --- 3th vector of i'th path from (0,0,0) to (tx,ty,tz)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prepare_move3
      implicit integer(i-z)
      parameter(nvec=416)
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/move31/v31(-15:15,-15:15,-15:15,6)
      common/move32/v32(-15:15,-15:15,-15:15,6)
      common/move33/v33(-15:15,-15:15,-15:15,6)
      common/move34/Np3(-15:15,-15:15,-15:15)
      common/nwmax1/nvecr

c      dimension v31(-15:15,-15:15,-15:15,3150)
c      dimension v32(-15:15,-15:15,-15:15,3150)
c      dimension v33(-15:15,-15:15,-15:15,3150)
c      dimension Np3(-15:15,-15:15,-15:15)

      do i=-15,15
         do j=-15,15
            do k=-15,15
               Np3(i,j,k)=0
            enddo
         enddo
      enddo


      num3=0
      max=0
      do 101 i=1,nvecr
         do 102 j=1,nvecr
            if(goodc(i,j))then
               do 103 k=1,nvecr
                  if(goodc(j,k))then
                     num3=num3+1
                     rx=vx(i)+vx(j)+vx(k)
                     ry=vy(i)+vy(j)+vy(k)
                     rz=vz(i)+vz(j)+vz(k)
                     Np3(rx,ry,rz)=Np3(rx,ry,rz)+1
                     v31(rx,ry,rz,Np3(rx,ry,rz))=i
                     v32(rx,ry,rz,Np3(rx,ry,rz))=j
                     v33(rx,ry,rz,Np3(rx,ry,rz))=k
                     if(max.le.Np3(rx,ry,rz)) max=Np3(rx,ry,rz)
                     write(*,*)i,j,k,max
                  endif
 103           continue
            endif
 102     continue
 101  continue

      write(*,*)'number of move3=',num3,max !14507376        3140
      stop
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c calculate contact order according to the length and secondary structure
c	  SIZE     H       E       H/E
c  	2   40  0.116   0.324      0.252
c   	3   60  0.119   0.357      0.230
c  	4   80  0.115   0.280      0.212
c   	5  100  0.105   0.259      0.198
c   	6  120  0.132   0.269      0.168
c   	7  140  0.105   0.272      0.176
c   	8  160  0.114   0.186      0.183
c   	9  180  0.116   0.197      0.160
c      10  200  0.107   0.184      0.134
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_acorder
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      common/seqe/seq(ndim),sec(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/lengths/Lch,Lch1,Lch2
      
c     data struct /' coil','helix',' turn',' beta','    ?'/

************* number of predicted structures *************
      n_H=0                     !number of Helix
      n_E=0                     !number of Extension
      do i=1,Lch
         if(sec(i).eq.2) n_H=n_H+1
         if(sec(i).eq.4) n_E=n_E+1
      enddo

      if(n_H+n_E.lt.2)then      !use alpha=H/E
         if(Lch.lt.50)then
            alph=0.252
         else if(Lch.lt.70)then
            alph=0.230
         else if(Lch.lt.90)then
            alph=0.212
         else if(Lch.lt.110)then
            alph=0.198
         else if(Lch.lt.130)then
            alph=0.168
         else if(Lch.lt.150)then
            alph=0.176
         else if(Lch.lt.170)then
            alph=0.183
         else if(Lch.lt.190)then
            alph=0.160
         else
            alph=0.134
         endif
      else                      !use alpha=aH+bE
         a1=float(n_H)/float(n_H+n_E)
         a2=float(n_E)/float(n_H+n_E)
         if(Lch.lt.50)then
            alph=0.116*a1+0.324*a2
         else if(Lch.lt.70)then
            alph=0.119*a1+0.357*a2
         else if(Lch.lt.90)then
            alph=0.115*a1+0.280*a2
         else if(Lch.lt.110)then
            alph=0.105*a1+0.259*a2
         else if(Lch.lt.130)then
            alph=0.132*a1+0.269*a2
         else if(Lch.lt.150)then
            alph=0.105*a1+0.272*a2
         else if(Lch.lt.170)then
            alph=0.114*a1+0.186*a2
         else if(Lch.lt.190)then
            alph=0.116*a1+0.197*a2
         else
            alph=0.107*a1+0.184*a2
         endif
      endif

      acorder=alph*Lch
c^^^^^^^^^^^^^^^^^^ contact order done ^^^^^^^^^^^^^^^^^^
      return
      end

*****************************************************************
*     read initial (x,y,z) from init.dat
*****************************************************************
      subroutine read_initial
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/initialinput/switch,k_cycle,k_phot,N_ann
      character c1*4,c2*2,c3*3,aaaaa*10,text,text1*22
      character*80 head
      character sign(ndim)
      common/initialstr/cx(ndim),cy(ndim),cz(ndim)
      common/seqe/seq(ndim),sec(ndim)
      common/looks/exc,exc1,exc2

      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      dimension ax(ndim),ay(ndim),az(ndim)
      dimension x0(ndim),y0(ndim),z0(ndim)
      dimension M_i(ndim),M_f(ndim)

      common/chain0/ras(ndim),nfl
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6

      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/echain2/egx(ndim),egy(ndim),egz(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C

      common/chainm/mv(ndim)
      dimension q_bk(ndim)

      character*3 sequ
      common/aminoacid/sequ(ndim)
      common/bigbond/i_bigbond,teco
      common/rmsdder/i_chunk(ndim),ex0(ndim),ey0(ndim),ez0(ndim)

      common/consensus1/cx0(0:ndim),cy0(0:ndim),cz0(0:ndim),q(0:ndim)
      dimension n_i(ndim),n_f(ndim)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      common/stick1/nstick,astick,nrmsd,ermsd
      common/stick2/iq(ndim,nrep)
      common/stick3/ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
      common/nwmax1/nvecr

      i_rep=0                   !replicas

 104  rewind(24)
      read(24,*)n_thr           !Number of real templates in 'init.dat'
      do 102 i_init=1,n_thr

****************************************************
c     read q(i), cx0(i) from 'init.dat'
****************************************************
         do i=1,Lch
            q(i)=0
         enddo
         read(24,*)L_ali
         do i=1,L_ali
            read(24,1237)text,ii,text,a1,a2,a3
            cx0(ii)=a1
            cy0(ii)=a2
            cz0(ii)=a3
            q(ii)=1
         enddo
         read(24,*)text         !TER
 1237    format(A22,I4,A4,3F8.3)
c^^^^^^^^^^q(i) and cx0(i) obtained ^^^^^^^^^^^^^^^^^^^^

********************************************************
c     remove small segments, so that we have a better random walk
********************************************************
         M_a=0                  !number of aligned segments
         q(0)=0
         do i=1,Lch
            if(q(i).ne.q(i-1).and.q(i).eq.1)then
               M_a=M_a+1
               M_i(M_a)=i
            endif
            if(q(i).eq.1)M_f(M_a)=i
         enddo
         do i=1,M_a
            L_a=M_f(i)-M_i(i)+1
            if(L_a.le.L_cut)then
               do j=M_i(i),M_f(i)
                  q(j)=0
               enddo
            endif
         enddo
c^^^^^^^^^^^remove small segment finished ^^^^^^^^^^^^^

********************************************************
c     check and find GAP, only according to q(i)
********************************************************
         N_g=0                  !number of gaps
         q(0)=1
         do i=1,Lch
            if(q(i).ne.q(i-1).and.q(i).eq.0)then
               N_g=N_g+1
               n_i(N_g)=i       !initial point of the gap
            endif
            if(q(i).eq.0)n_f(N_g)=i !final point of the gap
         enddo
c^^^^^^^^^^^^^^^check gap finished ^^^^^^^^^^^^^^^^^

********************************************************
c     fill GAP: cx0(i) -> cx(i)
********************************************************
         n_walk=0               !for number of reject by excluded volumn
         exc_eff=exc            !for excluded volumn
 70      n_walk=n_walk+1
         if(n_walk.gt.100)then
            exc_eff=exc_eff*0.99
         endif
         if(n_walk.gt.10000)then
            write(20,*)'unsolvable structure',i
            write(*,*)'unsolvable structure',i
            stop
         endif
         do i=1,Lch
            if(q(i).eq.1)then
               cx(i)=cx0(i)
               cy(i)=cy0(i)
               cz(i)=cz0(i)
            else
               cx(i)=1000000.   !for checking excluded volumn
               cy(i)=1000000.
               cz(i)=1000000.
            endif
         enddo
         do i=2,n_g
            i1=n_i(i)
            i2=n_f(i)
            call connect(i1,i2,pass) !fill missed cooridinates, no move others
            if(pass.ge.3)goto 70 !re-walk
         enddo
         if(n_g.ge.1)then
            i1=n_i(1)
            i2=n_f(1)
            call connect(i1,i2,pass)
            if(pass.ge.3)goto 70 !re-walk
         endif
*^^^^^^^^^^^^^^Fill gap done, cx(i) is continuous ^^^^^^^^^^^^^^^^^^^^^^^

*************************************************************************
c     project chain onto lattices, decide (x,y,z):
*************************************************************************
***   ax(i)=cx(i)/0.87, for the decision of (x,y,z):
         do i=1,Lch
            ax(i)=cx(i)/0.87    !C_alpha scaled by 0.87
            ay(i)=cy(i)/0.87
            az(i)=cz(i)/0.87
         enddo
c     (ax,ay,az) into (x,y,z)----------------------->
         x(1)=nint(ax(1))
         y(1)=nint(ay(1))
         z(1)=nint(az(1))
         xm=x(1)
         ym=y(1)
         zm=z(1)
         do 101 i=2,Lch
            dis_min=100000000000. !minimum distance between c_alpha and lattice
            j_ch=0              !chosen vector
            do 100 j=1,nvecr
               if(i.gt.2)then   !check good neighbor
                  if(.not.goodc(jm,j))goto 100
               endif
               x_tmp=xm+vx(j)
               y_tmp=ym+vy(j)
               z_tmp=zm+vz(j)
c     check excluded volumn---->
               do m=1,i-3       !!!!! on-lattice part.
                  disaa=(x_tmp-x(m))**2+(y_tmp-y(m))**2+(z_tmp-z(m))**2
                  if(disaa.lt.exc_eff) goto 100
               enddo
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^
               dis=(ax(i)-x_tmp)**2+(ay(i)-y_tmp)**2+(az(i)-z_tmp)**2
               if(dis.lt.dis_min)then
                  j_ch=j
                  x_ch=x_tmp
                  y_ch=y_tmp
                  z_ch=z_tmp
                  dis_min=dis
               endif
 100        continue
            if(j_ch.lt.1)goto 70 !refill the gaps
            x(i)=x_ch           !get (x(i),y(i),z(i)) here
            y(i)=y_ch
            z(i)=z_ch
            jm=j_ch
            xm=x(i)
            ym=y(i)
            zm=z(i)
 101     continue
c^^^^^^^^^^^^^^^^^^^project lattice done ^^^^^^^^^^^^^^^^^^^^^^^

***   record the initial conformation of k'th replica--->
         i_rep=i_rep+1
         do i=1,Lch
            xrep(i,i_rep)=x(i)  !ica will be calculated in set_current
            yrep(i,i_rep)=y(i)
            zrep(i,i_rep)=z(i)
            iq(i,i_rep)=q(i)
            if(q(i).eq.1)then
               ax00(i,i_rep)=cx0(i)/0.87
               ay00(i,i_rep)=cy0(i)/0.87
               az00(i,i_rep)=cz0(i)/0.87
            endif
         enddo
         if(i_rep.ge.N_rep)goto 105
 102  continue
      if(i_rep.lt.N_rep)goto 104

******prepare movement for normal movement *****************
 105  nfl=Lch
      do i=1,nfl
         ras(i)=i
      enddo
      call move_point           !decide movement point for notmal run
      nfr=0                     !number of frozen fragments.

c      do i=1,nfl
c         write(*,*)i,ras(i),ras2(i),ras3(i),ras4(i),ras5(i)
c      enddo
c      open(unit=70,file='initial.pdb',status='unknown')
c      write(70,*)N_rep
c      do j=1,N_rep
c         write(70,*)Lch
c         do i=1,Lch
c            write(70,1037)i,sequ(i),i,xrep(i,j)*0.87,yrep(i,j)*0.87,zrep(i,j)*0.87
c         enddo
c         write(70,*)'TER'
c      enddo
c 1037 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
c      close(70)
c      stop

c^^^^^^^^^^^^ read initial chain finished ^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   produce initial structures randomly:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine random_initial
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
                parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      dimension ip(ndim)
      common/sw3/icarep(ndim,nrep)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/chain0/ras(ndim),nfl
      common/looks/exc,exc1,exc2
      common/ranzy/nozy
      common/nwmax1/nvecr

      do 102 k=1,N_REP
 88      x(1)=0
         y(1)=0
         z(1)=0
         m=0
         do 101 i=2,Lch
 99         ip(i)=int(aranzy(nozy)*nvecr)+1
            m=m+1
            if(m.gt.1000000)then
               write(*,*) 'UNSOLVABLE STERIC PROBLEM > EXIT_2'
               goto 88          !unsolvable steric problem
            endif
            if(i.gt.2)then      !check neighbor
               if(.not.goodc(ip(i-1),ip(i)))goto 99
            endif
            x(i)=x(i-1)+vx(ip(i))
            y(i)=y(i-1)+vy(ip(i))
            z(i)=z(i-1)+vz(ip(i))
            do j=1,i-1          !check excluded volumn for Ca
               ir=(x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2
               if(ir.lt.exc) goto 99
            enddo
 101     continue

         do i=1,Lch
            xrep(i,k)=x(i)
            yrep(i,k)=y(i)
            zrep(i,k)=z(i)
         enddo
 102  continue

******prepare movement for normal movement *****************
      nfl=Lch
      do i=1,nfl
         ras(i)=i
      enddo
      call move_point           !decide movement point for notmal run
      nfr=0                     !number of frozen fragments.

c^^^^^^^^^^^^ initial chains finished ^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   displace some of initial structures using threading structures
c     cx0(i): gapped
c     cx(i):  gap-filled
c     ax(i):  ax=cx/0.87
c     x(i):   on-lattice
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     The procedure of this subroutine is following:
c     1, read coordinates of templates, i.e. cx0(i);
c     2, find gaps that no coordinates given in template (n_i,n_f,n_g);
c     3, fill gaps by random walk. when gap too big, keep the last big-bond. 
c        (cx0->cx)
c     4, dicide movable points, +-2 residues around gap, nfl,ras(i)
c     5, find big-bond, add movable points acconding to width;
c     6, add ending points, cx->ax;
c     7, project ax onto lattices, ax->x;
c     8, decide number of bond-movements, nf2,ras2(i).
c     9, decide number of frozen-fragments, nfr,nfr_i(i),nfr_f(i)
cTTTccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine template_initial
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/xyzrs/xrs(ndim,40,40),yrs(ndim,40,40),zrs(ndim,40,40)
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/rs/i_thr0
      common/initialinput/switch,k_cycle,k_phot,N_ann
      character c1*4,c2*2,c3*3,aaaaa*10,text,text1*22
      character*80 head
      character sign(ndim)
      common/initialstr/cx(ndim),cy(ndim),cz(ndim)
      common/seqe/seq(ndim),sec(ndim)
      common/looks/exc,exc1,exc2

      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/distres/er4,es3c
      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)

      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep)
      common/excluded/vvv(ndim,ndim)

      common/consensus1/cx0(0:ndim),cy0(0:ndim),cz0(0:ndim),q(0:ndim)

      common/mng/m_g(100)
      dimension n_i(ndim),n_f(ndim)
      dimension cx0b(ndim),cy0b(ndim),cz0b(ndim)
      dimension cx0g(ndim),cy0g(ndim),cz0g(ndim)
      dimension ax(ndim),ay(ndim),az(ndim)
      dimension x0(ndim),y0(ndim),z0(ndim)

      common/beta1/axalf(0:19),ayalf(0:19),azalf(0:19)
      common/beta2/axbet(0:19),aybet(0:19),azbet(0:19)
      common/beta3/bxalf(0:19),byalf(0:19),bzalf(0:19)
      common/beta4/bxbet(0:19),bybet(0:19),bzbet(0:19)

      dimension M_i(ndim),M_f(ndim)

      dimension mf(0:ndim)

      common/chain0/ras(ndim),nfl
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/echain2/egx(ndim),egy(ndim),egz(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/ssp/ssp

      common/chainm/mv(ndim)
      dimension q_bk(ndim)

      character*3 sequ
      common/aminoacid/sequ(ndim)
      common/bigbond/i_bigbond,teco
      common/rmsdder/i_chunk(ndim),ex0(ndim),ey0(ndim),ez0(ndim)

      common/nwmax1/nvecr

      dimension loop_i(ndim),loop_f(ndim),loop_p(ndim)
      common/mloopf/mloop
      dimension q_old(ndim)

******************************************************
c     set-up the template to use
******************************************************
      exc_eff=exc               !for excluded volumn
      M_consensus=0             !not use consensus
      if(i_thr0.eq.0)then       !means use consensus of top-2 template
         i_thr0=1               !if failed, using the first template
         M_consensus=1          !use consensus in the following
      endif

****************************************************
c     read q(i), cx0(i) from 'init.dat'
****************************************************
      rewind(24)
      read(24,*)n_thr           !Number of real templates in 'init.dat'
      if(i_thr0.gt.n_thr)then
         write(20,*)'without threading structure at this i_thr0!'
         write(*,*)'without threading structure at this i_thr0!'
         stop
      endif
      do k=1,i_thr0
         do i=1,Lch
            q(i)=0
         enddo
         read(24,*)Nal
         do i=1,Nal
            read(24,1237)text,ii,text,a1,a2,a3
            cx0(ii)=a1
            cy0(ii)=a2
            cz0(ii)=a3
            q(ii)=1
         enddo
         read(24,*)text
      enddo
 1237 format(A22,I4,A4,3F8.3)
      if(M_consensus.eq.1)call get_consensus !get q(i),cx0(i) from consensus
c^^^^^^^^^^q(i) and cx0(i) obtained ^^^^^^^^^^^^^^^^^^^^

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     cut 5 loops if the template is full-length
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mloop=5
c     read(19,*) mloop          !how many broken points if full-length
      write(20,*)'mloop=',mloop
      do i=1,Lch
         q_old(i)=q(i)
      enddo
***** check whether the template is full ------------->
      nq0=0                     !number of un-aligned regions
      nbigb=0                   !number of big-bonds
      do i=1,Lch
         if(q(i).eq.0)then
            nq0=nq0+1
         endif
         if(i.ge.2)then
            if(q(i-1).eq.1.and.q(i).eq.1)then
               dis=di(cx0(i-1),cy0(i-1),cz0(i-1),cx0(i),cy0(i),cz0(i))
               if(dis.gt.5.0)then
                  nbigb=nbigb+1
               endif
            endif
         endif
      enddo
      write(20,*)'number of un-aligned regions, nq0=',nq0
      write(20,*)'number of big-bonds, nbigb=',nbigb
      if(nq0.eq.0.and.nbigb.eq.0)then !it is a full length model
ccc   check loop position & loop number
         nloop=0
         if(sec(1).eq.1)then
            nloop=nloop+1
            loop_i(nloop)=1
            loop_f(nloop)=1
         endif
         do i=2,Lch
            if(sec(i).eq.1.and.sec(i-1).ne.sec(i))then
               nloop=nloop+1
               loop_i(nloop)=i
            endif
            if(sec(i).eq.1)then
               loop_f(nloop)=i
            endif
         enddo
         do i=1,nloop
            loop_p(i)=int(float(loop_f(i)+loop_i(i))/2.0)
         enddo
ccc   cut mloop by setting q(i)=0
         if(mloop.gt.0)then
            delt=float(nloop)/float(mloop+1)
            if(delt.lt.1)delt=1
            kk=1
            do i=1,nloop
               if(i.ge.delt*kk)then
                  m=loop_p(i)
                  q(m)=0
                  if(kk.ge.mloop)goto 50
                  kk=kk+1
               endif
            enddo
         endif
      endif
 50   continue
c^^^^^^^^^^^^^^^^ cut loop is finished ^^^^^^^^^^^^^^^^^^^

c      do i=1,Lch
c         write(*,*)i,sec(i),q_old(i),q(i)
c      enddo
c      stop
**********************************************************
c     using predicted secondary structures
**********************************************************
      if(ssp.eq.1)then
         call secondary         !redifine q(i) and cx0(i)
      endif

********************************************************
c     remove small segments, so that we have a better random walk
********************************************************
      M_a=0                     !number of aligned segments
      q(0)=0
      do i=1,Lch
         if(q(i).ne.q(i-1).and.q(i).eq.1)then
            M_a=M_a+1
            M_i(M_a)=i
         endif
         if(q(i).eq.1)M_f(M_a)=i
      enddo
      do i=1,M_a
         L_a=M_f(i)-M_i(i)+1
         if(L_a.le.L_cut)then
            do j=M_i(i),M_f(i)
               q(j)=0
            enddo
         endif
      enddo
c^^^^^^^^^^^ remove small segment finished ^^^^^^^^^^^^^

********************************************************
c     check and find GAP, only according to q(i)
********************************************************
      N_g=0                     !number of gaps
      q(0)=1
      do i=1,Lch
         if(q(i).ne.q(i-1).and.q(i).eq.0)then
            N_g=N_g+1
            n_i(N_g)=i          !initial point of the gap
         endif
         if(q(i).eq.0)n_f(N_g)=i !final point of the gap
      enddo
c^^^^^^^^^^^^^^^ check gap finished ^^^^^^^^^^^^^^^^^

      n_check_bond=0
 103  continue                  !re-walk and try to make dis(i,i+1)<7A
********************************************************
c     fill GAP: cx0(i) -> cx(i)
********************************************************
      n_walk=0                  !for number of reject by excluded volumn
 70   n_walk=n_walk+1
      if(n_walk.gt.100)then
         exc_eff=exc_eff*0.99
      endif
      if(n_walk.gt.10000)then
         write(20,*)'unsolvable structure',pass,i,exc_eff
         write(*,*)'unsolvable structure',pass,i,exc_eff
         stop
      endif
      do i=1,Lch
         if(q(i).eq.1)then
            cx(i)=cx0(i)
            cy(i)=cy0(i)
            cz(i)=cz0(i)
         else
            cx(i)=1000000.      !for checking excluded volumn
            cy(i)=1000000.
            cz(i)=1000000.
         endif
      enddo
      do i=2,n_g
         i1=n_i(i)
         i2=n_f(i)
         call connect(i1,i2,pass) !fill missed cooridinates, no move others
         if(pass.ge.3)goto 70             !re-walk
      enddo
      if(n_g.ge.1)then
         i1=n_i(1)
         i2=n_f(1)
         call connect(i1,i2,pass)
         if(pass.ge.3)goto 70   !re-walk
      endif
*^^^^^^^^^^^^^^Fill gap done, cx(i) is continuous ^^^^^^^^^^^^^^^^^^^^^^^

**********************************************************************
c     decide moveable points nfl, ras(i), according to 3 factors:
**********************************************************************
***   1: ras => gap +- 2, smallest moveable-length is 3
      nfl=0                     !number of moveable points
      do i=1,Lch
         ras(i)=0               !residue name of i'th moveable point.
      enddo
      do i=1,n_g
         i1=n_i(i)-2
         i2=n_f(i)+2
         if(i1.lt.1)i1=1
         if(i2.gt.Lch)i2=Lch
         do ii=i1,i2
            nfl=nfl+1           !number of moveable residues
            ras(nfl)=ii         !locations
         enddo
      enddo
      if(nfl.gt.0)call sort_ras
      if(i_bigbond.eq.1)then    !only for medium
***   2a: ras => bigbond on  +-2 for templates:
         do i=2,Lch
            if(q(i-1).eq.1.and.q(i).eq.1)then
               bondlen=di(cx(i-1),cy(i-1),cz(i-1),cx(i),cy(i),cz(i))
               if(bondlen.gt.4.6)then
                  ik_f=i+1      !ending point
                  ik_i=i-2      !starting point of big-bond
                  if(ik_i.lt.1)ik_i=1
                  if(ik_f.gt.Lch)ik_f=Lch
                  do ik=ik_i,ik_f
                     nfl=nfl+1
                     ras(nfl)=ik
                  enddo
               endif
            endif
         enddo
         if(nfl.gt.0)call sort_ras
      else
***   2b: ras => bigbond +-1
         bond_max=0
         do i=2,Lch
            j=i-1               !starting point of big-bond
            bondlen=di(cx(i),cy(i),cz(i),cx(j),cy(j),cz(j))
            if(bondlen.gt.bond_max) bond_max=bondlen
            if(bondlen.gt.5.0)then
               do t=i,Lch
                  dh=di(cx(j),cy(j),cz(j),cx(t),cy(t),cz(t))
                  adh=dh/float(t-j)
                  if(adh.lt.3.498)goto 23 !can be walked to
               enddo
 23            ik_f=t+1         !ending point + 1
               ik_i=j-1         !starting point of big-bond - 1
               if(ik_i.lt.1)ik_i=1
               if(ik_f.gt.Lch)ik_f=Lch
               do ik=ik_i,ik_f
                  nfl=nfl+1
                  ras(nfl)=ik
               enddo
            endif
         enddo
         if(nfl.gt.0)call sort_ras
      endif
***   3: ras0 => smallpiece:
      do i=1,Lch
         mf(i)=0                !frozen
      enddo
      do i=1,nfl
         mf(ras(i))=1           !moveable points
      enddo
      nfr=0                     !number of frozen fragments
      mf(0)=1
      do i=1,Lch
         if(mf(i).ne.mf(i-1).and.mf(i).eq.0)then
            nfr=nfr+1
            nfr_i(nfr)=i        !starting point of nfr'th fragment
         endif
         if(mf(i).eq.0)nfr_f(nfr)=i !ending point of nfr'th fragment
      enddo
c     'nfr' frozen fragments in [nfr_i(i),nfr_f(i)], convert small piece ---->
      do i=1,nfr
         l_fr=nfr_f(i)-nfr_i(i)+1 !length of the frozen fragment
         if(l_fr.le.L_cut)then
            do j=nfr_i(i),nfr_f(i)
               nfl=nfl+1
               ras(nfl)=j
            enddo
         endif
      enddo
      if(nfl.gt.0)call sort_ras
*^^^^^^^^^^^^^^^^^^^^ nfl, ras(i) finished ^^^^^^^^^^^^^^^^^^^^^^

*****************************************************************
c     decide nfl2, ras2(i), mv(i):
*****************************************************************
      if(nfl.lt.1)then
         write(20,*)'without moveable points in this template!'
         write(*,*)'without moveable points in this template!'
      endif
      call move_point           !decide movement set from ras(i)
*^^^^^^^^^^ decision of movement set finished ^^^^^^^^^^^^^^^^^^^

****************************************************************************
c     decide nfr_i(i), nfr_f(i), i=1,...,nfr, according to nfl,ras(nfl):
****************************************************************************
      do i=1,Lch
         mf(i)=0
         sign(i)="f"            !frozen
      enddo
      do i=1,nfl
         mf(ras(i))=1
         sign(ras(i))="m"       !moveable points
      enddo
      mf(0)=1
      nfr=0                     !number of frozen fragments
      do i=1,100
         nfr_i(i)=0
         nfr_f(i)=0
      enddo
      do i=1,Lch
         if(mf(i).ne.mf(i-1).and.mf(i).eq.0)then
            nfr=nfr+1
            nfr_i(nfr)=i         !starting point of nfr'th fragment
         endif
         if(mf(i).eq.0)nfr_f(nfr)=i !ending point of nfr'th fragment
      enddo
c^^^^^^^^^^^^^^^^^^^ nfr, nfr_i(i), nfr_f(i) finished ^^^^^^^^^^^^^^^^^^

****************************************************************
c     redefine mv(i) for frozen fragments (for excluded volumn)
****************************************************************
      do i=1,nfr
         do j=nfr_i(i),nfr_f(i)
            mv(j)=-i
         enddo
      enddo
c^^^^^^^^^^^^^^^^^ mv(i) finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

****************************************************************
c     translation distance d_xyz00(i), rotation angle00(i)
      do i=1,nfr
         siz=nfr_f(i)-nfr_i(i)+1
         angle00(i)=angle0*10/float(siz)
         d_xyz00(i)=d_xyz0*10/float(siz)
         if(siz.gt.40)d_xyz00(i)=d_xyz00(i)/2.0
      enddo
c^^^ d_xyz00(i), angle00(i) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

****************************************************************
c     Record frozen parts:
c     there are 4 arraies needed to record and transfer:
c     common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)    CA
c     common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) SG
c     common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) cc
c     common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) bb
c     common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep) CB
****************************************************************
      cx0(0)=cx(1)+(cx(2)-cx(3))
      cy0(0)=cy(1)+(cy(2)-cy(3))
      cz0(0)=cz(1)+(cz(2)-cz(3))
      cx0(Lch+1)=cx(Lch)+(cx(Lch1)-cx(Lch2))
      cy0(Lch+1)=cy(Lch)+(cy(Lch1)-cy(Lch2))
      cz0(Lch+1)=cz(Lch)+(cz(Lch1)-cz(Lch2))
      do itemp=1,N_rep
         do i=1,Lch
            if(mv(i).lt.0)then  !frozen points
c     Ca:
               exrep(i,itemp)=cx0(i)/0.87
               eyrep(i,itemp)=cy0(i)/0.87
               ezrep(i,itemp)=cz0(i)/0.87
c     uniform vector:
               amx=cx0(i)-cx0(i-1)
               amy=cy0(i)-cy0(i-1)
               amz=cz0(i)-cz0(i-1)
               aaa=sqrt(amx**2+amy**2+amz**2)
               amx=amx/aaa      !ax(i)-ax(i-1)
               amy=amy/aaa
               amz=amz/aaa
               apx=cx0(i+1)-cx0(i)
               apy=cy0(i+1)-cy0(i)
               apz=cz0(i+1)-cz0(i)
               aaa=sqrt(apx**2+apy**2+apz**2)
               apx=apx/aaa      !ax(i+1)-ax(i)
               apy=apy/aaa
               apz=apz/aaa
               ang=acos(-(amx*apx+amy*apy+amz*apz))*180/3.1415926 !
               aaax=amx+apx
               aaay=amy+apy
               aaaz=amz+apz
               aaa=sqrt(aaax**2+aaay**2+aaaz**2)
               aaax=aaax/aaa
               aaay=aaay/aaa
               aaaz=aaaz/aaa
c     cc:
               ccx=amx-apx
               ccy=amy-apy
               ccz=amz-apz
               aaa=sqrt(ccx**2+ccy**2+ccz**2)
               ccx=ccx/aaa
               ccy=ccy/aaa
               ccz=ccz/aaa
               ecxrep(i,itemp)=ccx !cc
               ecyrep(i,itemp)=ccy
               eczrep(i,itemp)=ccz
c     Hb:
               bbx=amy*apz-amz*apy
               bby=amz*apx-amx*apz
               bbz=amx*apy-amy*apx
               aaa=sqrt(bbx**2+bby**2+bbz**2)
               bbx=bbx/aaa
               bby=bby/aaa
               bbz=bbz/aaa
               ebxrep(i,itemp)=bbx !Hb
               ebyrep(i,itemp)=bby
               ebzrep(i,itemp)=bbz
c     CB:
               k=seq(i)
               if(ang.lt.105)then !alpha
                  dx=(bxalf(k)*aaax+byalf(k)*bbx+bzalf(k)*ccx)/0.87
                  dy=(bxalf(k)*aaay+byalf(k)*bby+bzalf(k)*ccy)/0.87
                  dz=(bxalf(k)*aaaz+byalf(k)*bbz+bzalf(k)*ccz)/0.87
               else
                  dx=(bxbet(k)*aaax+bybet(k)*bbx+bzbet(k)*ccx)/0.87
                  dy=(bxbet(k)*aaay+bybet(k)*bby+bzbet(k)*ccy)/0.87
                  dz=(bxbet(k)*aaaz+bybet(k)*bbz+bzbet(k)*ccz)/0.87
               endif
               etxrep(i,itemp)=exrep(i,itemp)+dx
               etyrep(i,itemp)=eyrep(i,itemp)+dy
               etzrep(i,itemp)=ezrep(i,itemp)+dz
c     SG:
               k=seq(i)
               if(ang.lt.105)then !alpha
                  dx=(axalf(k)*aaax+ayalf(k)*bbx+azalf(k)*ccx)/0.87
                  dy=(axalf(k)*aaay+ayalf(k)*bby+azalf(k)*ccy)/0.87
                  dz=(axalf(k)*aaaz+ayalf(k)*bbz+azalf(k)*ccz)/0.87
               else
                  dx=(axbet(k)*aaax+aybet(k)*bbx+azbet(k)*ccx)/0.87
                  dy=(axbet(k)*aaay+aybet(k)*bby+azbet(k)*ccy)/0.87
                  dz=(axbet(k)*aaaz+aybet(k)*bbz+azbet(k)*ccz)/0.87
               endif
               egxrep(i,itemp)=exrep(i,itemp)+dx
               egyrep(i,itemp)=eyrep(i,itemp)+dy
               egzrep(i,itemp)=ezrep(i,itemp)+dz
            endif
         enddo
      enddo
c     current ex(i) for check of excluded volumn:
      do i=1,Lch
         if(mv(i).lt.0)then     !frozen points
            ex(i)=cx0(i)/0.87
            ey(i)=cy0(i)/0.87
            ez(i)=cz0(i)/0.87
         endif
      enddo
c^^^^^^^^^^^^^^^ record of frozen-fragment are set^^^^^^^^^^^^^^^^^^^^^^

*************************************************************************
c     project chain onto lattices, decide (x,y,z):
*************************************************************************
***   ax(i)=cx(i)/0.87, for the decision of (x,y,z):
      do i=1,Lch
         ax(i)=cx(i)/0.87       !C_alpha scaled by 0.87
         ay(i)=cy(i)/0.87
         az(i)=cz(i)/0.87
      enddo
c     (ax,ay,az) into (x,y,z)----------------------->
      x(1)=nint(ax(1))
      y(1)=nint(ay(1))
      z(1)=nint(az(1))
      xm=x(1)
      ym=y(1)
      zm=z(1)
      do 101 i=2,Lch
         dis_min=100000000000.  !minimum distance between c_alpha and lattice
         j_ch=0                 !chosen vector
         do 100 j=1,nvecr
            if(i.gt.2)then      !check good neighbor
               if(.not.goodc(jm,j))goto 100
            endif
            x_tmp=xm+vx(j)
            y_tmp=ym+vy(j)
            z_tmp=zm+vz(j)
c     check excluded volumn---->
            if(mv(i).gt.0)then  !on-lattice
               do m=1,Lch       !!!!! off-lattice part.
                  if(abs(i-m).ge.3.and.mv(m).lt.0)then
                     disaa=(x_tmp-ex(m))**2+(y_tmp-ey(m))**2+
     $                    (z_tmp-ez(m))**2
                     if(disaa.lt.exc_eff) goto 100
                  endif
               enddo
               do m=1,i-3       !!!!! on-lattice part.
                  if(mv(m).gt.0)then
                     disaa=(x_tmp-x(m))**2+(y_tmp-y(m))**2+
     $                    (z_tmp-z(m))**2
                     if(disaa.lt.exc_eff) goto 100
                  endif
               enddo
            endif
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^
            dis=(ax(i)-x_tmp)**2+(ay(i)-y_tmp)**2+(az(i)-z_tmp)**2
            if(dis.lt.dis_min)then
               j_ch=j
               x_ch=x_tmp
               y_ch=y_tmp
               z_ch=z_tmp
               dis_min=dis
            endif
 100     continue
         if(j_ch.lt.1)goto 70   !refill the gaps
         x(i)=x_ch              !get (x(i),y(i),z(i)) here
         y(i)=y_ch
         z(i)=z_ch
         jm=j_ch
         xm=x(i)
         ym=y(i)
         zm=z(i)
 101  continue
c^^^^^^^^^^^^^^^^^^^project lattice done ^^^^^^^^^^^^^^^^^^^^^^^

***********************************************************
*     check the bond-length
***********************************************************
      do i=1,Lch1
         r2=
     $        (nint(aax(i))-nint(aax(i+1)))**2+
     $        (nint(aay(i))-nint(aay(i+1)))**2+
     $        (nint(aaz(i))-nint(aaz(i+1)))**2
         if(r2.gt.65)then       !dis>7A
            n_check_bond=n_check_bond+1
            if(n_check_bond.lt.100)goto 103
         endif
c         write(*,*)i,r2,sqrt(float(r2))*.87
      enddo
cccc  

*************************************************************************
c     record replcas of (x,y,z), icarep:
*************************************************************************
      do itemp=1,N_rep
         do i=1,Lch
            xrep(i,itemp)=x(i)  !x(i)
            yrep(i,itemp)=y(i)
            zrep(i,itemp)=z(i)
            if(i.lt.Lch)then
               wx=x(i+1)-x(i)   !ica(i)
               wy=y(i+1)-y(i)
               wz=z(i+1)-z(i)
               icarep(i,itemp)=vector(wx,wy,wz)
            endif
         enddo
         icarep(Lch,itemp)=icarep(Lch2,itemp) !just for unity
      enddo

      call get_vvv !get vvv(i,j). vvv(i,j)>0 check; vvv(i,j)<0 not check

***********************************************************************
c      decide the frequence of bulk movements: f=1/n_fra
***********************************************************************
      n_fra=nint(8*sqrt(float(Lch)/100.0)) !move bulk once in n_fra local move
cccccc n_fra decrease with number of frozen residues
      nnn=0
      do i=1,Lch
        if(mv(i).lt.0)nnn=nnn+1
      enddo
      n_fra=nint(n_fra*(float(Lch)/(2.0*float(nnn)+0.001)))
cccccc n_fra decrease with number of frozen fragments:
      np=1
      if(nfr.eq.1)np=2
      n_fra=nint(n_fra*(2.8/(float(nfr)+0.001))**np)
      if(n_fra.lt.1) n_fra=1 !nfr big and nnn big
      if(n_fra.gt.20) n_fra=20 !nfr small and nnn small
c^^^^^^^^^^^ n_fra decided ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

********************************************************************
*     output the threading results:
********************************************************************
      write(20,*)'Lch=',Lch
      write(20,*)"#align=",Nal,'   #unalign=',Lch-Nal,"n_g=",n_g
      write(20,*)'#frozen=',Lch-nfl,'    #moveable=',nfl
      write(20,*)
      write(20,*)'number of frozen pieces: nfr=',nfr
      do i=1,nfr
         siz=nfr_f(i)-nfr_i(i)+1
         write(20,42)i,siz,nfr_i(i),nfr_f(i),d_xyz00(i),angle00(i)*57.3
      enddo
 42   format(i5,i5,' [',i3,',',i3,']',f8.3,f8.3)
      write(20,*)'---number of local movement for each bulk move:',n_fra
      write(20,*)
      write(20,*)'#chain fix/move  SEC   SEQ'
      do i=1,Lch
         write(20,41)i,sign(i),sec(i),sequ(i)
      enddo
 41   format(i7,a5,i8,a8)
      do i=1,Lch
         do j=1,Lch
            if(vvv(i,j).eq.-2)then
               write(20,*)i,j,'  not be checked for excluded volumn'
            endif
         enddo
      enddo
      write(20,*)'bond_max=',bond_max

*******for calculation of RMSD of fragment------------------------>
      do i=1,Lch
         i_chunk(i)=-1
      enddo
      do i=1,nfr
         do j=nfr_i(i),nfr_f(i)
            i_chunk(j)=i
            ex0(j)=ex(j)
            ey0(j)=ey(j)
            ez0(j)=ez(j)
         enddo
      enddo
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

c      open(23,file='tmp',status='unknown')
c      do i=1,Lch
c         write(*,*)i,mv(i),q(i)
c         if(mv(i).le.0)then
c            write(23,1037)i,sequ(i),i,aax(i)*0.87,aay(i)*0.87,aaz(i)*0.87
c         endif
c      enddo
c 1037 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
c      close(23)

c      stop
c^^^^^^^^^^^^^^^^^^^^^^^^^ template_initial finished ^^^^^^^^^^^^^^^^^^^^^
      return
      end

***********************************************************************
c     decided the pairs that excluded-volumn should not be checked:
c     1, pairs in same frozen fragments;
c     2, pairs in different frozen fragments but initially violated.
***********************************************************************
      subroutine get_vvv
      implicit integer (i-z)
      parameter(ndim=1999)      !maximum length of chain-length
      common/excluded/vvv(ndim,ndim)
      common/chainm/mv(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/looks/exc,exc1,exc2

      do i=1,Lch
         do j=1,Lch
            vvv(i,j)=1          !every pair should be checked
            if(mv(i).lt.0.and.mv(j).lt.0)then
               if(mv(i).eq.mv(j))then !belong to the same fragment
                  vvv(i,j)=-1   !Not be checked
               else
                  dist2=(aax(i)-aax(j))**2+(aay(i)-aay(j))**2
     $                 +(aaz(i)-aaz(j))**2
                  if(dist2.lt.exc)then !initially violated
                     vvv(i,j)=-2
                  endif
               endif
            endif
c            write(*,*)i,j,vvv(i,j)
         enddo
      enddo

c^^^^^^^^^^^^^^^^^^^ vvv(i,j) finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     get q(i), cx0(i) for consensus segments of top-2 templates.
c     if can not find consensus segment, back to top-1.
cSSSccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine secondary
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      character*3 sequ
      common/lengths/Lch,Lch1,Lch2
      common/consensus1/cx0(0:ndim),cy0(0:ndim),cz0(0:ndim),q(0:ndim)
      common/aminoacid/sequ(ndim)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension sec_high(0:ndim)
      dimension nse_i(200),nse_f(200),nse_type(200)

      dimension cx00(0:ndim),cy00(0:ndim),cz00(0:ndim)
      dimension ax(ndim),ay(ndim),az(ndim)

ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   

ccc   read SSP --------------------------------->
      open(unit=23,file='seq.dat',status='old') !ssp from native
      do i=1,Lch
         read(23,707)k,name,sec_high(i)
c     write(*,*)i,k,name,sec_high(i)
      enddo
      close(23)
 707  format(i5,3x,a3,2i5)

ccc   find continuous SSP pieces -------------------->
      nse=0                     !number of fragments
      do i=1,200
         nse_i(i)=0
         nse_f(i)=0
      enddo
      sec_high(0)=1
      do i=1,Lch
         if(sec_high(i).ge.2.and.sec_high(i).ne.sec_high(i-1))then
            nse=nse+1
            nse_i(nse)=i        !initial location of nse'th ssp
            nse_type(nse)=sec_high(i) !type of ssp
         endif
         if(sec_high(i).ne.1)nse_f(nse)=i !final location of nse'th ssp
      enddo

ccc   remove small secondary structures -------------->
      nse_old=nse
      nse=0
      do i=1,nse_old
         len=nse_f(i)-nse_i(i)+1
         if(len.ge.L_cut)then
            nse=nse+1
            nse_i(nse)=nse_i(i)
            nse_f(nse)=nse_f(i)
            nse_type(nse)=nse_type(i)
         endif
      enddo

c      do i=1,nse
c         write(*,*)i,nse_i(i),nse_f(i),nse_f(i)-nse_i(i)+1,nse_type(i)
c      enddo

ccc   generate new (cx0,cy0,cz0) --------------------------->
      do i=1,Lch
         if(q(i).eq.1)then
            cx00(i)=cx0(i)
            cy00(i)=cy0(i)
            cz00(i)=cz0(i)
         endif
      enddo
      do 10 i=1,nse             !nse of SSP
         len=nse_f(i)-nse_i(i)+1 !length of SSP

ccc   generate standard fragments --------------------->
         if(nse_type(i).eq.2)then
            call alpha_helix(ax,ay,az,len) !alpha-helix
         else
            call beta_sheet(ax,ay,az,len) !beta-sheet
         endif
         
ccc   rotation matrix based on initial template ------->
*1: aligned
         n=0
         do j=nse_i(i),nse_f(i)
            if(q(j).eq.1)then
               n=n+1
               k=j-nse_i(i)+1
               r_1(1,n)=ax(k)
               r_1(2,n)=ay(k)
               r_1(3,n)=az(k)
               r_2(1,n)=cx00(j)
               r_2(2,n)=cy00(j)
               r_2(3,n)=cz00(j)
            endif
         enddo
*2: gapped at middle
         if(n.lt.2)then         !it is a gap at ssp region
            r_1(1,1)=ax(1)
            r_1(2,1)=ay(1)
            r_1(3,1)=az(1)
            r_1(1,2)=ax(nse_f(i)-nse_i(i)+1)
            r_1(2,2)=ay(nse_f(i)-nse_i(i)+1)
            r_1(3,2)=az(nse_f(i)-nse_i(i)+1)
            n=0
            j=nse_i(i)
 11         j=j-1
            if(q(j).ne.1.and.j.gt.1)goto 11
            if(q(j).eq.1)then
               n=n+1
               ax1=cx00(j)
               ay1=cy00(j)
               az1=cz00(j)
               j0=j
               mm=1
            endif
            j=nse_f(i)
 12         j=j+1
            if(q(j).ne.1.and.j.lt.Lch)goto 12
            if(q(j).eq.1)then
               n=n+1
               ax2=cx00(j)
               ay2=cy00(j)
               az2=cz00(j)
               mm=2
             endif
             if(n.eq.2)then
               aaa=sqrt((ax2-ax1)**2+(ay2-ay1)**2+(az2-az1)**2)
               bbb=sqrt((ax(len)-ax(1))**2+(ay(len)-ay(1))**2+
     &              (az(len)-az(1))**2)
               ccc=abs(nse_i(i)-j0)*1
               r_2(1,1)=ax1+(ax2-ax1)/aaa*ccc
               r_2(2,1)=ay1+(ay2-ay1)/aaa*ccc
               r_2(3,1)=az1+(az2-az1)/aaa*ccc
               r_2(1,2)=ax1+(ax2-ax1)/aaa*(ccc+bbb)
               r_2(2,2)=ay1+(ay2-ay1)/aaa*(ccc+bbb)
               r_2(3,2)=az1+(az2-az1)/aaa*(ccc+bbb)
             endif
         endif
*3: gap at terminal
         if(n.lt.2)then         !n=1, the ssp gap is at terminal
            n=2
            if(mm.eq.1)then     !C-terminal
*3a: C-terminal
               m=0
               j=nse_i(i)
 13            j=j-1
               if(q(j).eq.1)m=m+1
               if(m.eq.1.and.q(j).eq.1)then
                  j0=j
                  ax2=cx00(j)
                  ay2=cy00(j)
                  az2=cz00(j)
               endif
               if(m.lt.5)goto 13 !the minimum alignment length is 5
               if(m.eq.5.and.q(j).eq.1)then
                  ax1=cx00(j)
                  ay1=cy00(j)
                  az1=cz00(j)
               endif
               aaa=sqrt((ax2-ax1)**2+(ay2-ay1)**2+(az2-az1)**2)
               bbb=sqrt((ax(len)-ax(1))**2+(ay(len)-ay(1))**2+
     &              (az(len)-az(1))**2)
               ccc=abs(nse_i(i)-j0)*2
               r_2(1,1)=ax2+(ax2-ax1)/aaa*ccc
               r_2(2,1)=ay2+(ay2-ay1)/aaa*ccc
               r_2(3,1)=az2+(az2-az1)/aaa*ccc
               r_2(1,2)=ax2+(ax2-ax1)/aaa*(ccc+bbb)
               r_2(2,2)=ay2+(ay2-ay1)/aaa*(ccc+bbb)
               r_2(3,2)=az2+(az2-az1)/aaa*(ccc+bbb)
            else
*3a: N-terminal
               m=0
               j=nse_f(i)
 14            j=j+1
               if(q(j).eq.1)m=m+1
               if(m.eq.1.and.q(j).eq.1)then
                  j0=j
                  ax1=cx00(j)
                  ay1=cy00(j)
                  az1=cz00(j)
               endif
               if(m.lt.5)goto 14 !the minimum alignment length is 5
               if(m.eq.5.and.q(j).eq.1)then
                  ax2=cx00(j)
                  ay2=cy00(j)
                  az2=cz00(j)
               endif
               aaa=sqrt((ax2-ax1)**2+(ay2-ay1)**2+(az2-az1)**2)
               bbb=sqrt((ax(len)-ax(1))**2+(ay(len)-ay(1))**2+
     &              (az(len)-az(1))**2)
               ccc=abs(j0-nse_f(i))*2
               r_2(1,2)=ax1+(ax1-ax2)/aaa*ccc
               r_2(2,2)=ay1+(ay1-ay2)/aaa*ccc
               r_2(3,2)=az1+(az1-az2)/aaa*ccc
               r_2(1,1)=ax1+(ax1-ax2)/aaa*(ccc+bbb)
               r_2(2,1)=ay1+(ay1-ay2)/aaa*(ccc+bbb)
               r_2(3,1)=az1+(az1-az2)/aaa*(ccc+bbb)
            endif
         endif
c         write(*,*)'n=',n,'  i=',i
c         do j=1,n
c            write(*,*)j,r_2(1,j),r_2(2,j),r_2(3,j)
c         enddo
         call u3b(w,r_1,r_2,n,1,rms,u,t,ier) !u rotate r_1 to r_2
ccc   rotate ssp onto initial template ------------>
         do j=nse_i(i),nse_f(i)
            k=j-nse_i(i)+1
            cx0(j)=t(1)+u(1,1)*ax(k)+u(1,2)*ay(k)+u(1,3)*az(k)
            cy0(j)=t(2)+u(2,1)*ax(k)+u(2,2)*ay(k)+u(2,3)*az(k)
            cz0(j)=t(3)+u(3,1)*ax(k)+u(3,2)*ay(k)+u(3,3)*az(k)
         enddo
 10   continue

ccc   re-define q(i) ------------------------------>
      do i=1,Lch
         q(i)=0
      enddo
      do i=1,nse
         do j=nse_i(i),nse_f(i)
            q(j)=1
         enddo
      enddo

c      open(23,file='tmp',status='unknown')
c      do i=1,Lch
c         if(q(i).eq.1)then
c            write(23,1237)i,sequ(i),i,cx0(i),cy0(i),cz0(i)
c         endif
c      enddo
c 1237 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
c      close(23)

c      stop
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     generate standard alpha-helix:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine alpha_helix(x,y,z,n)
      parameter(ndim=1999)
      dimension x(ndim),y(ndim),z(ndim)

      rad=2.3                   !redius of helix
      do i=1,n
         angle=100*3.1415926/180*(i-1) !100 degree per residues
         x(i)=rad*cos(angle)
         y(i)=rad*sin(angle)
         z(i)=1.5*(i-1)            !increase 1.5 per residue
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     generate standard beta-sheet:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine beta_sheet(x,y,z,n)
      parameter(ndim=1999)
      dimension x(ndim),y(ndim),z(ndim)

      dist2=6.70098             !distance between i and i+2
      high=1.792820             !distance between i+1 and middle of (i,i+2)
      do i=1,n
         x(i)=dist2/2*(i-1)
         if(int(i/2)*2.eq.i)then
            y(i)=0
         else
            y(i)=high
         endif
         z(i)=0
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     get q(i), cx0(i) for consensus segments of top-2 templates.
c     if can not find consensus segment, back to top-1.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_consensus
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/consensus1/cx0(0:ndim),cy0(0:ndim),cz0(0:ndim),q(0:ndim)
      dimension cx1(ndim),cy1(ndim),cz1(ndim),q1(ndim),ip1(ndim)
      dimension cx2(ndim),cy2(ndim),cz2(ndim),q2(ndim),ip2(ndim)
      dimension qq(ndim)
      real xx,yy,zz
ccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   

******************************************************************
***   read templates------------->
      rewind(24)
      read(24,*)N_tmp
      if(N_tmp.lt.2)then
         write(20,*)'There is only one template'
         write(20,*)'Top-1 template is used'
         return                 !use the first template
      else
***   read first template:
         read(24,*)N_al
         do i=1,N_al
            read(24,1237)text,ii,text,a1,a2,a3
            cx1(ii)=a1
            cy1(ii)=a2
            cz1(ii)=a3
            q1(ii)=1
         enddo
         read(24,*)text
***   read second template:
         read(24,*)N_al
         do i=1,N_al
            read(24,1237)text,ii,text,a1,a2,a3
            cx2(ii)=a1
            cy2(ii)=a2
            cz2(ii)=a3
            q2(ii)=1
         enddo
         read(24,*)text
***
      endif
 1237 format(A22,I4,A4,3F8.3)

******************************************************************
***   decided qq(i)------------->
      do i=1,Lch
         qq(i)=0
      enddo
      k=0
      do i=1,Lch
         if(q1(i).eq.1.and.q2(i).eq.1)then
            k=k+1
            ip1(k)=i
            r_1(1,k)=cx1(i)
            r_1(2,k)=cy1(i)
            r_1(3,k)=cz1(i)
            r_2(1,k)=cx2(i)
            r_2(2,k)=cy2(i)
            r_2(3,k)=cz2(i)
         endif
      enddo
      write(*,*)"#common aligned=",k
      if(k.lt.10)then
         write(20,*)'There is less than 10 common aligned residues'
         write(20,*)'Top-1 template is used'
         return                 !no common aligned points, using template1
      endif
      call u3b(w,r_1,r_2,k,1,rms,u,t,ier) !u rotate r_1 to r_2
      armsd=dsqrt(rms/k)        !RMSD12
      write(20,*)'RMSD1=',armsd,k !RMSD between template1,2 for common align
      kk=0
      do j=1,k
         xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
         yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
         zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
         dis=sqrt((xx-r_2(1,j))**2+(yy-r_2(2,j))**2+(zz-r_2(3,j))**2)
c         write(*,*)j,ip1(j),dis
         if(dis.lt.5)then
            kk=kk+1
            qq(ip1(j))=1
         endif
      enddo
      if(kk.lt.10)then
         write(20,*)'There is less than 10 close common residues'
         write(20,*)'Top-1 template is used'
         return                 !no consensus points, using template1
      endif

******************************************************************
***   q(i)=qq(i), cx0(i)=(cx1+cx2)/2--------------->
      do i=1,Lch
         q(i)=qq(i)
         cx0(i)=1000000.        !for checking excluded volumn
         cy0(i)=1000000.
         cz0(i)=1000000.
      enddo
      k=0
      do i=1,Lch
         if(q(i).eq.1)then
            k=k+1
            ip2(k)=i
            r_1(1,k)=cx1(i)
            r_1(2,k)=cy1(i)
            r_1(3,k)=cz1(i)
            r_2(1,k)=cx2(i)
            r_2(2,k)=cy2(i)
            r_2(3,k)=cz2(i)
         endif
      enddo
      call u3b(w,r_1,r_2,k,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,k
         xx=t(1)+u(1,1)*r_1(1,j)+u(1,2)*r_1(2,j)+u(1,3)*r_1(3,j)
         yy=t(2)+u(2,1)*r_1(1,j)+u(2,2)*r_1(2,j)+u(2,3)*r_1(3,j)
         zz=t(3)+u(3,1)*r_1(1,j)+u(3,2)*r_1(2,j)+u(3,3)*r_1(3,j)
         cx0(ip2(j))=(xx+r_2(1,j))/2
         cy0(ip2(j))=(yy+r_2(2,j))/2
         cz0(ip2(j))=(zz+r_2(3,j))/2
      enddo
********************************************************************
      L_cut=2                   !small segment to be shrowed away.

c^^^^^^^^^^q(i) and cx0(i) obtained ^^^^^^^^^^^^^^^^^^^^

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sort ras, recalculate nfl
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sort_ras
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/chain0/ras(ndim),nfl

      dimension ras0(ndim)

      ras_min=10000
      ras_max=-10000
      do i=1,nfl
         if(ras(i).gt.ras_max)ras_max=ras(i)
         if(ras(i).lt.ras_min)ras_min=ras(i)
      enddo

      nfl0=1
      ras0(nfl0)=ras_min
      do 1 while(ras0(nfl0).lt.ras_max)
         ras_min=10000
         do i=1,nfl
            if(ras(i).gt.ras0(nfl0))then
               if(ras(i).lt.ras_min)then
                  ras_min=ras(i)
               endif
            endif
         enddo
         nfl0=nfl0+1
         ras0(nfl0)=ras_min
 1    continue

      nfl=nfl0
      do i=1,nfl
         ras(i)=ras0(i)
      enddo

*^^^^^^^^^^^^^^ sort ras done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     (x,y,z) --> (r,thita,phi)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine triangle(bx,by,bz,br,bthita,bphi)
      parameter(bpi=3.1415926)

      br=sqrt(bx*bx+by*by+bz*bz)
      bthita=acos(bz/br)
      if(abs(bx).gt.0.00001)then
         bphi=atan(abs(by/bx))
      else
         bphi=0
      endif
      if(bx.gt.0)then
         if(by.lt.0)bphi=bphi+bpi*1.5
      else
         if(by.gt.0)then
            bphi=bphi+bpi*0.5
         else
            bphi=bphi+bpi
         endif
      endif

c^^^^^^^^^^^^^^^^^^ triangle done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccc decide movement point from ras(i) cccccccccccccccccccc
c     nfl----total number of moveable points on model-chain
c     nfl2---total number of 2-bond-movement
c     nfl3---total number of 3-bond-movement
c     nfl4---total number of 4-bond-movement
c     nfl5---total number of 5-bond-movement
c     ras(i)---position of i-th moveable residues
c     ras2(i)--initial position of i-th 2-bond movement.
c     ras3(i)--initial position of i-th 3-bond movement.
c     ras4(i)--initial position of i-th 4-bond movement.
c     ras5(i)--initial position of i-th 5-bond movement.
c     ras6(i)--initial position of i-th 6-bond movement.
cccccccccccccc
      subroutine move_point
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/chain0/ras(ndim),nfl
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/chainm/mv(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      
      dimension mv1(5000)
      
      do i=1,Lch+10
         mv1(i)=-1              !frozen point
      enddo
      
      do i=1,nfl
         mv1(ras(i))=1          !moveable point
      enddo
      
c     re-decide the movement range of tremendicy:
      Mend_N=0                  !Mend_N points from 1 are moveable
      k=0
      do i=1,Lch
         k=k+1
         if(k.gt.Mend.or.mv1(i).lt.0) goto 111
         Mend_N=Mend_N+1
      enddo
 111  continue
      Mend_C=0                  !Mend_C points from Lch are moveable
      k=0
      do i=Lch,1,-1
         k=k+1
         if(k.gt.Mend.or.mv1(i).lt.0) goto 222
         Mend_C=Mend_C+1
      enddo
 222  continue

c     find 6-bond-move position -------------->
      nfl6=0                    !total number of 6-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then
                  if(mv1(i+3).gt.0)then
                     if(mv1(i+4).gt.0)then
                        if(mv1(i+5).gt.0)then
                           if(mv1(i+6).gt.0)then !i+6: fixed border
                              nfl6=nfl6+1
                              ras6(nfl6)=i
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo

c     find 5-bond-move position -------------->
      nfl5=0                    !total number of 5-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then
                  if(mv1(i+3).gt.0)then
                     if(mv1(i+4).gt.0)then
                        if(mv1(i+5).gt.0)then !i+5: fixed border
                           nfl5=nfl5+1
                           ras5(nfl5)=i
                        endif
                     endif
                  endif
               endif
            endif
         endif
      enddo

c     find 4-bond-move position -------------->
      nfl4=0                    !total number of 4-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then
                  if(mv1(i+3).gt.0)then
                     if(mv1(i+4).gt.0)then !i+4: fixed border
                        nfl4=nfl4+1
                        ras4(nfl4)=i
                     endif
                  endif
               endif
            endif
         endif
      enddo

c     find 3-bond-move position -------------->
      nfl3=0                    !total number of 3-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then
                  if(mv1(i+3).gt.0)then !i+3 fixed border
                     nfl3=nfl3+1
                     ras3(nfl3)=i
                  endif
               endif
            endif
         endif
      enddo

c     find 2-bond-move position -------------->
      nfl2=0                    !total number of 2-bond-move points
      do i=1,Lch                !i: fixed border of moving block
         if(mv1(i).gt.0)then
            if(mv1(i+1).gt.0)then
               if(mv1(i+2).gt.0)then !i+2 fixed border
                  nfl2=nfl2+1
                  ras2(nfl2)=i
               endif
            endif
         endif
      enddo
      
      do i=1,Lch
         mv(i)=mv1(i)
      enddo
      
c^^^^^^^^^^^^^^^^^^^^ move_point finished ^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc connect gaps in [i1,i2] for (bx,by,bz):
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine connect(i1,i2,pass)
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/initialstr/cx(ndim),cy(ndim),cz(ndim)
      common/mcheck_dis/amcheck_dis

      bdis0=3.5
      amcheck_dis=3.1            !3.1*0.87=2.7
      pass=1
      if(i1.eq.1)then           !!!!!N-terminal or whole structure, random walk
         do i=i2,1,-1
            n_check=0
 10         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i+1)+axx
            cy(i)=cy(i+1)+ayy
            cz(i)=cz(i+1)+azz
            n_check=n_check+1
            if(n_check.gt.100000)then
               pass=3
               return
            endif
            if(mcheck(i,cx(i),cy(i),cz(i)).eq.3) goto 10 !excluded volumn
         enddo
      elseif(i2.eq.Lch)then     !!!!!!!!!!!!!!!!!!C_terminal,
         do i=i1,Lch
            n_check=0
 11         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i-1)+axx
            cy(i)=cy(i-1)+ayy
            cz(i)=cz(i-1)+azz
            n_check=n_check+1
            if(n_check.gt.100000)then
               pass=4
               return
            endif
            if(mcheck(i,cx(i),cy(i),cz(i)).eq.3) goto 11
         enddo
      else                      !!!!!!!!!!!!!interval
         diss=di(cx(i2+1),cy(i2+1),cz(i2+1),cx(i1-1),cy(i1-1),cz(i1-1))
         adis=diss/float((i2+1)-(i1-1))
***   linear connect for big gap---->
         adis0=3.5
         x_check=0
 15      x_check=x_check+1
         if(adis.ge.adis0)then
            dex=3.5*(cx(i1-1)-cx(i2+1))/diss !3.5*cos(thita)
            dey=3.5*(cy(i1-1)-cy(i2+1))/diss
            dez=3.5*(cz(i1-1)-cz(i2+1))/diss
            do j=i2,i1,-1
               cx(j)=cx(j+1)+dex
               cy(j)=cy(j+1)+dey
               cz(j)=cz(j+1)+dez
            enddo
            return              !end of connection
         endif
***   random walk from i1 to i2--------->
         m_check=0              !try 2000 times of whole walk
 13      m_check=m_check+1
         do 14 i=i1,i2
            n_check=0           !each point try 2000 times
 12         call get_bond(axx,ayy,azz,3.8)
            cx(i)=cx(i-1)+axx
            cy(i)=cy(i-1)+ayy
            cz(i)=cz(i-1)+azz
            n_check=n_check+1
            if(int(n_check/1000)*1000.eq.n_check)then 
               bdis0=bdis0*1.03
               if(bdis0.ge.3.8)bdis0=3.8
               amcheck_dis=2.8  !2.5*0.87=2.18
               if(i1.eq.i2)then !artificial gap, skip when initial model has error
                  amcheck_dis=2.0
               endif
            endif
            if(n_check.eq.4000)then
               goto 13
            endif
            if(m_check.ge.2000)then !can not pass the connection
               if(x_check.le.6)then
                  adis0=adis0*0.995
                  goto 15
               endif
               pass=5
               return
            endif
            if(mcheck(i,cx(i),cy(i),cz(i)).eq.3) goto 12 !check excluded V
            aaa=float(i2+1-i)
            bdis=di(cx(i),cy(i),cz(i),cx(i2+1),cy(i2+1),cz(i2+1))/aaa
            if(i.lt.i2.and.bdis.ge.bdis0) goto 12 !check remain steps
            if(i.eq.i2)then
               if(bdis.gt.4.2.or.bdis.lt.3.4) goto 12 !last step
            endif
            bdis0=3.5
            amcheck_dis=3.1
 14      continue
      endif

c^^^^^^^^^^^^^^^^ connect of gap [i1,i2] finished ^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  produce a vector of length=3.8A:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_bond(axx,ayy,azz,al)
      implicit integer(i-z)
      common/ranzy/nozy
      athita=acos(1.-2.*aranzy(nozy)) !thita angle in random, [0,pi]
      aphi=2.*3.1415926*aranzy(nozy) !phi angle in random, [0,2pi]
      axx=al*sin(athita)*cos(aphi)
      ayy=al*sin(athita)*sin(aphi)
      azz=al*cos(athita)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  check weak excluded volumn for cx():
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function mcheck(i0,bx0,by0,bz0)
      parameter(ndim=1999)
      implicit integer(i-z)
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/initialstr/cx(ndim),cy(ndim),cz(ndim)
      common/mcheck_dis/amcheck_dis

      mcheck=1
      do i=1,Lch
         if(i0.ne.i)then
            dis=di(bx0,by0,bz0,cx(i),cy(i),cz(i)) !distance
            if(dis.le.amcheck_dis)then
               mcheck=3
               return
            endif
         endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function di(x1,y1,z1,x2,y2,z2)
      di=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function di2(x1,y1,z1,x2,y2,z2)
      di2=(x1-x2)**2+(y1-y2)**2+(z1-z2)**2
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function aax(i)
      implicit integer(i-z)
      parameter(ndim=1999)
      common/chainm/mv(ndim)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      if(mv(i).gt.0)then
         aax=x(i)
      else
         aax=ex(i)
      endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function aay(i)
      implicit integer(i-z)
      parameter(ndim=1999)
      common/chainm/mv(ndim)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      if(mv(i).gt.0)then
         aay=y(i)
      else
         aay=ey(i)
      endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function aaz(i)
      implicit integer(i-z)
      parameter(ndim=1999)
      common/chainm/mv(ndim)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      if(mv(i).gt.0)then
         aaz=z(i)
      else
         aaz=ez(i)
      endif
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     optimize the initial conformations from threading by RS:
cTTTTTTTTTTccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine template_simulation
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/xyzrs/xrs(ndim,40,40),yrs(ndim,40,40),zrs(ndim,40,40)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      dimension x0(ndim),y0(ndim),z0(ndim)  !for E_min
      character*6 mname
      common/movename/mname(100)
      common/nrepfile/n_repf

      common/temperature/itemp,atemp

      common/chain0/ras(ndim),nfl
      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/mng/m_g(100)

      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/thrtem/aT_ann(100)
      common/eall/N_sum(100),energ_sum(100),energ_sum2(100),E_min
      common/acct/accept0
      common/moveretio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
      common/moveretio3/bh7a,bh7b,bhendn,bhendc

      common/nswap/bNSa(100),bNSt(100)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t

      dimension aNNa(100),aNNt(100),N_out(100) !for output annealing

      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep)

      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/echain2/egx(ndim),egy(ndim),egz(ndim)
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)
      common/echain6/etx(ndim),ety(ndim),etz(ndim)
      common/chainm/mv(ndim)
      common/outputxyz/fxyz(3,ndim)
      character*3 sequ
      common/aminoacid/sequ(ndim)
      character text
      common/E_defo/i_E_defo
      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)
      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs
      common/ranzy/nozy
      common/hours/hour_max

      do i=1,100
         bNSa(i)=0              !aceptance for swep
         bNSt(i)=0
         bNa(i)=0               !acceptance for move2,3,4,5,6,7
         bNt(i)=0
         bNNa(i)=0              !acceptance for different temperature.
         bNNt(i)=0
         N_sum(i)=0
         energ_sum(i)=0         !<E_tot>
         energ_sum2(i)=0        !<E_tot^2>
         armsd_sum(i)=0
         N_rmsd(i)=0
      enddo
      E_min=10000

      i_tr=0                    !order number of output trajectory
      i_fra=0
      mcycle=0
      do 1111 icycle=1,ncycle
         do 2222 itemp=1,N_rep  !iterate for all the replicas
            atemp=aT_rep(itemp) !current temperature
            aTs=aTs_rep(itemp)
            aTTs=aTTs_rep(itemp)
            call set_current_RS !get current (x,y,z,ica,T)
            call initial_move   !update center, axis, energy
ccc   
            do 3333 iphot=1,phot !N_swap, iterate at fixed temperature
               do 4444 i_nfl=1,nfl
                  fff=aranzy(nozy)
                  if(fff.le.bh2)then
                     call move2
                  elseif(fff.le.bh3s)then
                     if(nfl3.ge.1)then
                        call move3s
                     else
                        call move2
                     endif
                  elseif(fff.le.bh3d)then
                     if(nfl3.ge.1)then
                        call move3d
                     else
                        call move2
                     endif
                  elseif(fff.le.bh4s)then
                     if(nfl4.ge.1)then
                        call move4s
                     else
                        call move2
                     endif
                  elseif(fff.le.bh4d)then
                     if(nfl4.ge.1)then
                        call move4d
                     else
                        call move2
                     endif
                  elseif(fff.le.bh5s)then
                     if(nfl5.ge.1)then
                        call move5s
                     else
                        call move2
                     endif
                  elseif(fff.le.bh5d)then
                     if(nfl5.ge.1)then
                        call move5d
                     else
                        call move2
                     endif
                  elseif(fff.le.bh6)then
                     if(nfl6.ge.1)then
                        call move6
                     else
                        call move2
                     endif
                  elseif(fff.le.bhendn)then
                     if(ras(1).eq.1)then
                        call move_n_end
                     else
                        call move2
                     endif
                  else
                     if(ras(nfl).eq.Lch)then
                        call move_c_end
                     else
                        call move2
                     endif
                  endif
ccccccccccccccccccccfragment movement ccccccccccccccccccccc
                  if(switch.eq.3)goto 4441 !freeze the template
                  i_fra=i_fra+1
                  if(i_fra.ge.n_fra.and.nfr.gt.0)then !1 chunk move in nfr move
                     i_fra=0
                     ifr=int(aranzy(nozy)*nfr)+1 ![1,nfr]
                     if(switch.eq.2)then !rotation+translation
                        if(nfr_i(ifr).eq.1)then
                           call trot_N !translate+rotate N-terminal
                        elseif(nfr_f(ifr).eq.Lch)then
                           call trot_C !translate+rotate C-terminal
                        else
                           call trot_M !translate+rotate Middle-fragment
                        endif
                     elseif(switch.eq.4)then !translation
                        if(nfr_i(ifr).eq.1)then
                           call tran_N !translate N-terminal
                        elseif(nfr_f(ifr).eq.Lch)then
                           call tran_C !translate C-terminal
                        else
                           call tran_M !translate Middle-fragment
                        endif
                     elseif(switch.eq.5)then !rotation+translation+deformation
                        aran_num=aranzy(nozy)
                        if(nfr_i(ifr).eq.1)then
                           if(aran_num.gt.0.5)then
                              i_E_defo=3
                              call trot_N !translate+rotate N-terminal
                           else
                              i_E_defo=1
                              call defo_N !translate+deform N-terminal
                           endif
                        elseif(nfr_f(ifr).eq.Lch)then
                           if(aran_num.gt.0.5)then
                              i_E_defo=3
                              call trot_C !translate+rotate C-terminal
                           else
                              i_E_defo=1
                              call defo_C !translate+deform C-terminal
                           endif
                        else
                           if(aran_num.gt.0.5)then
                              i_E_defo=3
                              call trot_M !translate+rotate Middle-fragment
                           else
                              i_E_defo=1
                              call defo_M !translate+deform Middle-fragment
                           endif
                        endif
                     endif
                  endif
 4441             continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                  atime=second()/3600.0
                  if(atime.gt.hour_max)goto 901
 4444           continue
 3333         continue
ccc
ccccccccccrecord energy and (x,y,z) cccccccccccccccccccccccc
            E_rep(itemp)=energy_tot() !whole energy
            if(E_rep(itemp).lt.E_min) E_min=E_rep(itemp)
            do i=1,Lch
               icarep(i,itemp)=ica(i)
               if(mv(i).gt.0)then
                  xrep(i,itemp)=x(i)
                  yrep(i,itemp)=y(i)
                  zrep(i,itemp)=z(i)
               else
                  exrep(i,itemp)=ex(i) !Ca
                  eyrep(i,itemp)=ey(i)
                  ezrep(i,itemp)=ez(i)
                  egxrep(i,itemp)=egx(i) !SG
                  egyrep(i,itemp)=egy(i)
                  egzrep(i,itemp)=egz(i)
                  ecxrep(i,itemp)=ecx(i) !cc
                  ecyrep(i,itemp)=ecy(i)
                  eczrep(i,itemp)=ecz(i)
                  ebxrep(i,itemp)=ebx(i) !Hb
                  ebyrep(i,itemp)=eby(i)
                  ebzrep(i,itemp)=ebz(i)
                  etxrep(i,itemp)=etx(i) !CB
                  etyrep(i,itemp)=ety(i)
                  etzrep(i,itemp)=etz(i)
               endif
            enddo
c^^^^^^^^^^^^record done ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 2222    continue

ccccccccccccccccc<RMSD>, <E> cccccccccccccccccccccccccc
         do i=1,N_rep
            energ_sum(i)=energ_sum(i)+E_rep(i)
            energ_sum2(i)=energ_sum2(i)+E_rep(i)*E_rep(i)
            N_sum(i)=N_sum(i)+1
         enddo

ccccccccccccccccccccc snapshots of E(1), E(N_rep) ccccccccccccc
         if(icycle.eq.icycle/1*1)then
            i_tr=i_tr+1
            do k=1,n_repf
               write(30+k,401)Lch,E_rep(k),i_tr,icycle
               do i=1,Lch
                  if(mv(i).gt.0)then
                     abx=xrep(i,k)*0.87
                     aby=yrep(i,k)*0.87
                     abz=zrep(i,k)*0.87
                  else
                     abx=exrep(i,k)*0.87
                     aby=eyrep(i,k)*0.87
                     abz=ezrep(i,k)*0.87
                  endif
                  write(30+k,402)abx,aby,abz
               enddo
            enddo
         endif
 401     format(i8,1x,f10.1,2i8)
 402     format(f10.3,1x,f10.3,1x,f10.3)

ccccccccccccccccc swap replicas cccccccccccccccccccccccccccccccc
         if(icycle.eq.icycle/2*2)then
            do i=1,N_rep-1,2 !swap odd replicas
               call swap_RS(i,i+1)
            enddo
         else
            do i=2,N_rep-1,2
               call swap_RS(i,i+1) !swap even replicas
            enddo
          endif
          mcycle=mcycle+1
 1111 continue
 901  continue
c--------------------------Main cycle ended here!!---------------
      
cccccccccccccccccccccccc Na/Nt cccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*) 'RS-->i_move   move   Na(i)  Nt(i)   Na(i)/Nt(i)'
      do i=2,27
         if(bNt(i).gt.1)then
            write(20,5004) i,mname(i),bNa(i),bNt(i),bNa(i)/bNt(i)
         else
            write(20,5004) i,mname(i),bNa(i),bNt(i)
         endif
      enddo
 5004 format(I4,A9,2f15.1,f11.6)
      
ccccccccccccccccccccccccccE_final, NSa/NSt ccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'RS-->------------ E_final, Na_swap/Nt_swap ---------'
      WRITE(20,*) 'i  T(i) final_E(i)  Nsa(i)  Nst(i)  Nsa(i)/Nst(i)'
      do i=1, N_rep
         if(bNSt(i).gt.1)then
            write(20,5005) i,aT_rep(i),E_rep(i),
     $           bNSa(i),bNSt(i),bNSa(i)/bNSt(i)
         else
            write(20,5005) i,aT_rep(i),E_rep(i),bNSa(i),bNSt(i)
         endif
      enddo
 5005 format(I4,f7.2,f10.1,2f15.1,f11.6)
      energy_tot_tmp=energy_tot()
      write(20,*)'E_final=',energy_tot_tmp

ccccccccccccccccccccccc<E>, NNa/NNt ccccccccccccccccccccccccccc
      WRITE(20,*)
      WRITE(20,*)'RS---------- <energy>, Na(i)/Nt(i) ----------------'
      write(20,*)'i_rep   T    <E>   <RMSD>   c   Na   Nt   Na/Nt'
      do i=1,N_rep
        energ_sum(i)=energ_sum(i)/N_sum(i)
        energ_sum2(i)=energ_sum2(i)/N_sum(i)
        cheat=(energ_sum2(i)-energ_sum(i)**2)/(aT_rep(i)**2)
        write(20,5006)i,aT_rep(i),energ_sum(i),armsd_sum(i)
     &       /(N_rmsd(i)+0.00001),
     &       cheat,bNNa(i),bNNt(i),bNNa(i)/(bNNt(i)+0.00001)
      enddo
 5006 format(I4,f7.2,f10.1,f8.3,f8.1,2f12.1,f11.6)
      write(20,*)'E_min=',E_min

      write(20,*)
      write(20,*)'ncycle_max=',ncycle
      write(20,*)'ncycle_real=',mcycle
      write(20,*)
      write(20,*)'hour_max=',hour_max
      write(20,*)'hour_real=',atime
      write(20,*)
      write(20,*)'ending time: ',fdate()
      
      write(*,*)
      write(*,*)'ncycle_max=',ncycle
      write(*,*)'ncycle_real=',mcycle
      write(*,*)
      write(*,*)'hour_max=',hour_max
      write(*,*)'hour_real=',atime
      write(*,*)
      write(*,*)'ending time: ',fdate()
      
      stop
      return
      end

ccccccccccccccc print out initial parameters cccccccccccccccccc
      subroutine write_parameter
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      parameter(nrep=100)       !number of replicas
      common/commonuse1/random,ncycle
      common/commonuse2/atemp1,atemp2,N_rep,phot
      common/lengths/Lch,Lch1,Lch2
      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/distres/er4,es3c
      common/resnumber/Ncom,Ndis,accur
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/looks/exc,exc1,exc2
      common/pair1/eh2,eh1b,eh1c
      common/ranzy/nozy
      common/iter/n_run

      write(20,*)'Protein real length..................',Lch
      write(20,*)
      write(20,*)'Ncycle...............................',ncycle
      write(20,*)
      write(20,*)'Number of runs.......................',n_run
      write(20,*)
      write(20,*)'number of N_rep......................',N_REP
      write(20,*)'N_swap...............................',phot
      write(20,*)'ncycle*N_swap........................',ncycle*phot
      write(20,*)'local moves each replica:',ncycle*phot*float(Lch)
      write(20,*)'total moves in MC..',ncycle*phot*float(Lch)*N_rep
      write(20,*) 
      write(20,*)'.......N_dist........................',Ndis
      write(20,*)'.......N_comb........................',Ncom
      write(20,*)'.......N_dist/length........',Ndis/float(Lch)
      write(20,*)'.......N_comb/length........',Ncom/float(Lch)
      write(20,*) 
      write(20,*)'maximum temperture...................',atemp2
      write(20,*)'minimum temperture...................',atemp1
      write(20,*)
      write(20,*)'Excluded volumn parameter=',exc,sqrt(exc)*.87
      write(20,*)
      write(20,*)'................er1........',er1
      write(20,*)'................er3........',er3
      write(20,*)'................er4........',er4
      write(20,*)'................er5........',er5
      write(20,*)'................er6........',er6
      write(20,*)'................er7........',er7
      write(20,*)'................eh1c.......',eh1c
      write(20,*)
      write(20,*)'initial random number seed...........',random
      write(20,*)'contact order........................',acorder
      write(20,*)
      write(20,*)'the first arandom number.............',aranzy(nozy)
      write(20,*) 
      write(20,*)'*******************************************'
      if(switch.gt.1)then
         write(20,*)'This is RS running of threading structure!'
      else
         write(20,*)'This is a normal simulation.'
      endif
      write(20,*)'*******************************************'
      write(20,*) 

c ^^^^^^^^^^ write parameter finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccc set movement percentage cccccccccccccccccccc
cccccccccccccccccc CPU time of each movement ccccccccccccccccc
c                CPU        acceptance rate
c     move2      58s           11.5%
c     move3s     63s            7.5%
c     move3d     55s            6.0%
c     move4s     48s            3.3%
c     move4d     59s            2.6%
c     move4p     48s            3.3%
c     move5s     49s            1.8%
c     move5d     56s            1.3%
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_move_retio !set movement percentage
      implicit integer(i-z)
      common/commonuse2/atemp1,atemp2,N_rep,phot

      common/moveretio1/h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,h7,hend
      common/moveretio2/bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6
      common/moveretio3/bh7a,bh7b,bhendn,bhendc

      hsum=h2+h3s+h3d+h4s+h4d+h5s+h5d+h6+hend
      bh2=h2/hsum
      bh3s=bh2+h3s/hsum
      bh3d=bh3s+h3d/hsum
      bh4s=bh3d+h4s/hsum
      bh4d=bh4s+h4d/hsum
      bh5s=bh4d+h5s/hsum
      bh5d=bh5s+h5d/hsum
      bh6=bh5d+h6/hsum
      bhendn=bh6+hend/2./hsum
      bhendc=bhendn+hend/2./hsum

c      write(*,1)h2,h3s,h3d,h4s,h4d,h5s,h5d,h6,hend
c      write(*,1)bh2,bh3s,bh3d,bh4s,bh4d,bh5s,bh5d,bh6,
c     $     bhendn,bhendc
c      write(*,1)bh2,bh3s-bh2,bh3d-bh3s,bh4s-bh3d,bh4d-bh4s,bh5s-bh4d,
c     $     bh5d-bh5s,bh6-bh5d,bhendn-bh6,bhendc-bhendn
c 1    format(13f6.3)
c      stop

c ^^^^^^^^^^^^^^^^^^^^^^^ set retio finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set current (x,y,z), ica, T, from itemp's replica
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_current
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/echain2/egx(ndim),egy(ndim),egz(ndim)
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim)
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim)

*     put the itemp'th replica conformation as current conformation------>
*     moveable points:
      do i=1,Lch
         x(i)=xrep(i,itemp)     !initial coordinate of itemp's replica
         y(i)=yrep(i,itemp)
         z(i)=zrep(i,itemp)
      enddo

      do i=1,Lch1
         j=i+1
         wx=x(j)-x(i)
         wy=y(j)-y(i)
         wz=z(j)-z(i)
         ica(i)=vector(wx,wy,wz) !identify order number of each bond-vector
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

c^^^^^^^^^^^ set current (x,y,z,ica) finished ^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set current (x,y,z), ica, T, from itemp's replica
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine set_current_RS
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep) !Ca
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep) !SG
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep) !cc
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep) !hb
      common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep) !CB

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !Ca
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)

      do i=1,Lch
         x(i)=xrep(i,itemp)
         y(i)=yrep(i,itemp)
         z(i)=zrep(i,itemp)
         ica(i)=icarep(i,itemp)
         if(mv(i).lt.0)then
            ex(i)=exrep(i,itemp) !Ca
            ey(i)=eyrep(i,itemp)
            ez(i)=ezrep(i,itemp)
            egx(i)=egxrep(i,itemp) !SG
            egy(i)=egyrep(i,itemp)
            egz(i)=egzrep(i,itemp)
            ecx(i)=ecxrep(i,itemp) !cc
            ecy(i)=ecyrep(i,itemp)
            ecz(i)=eczrep(i,itemp)
            ebx(i)=ebxrep(i,itemp) !hb
            eby(i)=ebyrep(i,itemp)
            ebz(i)=ebzrep(i,itemp)
            etx(i)=etxrep(i,itemp) !CB
            ety(i)=etyrep(i,itemp)
            etz(i)=etzrep(i,itemp)
         endif
      enddo
      ica(0)=ica(2)             !only useful for moveable pieces
      ica(Lch)=ica(Lch2)        !only useful for moveable pieces

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, translate (x,y,z) to calculate (amx,amy,amz)
c     2, calculate initial energy, 
c     3, prepare initial parameters:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initial_move
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/seqe/seq(ndim),sec(ndim)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev

      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/chainm/mv(ndim)

      dimension axx(ndim),ayy(ndim),azz(ndim)

*     check excluded volumn
      call get_vvv
*     move coordinates of the chain to its center of mass ---->
      call get_center

******* get template ready *********************
      

*     calculate (amx,amy,amz) -------------->
c      call get_axis

cccccccc prepare all the initial parameters of this replica ccccccccc
      icnto=0                   !need when ISTAT>0
      sumcto=0                  !need when ISTAT>0
      do i=1,Lch
         nop(i)=0               !number of parallel contacts
         noa(i)=0
         nom(i)=0
      enddo

      energ=EHB(1,Lch,1)+ESHORT(1,Lch,10) !initialize all below parameters
ccc   eprofo,eprofn: calculated each time;
ccc   icnt: always start from icnto; (need initial)
ccc   sumct: always start from sumcto; (need initial)
ccc   nop(): start from 0, or nopp; (need initial)
ccc   noa(): start from 0, or noaa; (need initial)
ccc   nom(): start from 0, or nomm; (need initial)
ccc   codevsum: only used when istat=-1, i.e. in Enew.
ccc   didevsum: only used when istat=-1, i.e. in Enew.
ccc   afs(): keep in store.

cccc  backup all the parameters, calculated in EHB() and ESHORT() cccccccccc
      eprofo=0.0
      do k=1,Lch
         is=seq(k)
         ia=noa(k)              !number of antiparallel contact-apirs
         ip=nop(k)              !number of parallel contact-apirs ^^
         im=nom(k)              !number of orgonal contact-apirs
         eprofo=eprofo+envir(ia,im,ip,is,3)
      enddo

      icnto=icnt                !backup of contact-order
      sumcto=sumct              !backup of contact-order
      do i=1,Lch
         nopp(i)=nop(i)         !backup of number of contact-pair
         noaa(i)=noa(i)         !backup of number of contact-pair
         nomm(i)=nom(i)         !backup of number of contact-pair
      enddo
      codevsum=conew            !backup  panelity for total deviation of comb
      didevsum=dinew            !backup  panelity for total deviation of dist
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c^^^^^^^^^^^^^^^^^ initialization finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     record the centroid
c     record EE(3),HH(3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_center
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/initialinput/switch,k_cycle,k_phot,N_ann

      common/chain0/ras(ndim),nfl
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/chainm/mv(ndim)
      common/seqe/seq(ndim),sec(ndim)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      
      common/center/cex,cey,cez
      common/eigen/AA(3,3),EE(3),HH(3,3)

*********************************************************
*     record centriod------------------------->
      cex=0
      cey=0
      cez=0
      if(switch.eq.1)then
         do i=1,Lch
            cex=cex+x(i)
            cey=cey+y(i)
            cez=cez+z(i)
         enddo
      else                      !!!!!!!!template based run
         do i=1,Lch
            if(mv(i).gt.0)then
               cex=cex+x(i)
               cey=cey+y(i)
               cez=cez+z(i)
            else
               cex=cex+ex(i)
               cey=cey+ey(i)
               cez=cez+ez(i)
            endif
         enddo
      endif
      cex=cex/float(Lch)
      cey=cey/float(Lch)
      cez=cez/float(Lch)
      
*********************************************************************
*     calculate axis ----------------------->
      AA(1,1)=0
      AA(1,2)=0
      AA(1,3)=0
      AA(2,2)=0
      AA(2,3)=0
      AA(3,3)=0
      do i=1,Lch
         if(mv(i).gt.0)then
            agxi=x(i)+GX(ica(i-1),ica(i),seq(i))-cex
            agyi=y(i)+GY(ica(i-1),ica(i),seq(i))-cey
            agzi=z(i)+GZ(ica(i-1),ica(i),seq(i))-cez
         else
            agxi=egx(i)-cex
            agyi=egy(i)-cey
            agzi=egz(i)-cez
         endif
         AA(1,1)=AA(1,1)+agxi*agxi
         AA(1,2)=AA(1,2)+agxi*agyi
         AA(1,3)=AA(1,3)+agxi*agzi
         AA(2,2)=AA(2,2)+agyi*agyi
         AA(2,3)=AA(2,3)+agyi*agzi
         AA(3,3)=AA(3,3)+agzi*agzi
      enddo
      AA(1,1)=AA(1,1)/float(Lch)
      AA(1,2)=AA(1,2)/float(Lch)
      AA(1,3)=AA(1,3)/float(Lch)
      AA(2,2)=AA(2,2)/float(Lch)
      AA(2,3)=AA(2,3)/float(Lch)
      AA(3,3)=AA(3,3)/float(Lch)
      AA(2,1)=AA(1,2)
      AA(3,1)=AA(1,3)
      AA(3,2)=AA(2,3)

c     E(1)=<x^2>=Ax^2/5; E(2)=<y^2>=Ay^2/5; E(3)=<z^2>=Az^2/5
      call eigenvalue           !get rotation matrix TT

c^^^^^^^^^^^^^^^^^ (amx,amy,amz) finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c check whether passage of [jj,kk] overlap with other parts of chain.
c look =.ture., without overlap
c look =.false., with overlap
c only C_alpha is checked.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION LOOK(jj,kk)
      IMPLICIT INTEGER(I-Z)
      LOGICAL LOOK
      parameter(ndim=1999)
      parameter(nvec=416)
      common/seqe/seq(ndim),sec(ndim) 	 	
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/looks/exc,exc1,exc2
      common/chainm/mv(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/excluded/vvv(ndim,ndim)

c     Ca-Ca r2>13;   r>sqrt(13)*0.87=3.14A
c     Ca-Ca r2>14;   r>sqrt(14)*0.87=3.26A
c     Ca-Ca r2>15;   r>sqrt(15)*0.87=3.37A
c     Ca-Ca r2>16;   r>sqrt(16)*0.87=3.48A
c     Ca-Ca r2>17;   r>sqrt(17)*0.87=3.59A
c     Ca-Ca r2>18;   r>sqrt(18)*0.87=3.69A
c     Ca-Ca r2>19;   r>sqrt(19)*0.87=3.79A
c     Ca-Ca r2>20;   r>sqrt(20)*0.87=3.89A
c     Ca-Ca r2>21;   r>sqrt(21)*0.87=3.99A
c     Ca-Ca r2>22;   r>sqrt(22)*0.87=4.08A
c     Ca-Ca r2>23;   r>sqrt(23)*0.87=4.17A *
c     Ca-Ca r2>24;   r>sqrt(24)*0.87=4.26A
c     Ca-Ca r2>25;   r>sqrt(25)*0.87=4.35A
c     Ca-Ca r2>26;   r>sqrt(26)*0.87=4.44A
c     Ca-Ca r2>27;   r>sqrt(27)*0.87=4.52A

      if(switch.eq.1)then       !normal lattice running
         do k=jj,kk
            i1=k-3
            i2=max(k+3,kk+1)
            do i=1,Lch
               if(i.le.i1.or.i.ge.i2)then
                  ir2=(x(k)-x(i))**2+(y(k)-y(i))**2+(z(k)-z(i))**2
                  if(ir2.lt.exc) then
                     LOOK=.FALSE.
                     RETURN
                  endif
               endif
            enddo
         enddo
      else                      !off-lattice
         do k=jj,kk
            if(mv(k).gt.0)then  !off-lattice
               axk=x(k)
               ayk=y(k)
               azk=z(k)
            else                !on-lattice
               axk=ex(k)
               ayk=ey(k)
               azk=ez(k)
            endif
            i1=k-3
            i2=max(k+3,kk+1)
            do i=1,Lch
               if(i.le.i1.or.i.ge.i2)then
                  if(vvv(i,k).gt.0)then
                     if(mv(i).gt.0)then
                        axi=x(i)
                        ayi=y(i)
                        azi=z(i)
                     else
                        axi=ex(i)
                        ayi=ey(i)
                        azi=ez(i)
                     endif
                     ar2=(axk-axi)**2+(ayk-ayi)**2+(azk-azi)**2
                     if(ar2.lt.exc) then
                        LOOK=.FALSE.
                        RETURN
                     endif
                  endif
               endif
            enddo
         enddo
      endif
      LOOK=.TRUE.

c ^^^^^^^^^^ Look finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, pairwise interaction;
c     2, H-bond;
c     3, Electric;
c     4, contact number;
c     5, contact order.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION EHB(jjjj,kkkk,ISTAT)

      IMPLICIT INTEGER(I-Z)
      parameter(api=3.1415926)
      parameter(ndim=1999)
                parameter(nvec=416)
      common/seqe/seq(ndim),sec(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      COMMON/pair/ apa(ndim,ndim),app(ndim,ndim),apm(ndim,ndim)
      COMMON/size/ arla(0:19,0:19),arlm(0:19,0:19),arlp(0:19,0:19)
      COMMON/sizea/ ala(0:19,0:19),alm(0:19,0:19),alp(0:19,0:19)
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/icgg/ icg(ndim),EH6  		
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3 
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)	
      common/pair1/eh2,eh1b,eh1c
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,EHB5,EHB6
      common/distres/er4,es3c
      common/chainm/mv(ndim)
      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh
      common/ehbenergy1/EHB5a,EHB5b
      common/par/apar(ndim,ndim)
      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
      common/concuttu/npot1,concutU2(0:19,0:19),concutU(0:19,0:19),aw1
      common/concuttu4/npot4,aw4,Cr20
      common/concuttu1/dq1a(0:19,0:19),dq1b(0:19,0:19)
      common/concuttu2/dq1c(0:19,0:19),dq1d(0:19,0:19)
      
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !Ca
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/temp/npp

      EHB1=0                    !+1/r of Ca-SC
      EHB1a=0                   !+1/r for non-parallel contact of Ca-Ca
      EHB1b=0                   !excluded volume of SC-SC 
      EHB1c=0                   !pair potential of SC-SC, paor1.dat+ par.dat

      EHB2=0                    !quarsi3 for SC-SC, pair3.dat+quarsi3.comm
      EHB3=0                    !enhance good-piece contacts
      EHB4=0                    !-1/r for parallel contact of Ca-Ca

      EHB5a=0                    !H-bond energy for alpha-type
      EHB5b=0                    !H-bond energy for beta-type

      if(ISTAT.gt.0)THEN        !Enew
         ICNT=ICNTO
         SUMCT=SUMCTO
      ENDIF

c     coupling of secondary structure with pairwise interactions
c     included - thus expanded range
      jj=jjjj-1
      if(jj.lt.1)jj=1
      kk=kkkk+1
      if(kk.gt.Lch)kk=Lch
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*****************************************************************
********** circle for movement involved window ******************
*****************************************************************
      DO 1002 k=jj,kk
         kseq=seq(k)
         if(mv(k).gt.0)then     !moveable point
            axk=x(k)            !Ca
            ayk=y(k)
            azk=z(k)
            nv1=ica(k-1)
            nv2=ica(k)
            agxk=axk+GX(nv1,nv2,kseq) !Cg_k
            agyk=ayk+GY(nv1,nv2,kseq)
            agzk=azk+GZ(nv1,nv2,kseq)
            bxk=HBX(nv1,nv2)    !H-bond direction
            byk=HBY(nv1,nv2)
            bzk=HBZ(nv1,nv2)
            cxk=CAX(nv1,nv2)    !vertical vector
            cyk=CAY(nv1,nv2)
            czk=CAZ(nv1,nv2)	
         else
            axk=ex(k)           !Ca
            ayk=ey(k)
            azk=ez(k)
            agxk=egx(k)         !SG
            agyk=egy(k)
            agzk=egz(k)
            bxk=ebx(k)          !Hb
            byk=eby(k)
            bzk=ebz(k)
            cxk=ecx(k)          !cc
            cyk=ecy(k)
            czk=ecz(k)
         endif
         km2=k-2
         kp2=k+2
         if(km2.ge.1.and.kp2.le.Lch)then
            xxx=nint((aax(km2)-aax(kp2))**2+
     $           (aay(km2)-aay(kp2))**2+(aaz(km2)-aaz(kp2))**2) !real dist
            if(xxx.gt.500)xxx=500
            ek5=acops(km2,jbin(xxx)) !k2-residue, jbin. es<0.
         else
            ek5=0
         endif
*****************************************************************
********************** start i-cicle ****************************
*****************************************************************
         do 1001 i=1,Lch
            iend=max(k+1,kk)
            if(i.ge.k-1.and.i.le.iend)goto 1001 !to avoid repeat
            if(mv(i).gt.0)then
               axi=x(i)
               ayi=y(i)
               azi=z(i)
            else
               axi=ex(i)
               ayi=ey(i)
               azi=ez(i)
            endif
            axki=(axk-axi)      !r_k-r_i
            ayki=(ayk-ayi)
            azki=(azk-azi)
            Cr2=axki**2+ayki**2+azki**2+0.000001 !(r_k-r_i)**2
c            if(Cr2.lt.120)then  !120-->9.53, there is pairwise interaction
            if(Cr2.lt.Cr20)then !120-->9.53, there is pairwise interaction
               idist=iabs(i-k)
               iseq=seq(i)
               if(mv(i).gt.0)then
                  nv1=ica(i-1)
                  nv2=ica(i)
                  agxi=axi+GX(nv1,nv2,iseq) !Cg_k
                  agyi=ayi+GY(nv1,nv2,iseq)
                  agzi=azi+GZ(nv1,nv2,iseq)
                  bxi=HBX(nv1,nv2) !H-bond direction
                  byi=HBY(nv1,nv2)
                  bzi=HBZ(nv1,nv2)
                  cxi=CAX(nv1,nv2) !outside 
                  cyi=CAY(nv1,nv2)
                  czi=CAZ(nv1,nv2)	
               else
                  agxi=egx(i)   !SG
                  agyi=egy(i)
                  agzi=egz(i)
                  bxi=ebx(i)    !H-bond
                  byi=eby(i)
                  bzi=ebz(i)
                  cxi=ecx(i)    !bisector vector
                  cyi=ecy(i)
                  czi=ecz(i)
               endif

c     1/r excluded for Ca-SC pair ----------------->
               if(Cr2.lt.120)then  !120-->9.53, there is pairwise interaction
               if(kseq.gt.0) then !not GLY, we have SG
                  aks=(agxk-axi)**2+(agyk-ayi)**2+(agzk-azi)**2+.0001 !SG_kCa_i
                  if(aks.lt.36) then !36-->5.22A
                     if(aks.lt.13.0) aks=13.0 !13-->3.17A
                     EHB1=EHB1+((13.0/aks)-0.5)
                  endif
               endif
               if(iseq.gt.0) then
                  aks=(axk-agxi)**2+(ayk-agyi)**2+(azk-agzi)**2+.0001 !Ca_kSG_i
                  if(aks.lt.36.0) then
                     if(aks.lt.13.0) aks=13.0
                     EHB1=EHB1+((13.0/aks)-0.5)
                  endif
               endif
               endif
               
c     pair-wise potential of SC-SC from 'pair1.dat' and 'par.dat' ---------->
               Gr2=(agxi-agxk)**2+(agyi-agyk)**2+(agzi-agzk)**2 !SG_i,SG_k
c               if(Gr2.lt.concut2(iseq,kseq))then
c                  EHB1c=EHB1c+apar(i,k)
c               endif
cccc  
               if(Gr2.lt.concut2(iseq,kseq))then
                  EHB1c=EHB1c+apar(i,k) !apar=Energy, i.e., the lower the better
               endif
c     npp=npp+1
c     if(npp.le.4000)then
c     write(93,*)npp,npot4,sqrt(Gr2)*0.87,ag,i,k,k-i,jjjj,kkkk
c     endif
cccc  
               
c     quarsi3 for SC-SC pair, 1/r for Ca-Ca ----------------->
               cc=cxi*cxk+cyi*cyk+czi*czk !c_i*c_k
               IF(cc.gt.0.5)THEN !c_i//c_k
                  if(Cr2.lt.60)EHB4=EHB4-(30/(max(30.,Cr2))-0.5) !6.74A
                  if(Gr2.lt.alp(iseq,kseq))then
                     EHB3=EHB3-ek5*ei5(i,idist)
                     NOP(k)=NOP(k)+ISTAT
                     NOP(i)=NOP(i)+ISTAT
                     ICNT=ICNT+istat*idist
                     SUMCT=SUMCT+ISTAT
                     if(Gr2.gt.arlp(iseq,kseq))THEN
                        EHB2=EHB2+app(i,k) !quarsi3
                     else
                        EHB1b=EHB1b+1 !excluded volume
                     endif
                  endif
               ELSE
                  if(Cr2.lt.33)EHB1a=EHB1a+(16.0/Cr2-0.5) !5A
                  IF(cc.lt.-0.5) THEN !antiparallel pair-wise of (i,k)
                     if(Gr2.lt.ala(iseq,kseq))then
                        EHB3=EHB3-ek5*ei5(i,idist)
                        NOA(k)=NOA(k)+ISTAT
                        NOA(i)=NOA(i)+ISTAT
                        ICNT=ICNT+istat*idist
                        SUMCT=SUMCT+ISTAT
                        if(Gr2.gt.arla(iseq,kseq))THEN
                           EHB2=EHB2+apa(i,k)
                        else
                           EHB1b=EHB1b+1 !excluded volume
                        endif
                     endif
                  ELSE          !neither parallel nor antiparallel
                     if(Gr2.lt.alm(iseq,kseq)) then
                        EHB3=EHB3-ek5*ei5(i,idist)
                        NOM(k)=NOM(k)+ISTAT
                        NOM(i)=NOM(i)+ISTAT
                        ICNT=ICNT+istat*idist !distance of pairs
                        SUMCT=SUMCT+ISTAT !number of contact pairs
                        if(Gr2.gt.arlm(iseq,kseq))THEN
                           EHB2=EHB2+apm(i,k)
                        else
                           EHB1b=EHB1b+1 !excluded volume
                        endif
                     endif
                  ENDIF
               ENDIF
c^^^^^^^^^^^^^ interaction of pair-wise finished ^^^^^^^^^^^^^^^^^

ccccccccccccccccccccccccc Hydrogen-bond cccccccccccccccccccccccccc
***** 1, cc; 2, bb; 3, distance; 4, vij; 5, sec(i).
c               if(idist.eq.3)then
c                  alpha-H-bond
c               elseif(idist.lt.20)then
c                  anti-parallel-sheet-H-bond
c               else
c                  if(bb<0)then
c                     antiparallel sheet-H-bond
c                  else
c                     parallel sheet-H-bond
c                  endif
c               endif
*HHHHH***************** alpha ******************************************
***   bb   bb_an  cc    cc_an vv   vv_an hbl r2i  r2k  Cr2  Cr2m  Hdis
***   0.82  34.9  0.42  65.3  0.52  58.8 5.7 3.29 3.53 35.9 31.62 4.19 (alpha)
***   0.81  35.3  0.38  67.3  0.48  61.0 5.7 3.55 3.39 35.5 30.45 4.10 (alpha1)
***   0.85  31.9  0.50  59.5  0.58  54.0 5.7 3.03 2.92 34.1 30.29 4.13 (alpha2)
************************************************************************
*********************** beta *******************************************
***   bb    bb_an   cc  cc_an vv12 vv12_a vv12 vv12_a hbl r2i r2k Cr2 Cr2m
***   -0.81 144.8 0.79  36.7 -0.94 160.6 -0.95 162.4 5.3 8.38 1.37 30 30 (r)
***   -0.84 149.1 0.83  31.8 -0.92 161.1 -0.94 162.7 5.3 3.91 4.43 33 31 (r1)
***   -0.64 135.0 0.46  61.8 -0.75 139.6 -0.74 137.6 5.1 7.81 6.55 33 31 (r2)
***    0.98   8.6 0.92  22.2  0.92  22.5  0.93  18.8 5.3 2.49 0.99 31 30 (p)
************************************************************************
               if(Cr2.lt.120)then  !120-->9.53, there is pairwise interaction
               if(i.eq.1.or.i.eq.Lch.or.k.eq.1.or.k.eq.Lch)goto 1003
               if(idist.eq.3)then !!!!!!!!!!!!!!!!!!!!!!!!!!!
***   alpha-helix:
                  if(sec(k).ne.4.and.sec(i).ne.4)then
                  if(Cr2.lt.Cr2a.and.cc.gt.acut_cc)then
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk)
                  if(bb.gt.acut_bb)then
                  av11=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                  if(av11.gt.acut_vv)then
                  av22=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                  if(av22.gt.acut_vv)then
***   ->
                     fact=(1-abs(cc-0.4))*(1-abs(bb-0.815))
                     EHB5a=EHB5a+energyHBa(i,k,fact,
     $                    axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                  endif
                  endif
                  endif
                  endif
                  endif
               elseif(idist.gt.4.and.idist.lt.20)then !!!!!!!!!!!!!!!!!!!!!!!
***   antiparallel-sheet
                  if(sec(k).ne.2.and.sec(i).ne.2)then
                  if(Cr2.lt.Cr2b.and.cc.gt.bcut_cc)then
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk)
                  if(bb.lt.-bcut_bb)then !antiparallel
                  av12=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                  if(av12.lt.-bcut_vv)then
                  av21=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                  if(av21.lt.-bcut_vv)then
***   ->
                     fact=abs(bb)*cc !bb->-1,cc->1
                     EHB5b=EHB5b+energyHBb(i,k,fact,
     $                    axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                  endif
                  endif
                  endif
                  endif
                  endif
               elseif(idist.ge.20)then !!!!!!!!!!!!!!!!!!!!!!!!
                  if(sec(k).ne.2.and.sec(i).ne.2)then
                  if(Cr2.lt.Cr2b.and.cc.gt.bcut_cc)then
                  bb=bxi*bxk+byi*byk+bzi*bzk !cos(Hbi.Hbk) 
                  if(bb.lt.-bcut_bb)then !antiparallel
***   antiparallel-sheet:
                     av12=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                     if(av12.lt.-bcut_vv)then
                     av21=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                       aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                     if(av21.lt.-bcut_vv)then
***   ->
                        fact=abs(bb)*cc !bb->1,cc->1
                        EHB5b=EHB5b+energyHBb(i,k,fact,
     $                       axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                     endif
                     endif
                  elseif(bb.gt.bcut_bb)then !bb>0, parallel H-bond
***   parallel-sheet:
                     av11=avv(aax(i-1),aay(i-1),aaz(i-1),axi,ayi,azi,
     $                    aax(k-1),aay(k-1),aaz(k-1),axk,ayk,azk)
                     if(av11.gt.bcut_vv)then
                     av22=avv(axi,ayi,azi,aax(i+1),aay(i+1),aaz(i+1),
     $                       axk,ayk,azk,aax(k+1),aay(k+1),aaz(k+1))
                     if(av22.gt.bcut_vv)then
***   ->
                        fact=bb*cc !bb->1,cc->1
                        EHB5b=EHB5b+energyHBb(i,k,fact,
     $                       axki,ayki,azki,bxi,byi,bzi,bxk,byk,bzk)
                     endif
                     endif
                  endif
                  endif
                  endif
               endif            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 1003          continue
               endif            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c^^^^^^^^^^^^^H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            endif               !Ca-Ca<120
 1001    continue               !i -> [1,Lch]
 1002 continue                  !k -> [jj,kk]

c     ICNT/SUMCT: average distance of residues in each contact-pairs.
c     ICNT/SUMCT are calculated when call Enew
c     ICNTO/SUMCTO are not changed when call Eold
      if(istat.lt.0) then       !Eold
         b=abs(ICNT/(0.000001+float(SUMCT))-acorder)
         a=abs(ICNTO/(0.000001+float(SUMCTO))-acorder)
         d=abs(float(SUMCT)-contt) !deviation of contact number on new conform
         c=abs(float(SUMCTO)-contt) !deviation of contact number on new confor
         dord=en2*(b-a)+en3*(d-c) !not included in EHB
      endif
      
      EHB
     $     =eh1*EHB1            !+1/r of Ca-SC
     $     +eh1a*EHB1a          !+1/r for non-parallel of Ca-Ca
     $     +eh1b*EHB1b          !excluded volumn of SC-SC
     $     +eh1c*EHB1c          !pair-wise potential of SC-SC
     $     +eh2*EHB2            !quarsi3 for SC-SC
     $     +eh3*EHB3            !enhance good piece
     $     +eh4*EHB4            !-1/r for parallel contact of Ca-Ca
     $     +eh5a*EHB5a          !H-bond energy (alpha)
     $     +eh5b*EHB5b          !H-bond energy (beta)

c      write(*,*)EHB,eh1c,EHB1c,jjjj,kkkk

c ^^^^^^^^^^ EHB finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     H-bond for alpha-helix
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function energyHBa(i,k,a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk)
      implicit integer(i-z)
      parameter(ndim=1999)      !maximum length of chain-length
      parameter(nvec=416)       !number of vectors
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh

c     because helix is right-hand, bi=v(i-1)(x)v(i) is always
c     same direction as v(i---k), if (k>i);
c     reverse direction as v(i---k), if (k<i);

      energyHBa=0

      bxk=5.7*dxk
      byk=5.7*dyk
      bzk=5.7*dzk
      bxi=5.7*dxi
      byi=5.7*dyi
      bzi=5.7*dzi
      if(k.gt.i)then            !bki,bxk,axki are in same directory
         br=(bxk-axki)**2+(byk-ayki)**2+(bzk-azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2)) !br->3.4
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
         br=(bxi-axki)**2+(byi-ayki)**2+(bzi-azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2))
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
      else                      !bki,bxk, and axki are in reverse directory
         br=(bxk+axki)**2+(byk+ayki)**2+(bzk+azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2))
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
         br=(bxi+axki)**2+(byi+ayki)**2+(bzi+azki)**2
         if(br.lt.acut_hh)then
            b=1/(1+abs(br-3.2))
            energyHBa=energyHBa-EHBIJ(i,k)*a*b
         endif
      endif

c^^^^^^^^^^^^ H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     H-bond for beta-sheet
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function energyHBb(i,k,a,axki,ayki,azki,dxi,dyi,dzi,dxk,dyk,dzk)
      implicit integer(i-z)
      parameter(ndim=1999)      !maximum length of chain-length
      parameter(nvec=416)       !number of vectors
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      common/hba/eh5a,Cr2a,acut_bb,acut_cc,acut_vv,acut_hh
      common/hbb/eh5b,Cr2b,bcut_bb,bcut_cc,bcut_vv,bcut_hh

      energyHBb=0

      bxk=5.3*dxk
      byk=5.3*dyk
      bzk=5.3*dzk
      bxi=5.3*dxi
      byi=5.3*dyi
      bzi=5.3*dzi
      br=(bxk-axki)**2+(byk-ayki)**2+(bzk-azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)           !br->0
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif
      br=(bxk+axki)**2+(byk+ayki)**2+(bzk+azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif
      br=(bxi-axki)**2+(byi-ayki)**2+(bzi-azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif
      br=(bxi+axki)**2+(byi+ayki)**2+(bzi+azki)**2
      if(br.lt.bcut_hh) then
         b=2./(2.+br)
         energyHBb=energyHBb-EHBIJ(i,k)*a*b
      endif

c^^^^^^^^^^^^ H-bond energy finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate v(12).v(34)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function avv(ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3,ax4,ay4,az4)

      ax12=ax2-ax1
      ay12=ay2-ay1
      az12=az2-az1
      ax34=ax4-ax3
      ay34=ay4-ay3
      az34=az4-az3
      avv=(ax12*ax34+ay12*ay34+az12*az34)/(0.000001+
     $     sqrt((ax12**2+ay12**2+az12**2)*(ax34**2+ay34**2+az34**2)))

c^^^^^^^^^^^^^ v.v finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate ei5=acops(i,jbin(r15))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ei5(i,idist)
      implicit integer(i-z)
      parameter(ndim=1999)
      common/chainm/mv(ndim)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim)
      common/lengths/Lch,Lch1,Lch2
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)

      ei5=0
      if(idist.gt.5)then
         im2=i-2
         ip2=i+2
         if(im2.ge.1.and.ip2.le.Lch)then
            if(mv(im2).gt.0)then
               axm=x(im2)
               aym=y(im2)
               azm=z(im2)
            else
               axm=ex(im2)
               aym=ey(im2)
               azm=ez(im2)
            endif
            if(mv(ip2).gt.0)then
               axp=x(ip2)
               ayp=y(ip2)
               azp=z(ip2)
            else
               axp=ex(ip2)
               ayp=ey(ip2)
               azp=ez(ip2)
            endif
            xxx=nint((axm-axp)**2+(aym-ayp)**2+(azm-azp)**2)
            if(xxx.gt.500)xxx=500
            ei5=acops(im2,jbin(xxx)) !i2-residue, jbin. es<0.
         endif
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, centrosymmetric energy of (x,y,z) and Cg. (centro.comm).
c     2, bury potential of SC.
c     3, distmap and contact restrains.
c     4, bias to (prodicted) protein-like structure, panality on crumpling.
c     5, E13,E15,E15.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION ESHORT(iiii,jjjj,ISTAT)
      IMPLICIT INTEGER(I-Z)
      parameter(ndim=1999)	
      parameter(nrep=100)       !maximum number of replicas
      parameter(nvec=416)
      parameter(api=3.1415926)
      common/lengths/Lch,Lch1,Lch2
      common/seqe/seq(ndim),sec(ndim) 	  
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/hb/hbx(nvec,nvec),hby(nvec,nvec),hbz(nvec,nvec)
      common/bisec/cax(nvec,nvec),cay(nvec,nvec),caz(nvec,nvec)
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      COMMON/short/ IBIN(-300:300),asr(ndim,-12:12)
      COMMON/short1/ JBIN(0:500),bsr(ndim,16),acops(ndim,16)
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/fr/frga(ndim),frgb(ndim)
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/three/angle(nvec,nvec)
      common/sg/gx(nvec,nvec,0:19),gy(nvec,nvec,0:19),gz(nvec,nvec,0:19)
      common/cb/hx(nvec,nvec,0:19),hy(nvec,nvec,0:19),hz(nvec,nvec,0:19)
      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      COMMON/RCN/Mdis(ndim),kdis(ndim,100),dist(ndim,100),dev(ndim,100)
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/msichores/msicho
      common/distres/er4,es3c
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      common/eshortenergy1/ESHORT1,ESHORT2,ESHORT3,ESHORT4,ESHORT11
      common/eshortenergy2/ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,ESHORT5c
      common/eshortenergy3/ESHORT6,ESHORT7,ESHORT8,ESHORT9,ESHORT10
      common/eshortenergy4/ESHORT12
      common/hopp/eonehw(0:19)
      common/concutt/concut(0:19,0:19),concut2(0:19,0:19),concut_sc
      common/chainm/mv(ndim)

      common/hom1/asr2(ndim,4),asr3(ndim,4),asr4(ndim,14),asr5(ndim,8)
      common/hom2/ibb2(0:999),ibb3(0:999),ibb4(-999:999),ibb5(0:999)
      common/CAcontact/McomCA(ndim),KcomCA(ndim,100),aweigCA(ndim,ndim)
      common/CA8/McomCA8(ndim),KcomCA8(ndim,100),aweigCA8(ndim,ndim)
      common/CA8a/npot3,e_CA_cut,e_CA_cut2,e_CA_cutU,e_CA_cutU2,aw3
      common/CA8b/eq2a,eq2b,eq2c,eq2d
      common/CB6a/npot5,aw5,f_CB_cut,f_CB_cut2,f_CB_cutU,f_CB_cutU2
      common/CB6b/fq2a,fq2b,fq2c,fq2d
      common/CB6c/McomCB6(ndim),KcomCB6(ndim,100),aweigCB6(ndim,ndim)
      common/CB8a/npot6,aw6,g_CB_cut,g_CB_cut2,g_CB_cutU,g_CB_cutU2
      common/CB8b/gq2a,gq2b,gq2c,gq2d
      common/CB8c/McomCB8(ndim),KcomCB8(ndim,100),aweigCB8(ndim,ndim)
      
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !Ca
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/zscore/izscore
      common/freg/aweig(ndim,ndim)
      common/CAc1/npot2,d_CA_cut,d_CA_cut2,d_CA_cutU,d_CA_cutU2,aw2
      common/CAcontact1a/dq2a,dq2b,dq2c,dq2d
      COMMON/longdist/MdisL(ndim),kdisL(ndim,500),distL(ndim,500)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/rmsdder/i_chunk(ndim),ex0(ndim),ey0(ndim),ez0(ndim)
      common/E_defo/i_E_defo
      common/initialinput/switch,k_cycle,k_phot,N_ann
      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)
      common/expose/mp(20,ndim),area(ndim)
      common/res2/er14,er15,er16,er17

      common/sidechain/bgx(ndim),bgy(ndim),bgz(ndim)
      dimension cgx(ndim),cgy(ndim),cgz(ndim)

      common/center/cex,cey,cez
      common/eigen/AA(3,3),EE(3),HH(3,3)
      
      common/stick1/nstick,astick,nrmsd,ermsd
      common/stick2/iq(ndim,nrep)
      common/stick3/ax00(ndim,nrep),ay00(ndim,nrep),az00(ndim,nrep)
      common/stick6/bx00(ndim),by00(ndim),bz00(ndim),armsd_min0
      common/stick7/itemp0,icycle,icycle0
      common/trackn/n_tem(100)

      common/temperature/itemp,atemp
      
      common/svm1/Mcon(3,ndim),Kcon(3,ndim,300),awei(3,ndim,ndim),fw
      common/svm2/acc_cut,dist_svm2(15),nk_svm(3),ik_svm(3,15)
      common/svm3/er21,er22,er23
      common/svm6/npot,dist_svm(15),dist_svmU(15),dist_svmU2(15)
      common/concuttu/npot1,concutU2(0:19,0:19),concutU(0:19,0:19),aw1
      common/concuttu1/dq1a(0:19,0:19),dq1b(0:19,0:19)
      common/concuttu2/dq1c(0:19,0:19),dq1d(0:19,0:19)
      common/temp/npp
      common/dwell/dwell,ha(15),hb(15),hc(15),hd(15)
      
      dimension aij2_d(ndim)

cccc   RMSD:
      double precision r_1(3,ndim),r_2(3,ndim),r_3(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /ndim*1.0/
ccc   

c      ESHORT=0.0
      ESHORT2=0                 !bury potential for SC

      ESHORT3=0  !distance restrain from both threading (for C_a)
      ESHORT4=0  !contact restrain from threading (for SC)
      ESHORT4a=0 !deviation of contact restrain

      ESHORT5=0  !bias2,3: v(i)-v(i+4) anti/parallel; c(i)-c(i+2) anit/paralel
      ESHORT5a=0 !crumpling
      ESHORT5b=0 !bias4 to predicted alpha/beta structure.
      ESHORT5c=0 !bias1 to possible alpha/beta structure. 

      ESHORT6=0  !correlation of E13 of Ca
      ESHORT7=0  !correlation of E14, from both common and 2th specific data
      ESHORT8=0  !correlation of E15, from both common and 2th specific data

      ESHORT9=0  !contact restraints of combCA.dat
      ESHORT10=0 !Long-range distance restraints of CA

      ESHORT11=0 !RMSD deviation
      ESHORT12=0 !distance from template

      ESHORT13=0 !rmsd to template
c      ESHORT14=0                !contact restraints of comb8CA.dat
      
c      ESHORT17=0                !deviation to comb8CA.dat
      ESHORT21=0                !CA_SVM
      ESHORT22=0                !CB_SVM
      ESHORT23=0                !SG_SVM

      if(istat.lt.0)THEN       !Eold 
         diold=0.0
         coold=0.0
      else                      !Enew
         conew=0.0   !penality for incorrect contact (from contact restrain)
         dinew=0.0   !penality for incorrect distance (from distant restain)
      endif

c1ccccccccccccccccccccccccccccccccccccc
c     1, bury ellipsoid
c     2, 'dist.dat'
c     3, 'distL.dat'
cccccccccccccccccccccccccccccccccccccc
      do 1 i=iiii,jjjj
         iseq=seq(i)
         if(mv(i).gt.0)then
            axi=x(i)
            ayi=y(i)
            azi=z(i)
            nv1=ica(i-1)
            nv2=ica(i)
            agxi=axi+GX(nv1,nv2,iseq) !Cg_k
            agyi=ayi+GY(nv1,nv2,iseq)
            agzi=azi+GZ(nv1,nv2,iseq)
         else
            axi=ex(i)           !Ca
            ayi=ey(i)
            azi=ez(i)
            agxi=egx(i)         !SG
            agyi=egy(i)
            agzi=egz(i)
         endif
***** Bury energy of SG --------------------->
         gxn=HH(1,1)*(agxi-cex)+HH(1,2)*(agyi-cey)+HH(1,3)*(agzi-cez)
         gyn=HH(2,1)*(agxi-cex)+HH(2,2)*(agyi-cey)+HH(2,3)*(agzi-cez)
         gzn=HH(3,1)*(agxi-cex)+HH(3,2)*(agyi-cey)+HH(3,3)*(agzi-cez)
         fff=gxn**2/EE(1)+gyn**2/EE(2)+gzn**2/EE(3) !=5 for ellipsoid surphase
         if(fff.lt.2.5)then
           aaa=fff-2.5
            if(aaa.lt.-1)aaa=-1
            ESHORT2=ESHORT2-aaa*area(i) !area(i)=-1, buried; +1, exposed
         endif
c^^^^^^^^^^^^^^^^^^ bury energy finished ^^^^^^^^^^^^^^^^^^

c     deepth factor: afs=1, r<r0; [0.5,1], r0<r<2r0; 0.5, r>2r0
         ar=sqrt((axi-cex)**2+(ayi-cey)**2+(azi-cez)**2+0.01) !r
         afs(i)=acrit/ar
         if(afs(i).gt.1)afs(i)=1
         if(afs(i).lt.0.5)afs(i)=0.5
         
****  restrain from 'dist.dat'---------------------->
         N_dis=Mdis(i)          !number of distance restraints on 'i'
         do k=1,N_dis
            j=kdis(i,k)         !i-j restraints
            if(j.lt.i.or.j.gt.jjjj) then !to avoid repeat
               dij2=(axi-aax(j))**2+(ayi-aay(j))**2+(azi-aaz(j))**2
               dij=sqrt(dij2)
               err=abs(dij-dist(i,k)) !dist: predicted dis
               if(err.gt.dev(i,k)) then !dev: deviation for arca
                  ESHORT3=ESHORT3+1
               endif
            endif
         enddo
c^^^^^^^^^^^^^^^^^ distant restrain finished ^^^^^^^^^^^^^^^^^^^^^^^^^

****  long-range dist-restrain from 'distL.dat'---------------------->
         N_disL=MdisL(i)        !number of long-range dist restraints on 'i'
         do k=1,N_disL
            j=kdisL(i,k)        !i-j restraints
            if(j.lt.i.or.j.gt.jjjj) then !to avoid repeat
               dij2=(axi-aax(j))**2+(ayi-aay(j))**2+(azi-aaz(j))**2
               dij=sqrt(dij2)
               err=abs(dij-distL(i,k))
               if(err.lt.1)err=1
               ESHORT10=ESHORT10-1/err
            endif
         enddo
c^^^^^^^^^^^^^^^^^distant restrain finished ^^^^^^^^^^^^^^^^^^^^^^^^^
 1    continue
      
c2ccccccccccccccccccccccccccccccccccccc
c     4, 'comb.dat'
c     5, 'svmseqSG.dat' (not used)
cccccccccccccccccccccccccccccccccccccc
ccc   Contact restrain from 'comb.dat' -------------------------------->
       do 11 i=iiii,jjjj
          if(Mcom(i).ge.1.or.Mcon(3,i).ge.1)then
             iseq=seq(i)
             if(mv(i).gt.0)then
                agxi=x(i)+GX(ica(i-1),ica(i),iseq) !SG
                agyi=y(i)+GY(ica(i-1),ica(i),iseq) !SG
                agzi=z(i)+GZ(ica(i-1),ica(i),iseq) !SG
             else
                agxi=egx(i)     !SG
                agyi=egy(i)
                agzi=egz(i)
             endif
             
ccc   'comb.dat':
             do k=1,Mcom(i)
                j=Kcom(i,k)     !k'th contact with i
                if(j.lt.i.OR.j.gt.jjjj)then
                   jseq=seq(j)
                   if(mv(j).gt.0)then
                      agxj=x(j)+gx(ica(j-1),ica(j),jseq)
                      agyj=y(j)+gy(ica(j-1),ica(j),jseq)
                      agzj=z(j)+gz(ica(j-1),ica(j),jseq)
                   else
                      agxj=egx(j)
                      agyj=egy(j)
                      agzj=egz(j)
                   endif
                   aij2=(agxi-agxj)**2+(agyi-agyj)**2+(agzi-agzj)**2
                   if(aij2.gt.concut2(iseq,jseq))then
                      if(istat.lt.0)then
                         coold=coold+sqrt(aij2)-concut(iseq,jseq) !penality
                      else
                         conew=conew+sqrt(aij2)-concut(iseq,jseq) !panelity
                      endif
                   endif
cccc  
                   if(aij2.lt.concut2(iseq,jseq))then
                      ESHORT4=ESHORT4-aweig(i,j)
                   elseif(aij2.lt.concutU2(iseq,jseq))then
                      f1=(1-sin((sqrt(aij2)-dq1a(iseq,jseq))/
     &                     dq1b(iseq,jseq)*api))/2 !0<f1<1
                      ESHORT4=ESHORT4-aweig(i,j)*f1
                   elseif(aij2.lt.(80/0.87)**2)then
                      f1=(1+sin((sqrt(aij2)-dq1c(iseq,jseq))/
     &                     dq1d(iseq,jseq)*api))/2 !0<f1<1
                      ESHORT4=ESHORT4+aweig(i,j)*f1
                   else
                      ESHORT4=ESHORT4+aweig(i,j)
                   endif
c                   npp=npp+1
c                   if(npp.le.2000)then
c                      write(93,*)npp,npot,sqrt(aij2)*0.87,ag
c                   endif
cccc  
                endif
             enddo
             
ccc   'svmseqSG8.dat':
             goto 421           !without sequence-based SG contacts:
             do k=1,Mcon(3,i)   !number of SG-contacts on i
                j=Kcon(3,i,k)   !j is k'th SG-contact with i
                if(j.lt.i.OR.j.gt.jjjj)then
                   jseq=seq(j)
                   if(mv(j).gt.0)then
                      agxj=x(j)+gx(ica(j-1),ica(j),jseq)
                      agyj=y(j)+gy(ica(j-1),ica(j),jseq)
                      agzj=z(j)+gz(ica(j-1),ica(j),jseq)
                   else
                      agxj=egx(j)
                      agyj=egy(j)
                      agzj=egz(j)
                   endif
                   aij2=(agxi-agxj)**2+(agyi-agyj)**2+(agzi-agzj)**2
ccc   
                   if(aij2.lt.dist_svm2(3))then
                      ESHORT23=ESHORT23-awei(3,i,j)
                      mk=1
                      ag=-1
                   elseif(aij2.lt.dist_svmU2(3))then
                      f1=(1-sin(((sqrt(aij2)-ha(3))/hb(3))*api))/2 !0<f1<1
                      ESHORT23=ESHORT23-awei(3,i,j)*f1
                      mk=2
                      ag=-f1
                   elseif(aij2.lt.(80/0.87)**2)then
                      f1=(1+sin(((sqrt(aij2)-hc(3))/hd(3))*api))/2 !0<f1<1
                      ESHORT23=ESHORT23+awei(3,i,j)*f1*fw
                      mk=3
                      ag=f1
                   else
                      ESHORT23=ESHORT23+awei(3,i,j)*fw
                      mk=4
                      ag=1
                   endif
                 npp=npp+1
                 if(npp.le.4000)then
                    write(93,*)npp,sqrt(aij2)*0.87,ag,mk
                 endif
ccc   
                endif
             enddo
 421         continue
          endif
 11    continue
c^^^^^^^^^^^^^^^^^^ contact restrains finished ^^^^^^^^^^^^^^^^^^^^^

c3ccccccccccccccccccccccccccccccccccccc
c     6, 'svmseqCB.dat'
cccccccccccccccccccccccccccccccccccccc
ccc   Contact restrain from 'svmseqCB678.dat' ---------------------->
       do 12 i=iiii,jjjj
          if(Mcon(2,i).ge.1)then
             iseq=seq(i)
             if(mv(i).gt.0)then
                abxi=x(i)+HX(ica(i-1),ica(i),iseq) !CB
                abyi=y(i)+HY(ica(i-1),ica(i),iseq) !CB
                abzi=z(i)+HZ(ica(i-1),ica(i),iseq) !CB
             else
                abxi=etx(i)     !CB
                abyi=ety(i)
                abzi=etz(i)
             endif
             
ccc   'svmseqCB8.dat':
             do k=1,Mcon(2,i)   !number of CB-contacts on i
                j=Kcon(2,i,k)   !j is k'th CB-contact with i
                if(j.lt.i.OR.j.gt.jjjj)then
                   jseq=seq(j)
                   if(mv(j).gt.0)then
                      abxj=x(j)+hx(ica(j-1),ica(j),jseq)
                      abyj=y(j)+hy(ica(j-1),ica(j),jseq)
                      abzj=z(j)+hz(ica(j-1),ica(j),jseq)
                   else
                      abxj=etx(j)
                      abyj=ety(j)
                      abzj=etz(j)
                   endif
                   aij2=(abxi-abxj)**2+(abyi-abyj)**2+(abzi-abzj)**2
ccc   
                   if(aij2.lt.dist_svm2(2))then
                      ESHORT22=ESHORT22-awei(2,i,j)
                   elseif(aij2.lt.dist_svmU2(2))then
                      f1=(1-sin(((sqrt(aij2)-ha(2))/hb(2))*api))/2 !0<f1<1
                      ESHORT22=ESHORT22-awei(2,i,j)*f1
                   elseif(aij2.lt.(80/0.87)**2)then
                      f1=(1+sin(((sqrt(aij2)-hc(2))/hd(2))*api))/2 !0<f1<1
                      ESHORT22=ESHORT22+awei(2,i,j)*f1*fw
                   else
                      ESHORT22=ESHORT22+awei(2,i,j)*fw
                   endif
c                   npp=npp+1
c                   if(npp.le.2000)then
c                      write(93,*)npp,sqrt(aij2)*0.87,ag,mk,npot
c                   endif
ccc   
                endif
             enddo
             
          endif
 12    continue
c^^^^^^^^^^^^^^^^^^CB-contact restrains finished ^^^^^^^^^^^^^^^^^^^^^

c4ccccccccccccccccccccccccccccccccccccc
c     7, 'combCA.dat'
c     9, 'svmseqCA.dat'
cccccccccccccccccccccccccccccccccccccc
ccc   Contact restraint from 'combCA.dat', 'comb8CA.dat', 'svmseqCA.dat'------>
c      d_CA_cut2=(6.0/0.87)**2
      do 33 i=iiii,jjjj
         if(McomCA(i).ge.1.or.Mcon(1,i).ge.1)then
            if(mv(i).gt.0)then
               axi=x(i)
               ayi=y(i)
               azi=z(i)
            else
               axi=ex(i)
               ayi=ey(i)
               azi=ez(i)
            endif
            
ccc   'combCA.dat':
            do 41 k=1,McomCA(i)
               j=KcomCA(i,k)    !k'th contact with i
               if(j.lt.i.OR.j.gt.jjjj)then
                  aij2=(axi-aax(j))**2+(ayi-aay(j))**2+(azi-aaz(j))**2
ccc   
                  if(aij2.lt.d_CA_cut2)then
                     ESHORT9=ESHORT9-aweigCA(i,j)
                  endif
ccc   
               endif
 41         continue
            
ccc   'svmseqCA8.dat': (in fact, there is no svmseqCA8.dat anymore since all through CB)
           do 431 k=1,Mcon(1,i)
              j=Kcon(1,i,k)     !k'th contact with i
              if(j.lt.i.OR.j.gt.jjjj)then
                 aij2=(axi-aax(j))**2+(ayi-aay(j))**2+
     &                (azi-aaz(j))**2
ccc   
                 if(aij2.lt.dist_svm2(1))then
                    ESHORT21=ESHORT21-awei(1,i,j)
                 elseif(aij2.lt.dist_svmU2(1))then
                    f1=(1-sin(((sqrt(aij2)-ha(1))/hb(1))*api))/2 !0<f1<1
                    ESHORT21=ESHORT21-awei(1,i,j)*f1
                 elseif(aij2.lt.(80/0.87)**2)then
                    f1=(1+sin(((sqrt(aij2)-hc(1))/hd(1))*api))/2 !0<f1<1
                    ESHORT21=ESHORT21+awei(1,i,j)*f1*fw
                 else
                    ESHORT21=ESHORT21+awei(1,i,j)*fw
                 endif
c                 npp=npp+1
c                 if(npp.le.4000)then
c                    write(93,*)npp,sqrt(aij2)*0.87,ag,mk
c                 endif
ccc   
              endif
 431       continue
        endif
 33   continue
c^^^^^^^^^^^^^^^^^^ CAcontact restrains finished ^^^^^^^^^^^^^^^^^^^^^
      
*********E13 --------------------------------------->
      i1=max(iiii-1,1)
      i2=min(jjjj-1,Lch-2)
      do i=i1,i2
         ar13=(aax(i)-aax(i+2))**2+(aay(i)-aay(i+2))**2+
     $        (aaz(i)-aaz(i+2))**2
         if(ar13.lt.48) then    !6.03A
            ESHORT6=ESHORT6+csr(i,1)
         else
            ESHORT6=ESHORT6+csr(i,2)
         endif
      enddo
c^^^^^^^^^^^^^^ E_13 finished ^^^^^^^^^^^^^^^^^^^^^^^^^
      
*********E14, E15, proteinlike bias1 to possible alpha/beta --------------->
      i1=max(iiii-4,1)
      i2=min(jjjj,Lch-4)
      do 2 i=i1,i2
         if(mv(i).gt.0)then
            ax1=x(i)
            ay1=y(i)
            az1=z(i)
            agx1=x(i)+gx(ica(i-1),ica(i),seq(i))
            agy1=y(i)+gy(ica(i-1),ica(i),seq(i))
            agz1=z(i)+gz(ica(i-1),ica(i),seq(i))
         else
            ax1=ex(i)
            ay1=ey(i)
            az1=ez(i)
            agx1=egx(i)
            agy1=egy(i)
            agz1=egz(i)
         endif
         if(mv(i+1).gt.0)then
            ax2=x(i+1)
            ay2=y(i+1)
            az2=z(i+1)
            agx2=x(i+1)+gx(ica(i),ica(i+1),seq(i+1))
            agy2=y(i+1)+gy(ica(i),ica(i+1),seq(i+1))
            agz2=z(i+1)+gz(ica(i),ica(i+1),seq(i+1))
         else
            ax2=ex(i+1)
            ay2=ey(i+1)
            az2=ez(i+1)
            agx2=egx(i+1)
            agy2=egy(i+1)
            agz2=egz(i+1)
         endif
         if(mv(i+2).gt.0)then
            ax3=x(i+2)
            ay3=y(i+2)
            az3=z(i+2)
            agx3=x(i+2)+gx(ica(i+1),ica(i+2),seq(i+2))
            agy3=y(i+2)+gy(ica(i+1),ica(i+2),seq(i+2))
            agz3=z(i+2)+gz(ica(i+1),ica(i+2),seq(i+2))
         else
            ax3=ex(i+2)
            ay3=ey(i+2)
            az3=ez(i+2)
            agx3=egx(i+2)
            agy3=egy(i+2)
            agz3=egz(i+2)
         endif
         if(mv(i+3).gt.0)then
            ax4=x(i+3)
            ay4=y(i+3)
            az4=z(i+3)
            agx4=x(i+3)+gx(ica(i+2),ica(i+3),seq(i+3))
            agy4=y(i+3)+gy(ica(i+2),ica(i+3),seq(i+3))
            agz4=z(i+3)+gz(ica(i+2),ica(i+3),seq(i+3))
         else
            ax4=ex(i+3)
            ay4=ey(i+3)
            az4=ez(i+3)
            agx4=egx(i+3)
            agy4=egy(i+3)
            agz4=egz(i+3)
         endif
         if(mv(i+4).gt.0)then
            ax5=x(i+4)
            ay5=y(i+4)
            az5=z(i+4)
            agx5=x(i+4)+gx(ica(i+3),ica(i+4),seq(i+4))
            agy5=y(i+4)+gy(ica(i+3),ica(i+4),seq(i+4))
            agz5=z(i+4)+gz(ica(i+3),ica(i+4),seq(i+4))
         else
            ax5=ex(i+4)
            ay5=ey(i+4)
            az5=ez(i+4)
            agx5=egx(i+4)
            agy5=egy(i+4)
            agz5=egz(i+4)
         endif
ccccccE14:
         ax=ax2-ax1
         ay=ay2-ay1
         az=az2-az1
         bx=ax3-ax2
         by=ay3-ay2
         bz=az3-az2
         cx=ax4-ax3
         cy=ay4-ay3
         cz=az4-az3
         abx=ay*bz-az*by
         aby=az*bx-ax*bz
         abz=ax*by-ay*bx
         hand=abx*cx+aby*cy+abz*cz !a(x)b.c, chirality, >0, right-hand
         ar14=(ax1-ax4)**2+(ay1-ay4)**2+(az1-az4)**2
         r14=nint(ar14)
         if(r14.gt.300) r14=300
         if(hand.lt.0) r14=-r14 !<0, left-hand three-bond.
         ESHORT7=ESHORT7+asr(i,ibin(r14)) !asr(i,dis(i,i+4)) from r14*.dat
cccccccccE15:
         ar15=(ax1-ax5)**2+(ay1-ay5)**2+(az1-az5)**2
         r15=nint(ar15)
         if(r15.gt.500) r15=500
         ESHORT8=ESHORT8+bsr(i,jbin(r15))

cccccccccbias1: encourage helix/sheet-like structure
         if(ar15.lt.75)THEN     !7.53A: helix
            if(sec(i+1).ne.4.AND.sec(i+2).ne.4.AND.sec(i+3).ne.4)then
               if(ibin(r14).gt.4.AND.ibin(r14).lt.8)then !right-hand,3A->8A
                  dot13=(ax4-ax3)*(ax2-ax1)
     $                 +(ay4-ay3)*(ay2-ay1)+(az4-az3)*(az2-az1)
                  if(dot13.lt.0)then
                     dot24=(ax5-ax4)*(ax3-ax2)
     $                    +(ay5-ay4)*(ay3-ay2)+(az5-az4)*(az3-az2)
                     if(dot24.lt.0)then
                        dot14=(ax5-ax4)*(ax2-ax1)
     $                       +(ay5-ay4)*(ay2-ay1)+(az5-az4)*(az2-az1)
                        if(dot14.gt.0)then
                           ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2
                           ESHORT5c=ESHORT5c-2-ff
                        endif
                     endif
                  endif
               endif
            endif
         elseif(ar15.gt.160)then !11A, beta
            if(sec(i+1).ne.2.AND.sec(i+2).ne.2.AND.sec(i+3).ne.2)then
               if(mv(i+1).gt.0)then
                  bx2=hbx(ica(i),ica(i+1))
                  by2=hby(ica(i),ica(i+1))
                  bz2=hbz(ica(i),ica(i+1))
               else
                  bx2=ebx(i+1)  !Hb
                  by2=eby(i+1)
                  bz2=ebz(i+1)
               endif
               if(mv(i+3).gt.0)then
                  bx4=hbx(ica(i+2),ica(i+3))
                  by4=hby(ica(i+2),ica(i+3))
                  bz4=hbz(ica(i+2),ica(i+3))
               else
                  bx4=ebx(i+3)
                  by4=eby(i+3)
                  bz4=ebz(i+3)
               endif
               dot24=bx2*bx4+by2*by4+bz2*bz4 !cos(H1.H3)
               if(dot24.gt.0.71)then !angle<45
                  if(mv(i+2).gt.0)then
                     bx3=hbx(ica(i+1),ica(i+2))
                     by3=hby(ica(i+1),ica(i+2))
                     bz3=hbz(ica(i+1),ica(i+2))
                  else
                     bx3=ebx(i+2)
                     by3=eby(i+2)
                     bz3=ebz(i+2)
                  endif
                  dot23=bx2*bx3+by2*by3+bz2*bz3 !cos(H1.H2)
                  if(dot23.lt.-0.71)then !angle >135
                     ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2
c                     ESHORT5c=ESHORT5c-2-ff
                     ESHORT5c=ESHORT5c-2-ff*2
                  endif
               endif
            endif
         endif
 2    continue
c^^^^^^^^^^E14, E_15, bias1 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^

****  Protein-like bias2: b(i), b(i+4) parallel------->
      i1=max(iiii-4,1)
      i2=min(jjjj-1,Lch-5)
      do i=i1,i2
         ff=(afs(i)*afs(i+1)+afs(i+3)*afs(i+4))/2 !ff=1, r<r0; =0.25, r>2r0.
         if(mv(i).gt.0)then
            bx1=hbx(ica(i-1),ica(i))
            by1=hby(ica(i-1),ica(i))
            bz1=hbz(ica(i-1),ica(i))
         else
            bx1=ebx(i)
            by1=eby(i)
            bz1=ebz(i)
         endif
         if(mv(i+4).gt.0)then
            bx5=hbx(ica(i+3),ica(i+4))
            by5=hby(ica(i+3),ica(i+4))
            bz5=hbz(ica(i+3),ica(i+4))
         else
            bx5=ebx(i+4)
            by5=eby(i+4)
            bz5=ebz(i+4)
         endif
         b15=bx1*bx5+by1*by5+bz1*bz5
         if(seq(i).eq.2.and.seq(i+2).eq.2.and.seq(i+4).eq.2)then
            if(b15.gt.0.9)then
               ESHORT5=ESHORT5-ff !alpha
            endif
         else
            if(b15.gt.0.5.or.b15.lt.-0.3)then !beta or turn
c               ESHORT5=ESHORT5-ff
               ESHORT5=ESHORT5-ff*2
            endif
         endif
      enddo
c^^^^^^^^^^^^^^^^^^ bias2 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

****  Protein-like bias3: c(i),c(i+2) anti/parallel:
      i1=max(iiii-2,1)
      i2=min(jjjj,Lch-2)
      do i=i1,i2
         ff=(afs(i)+afs(i+1)+afs(i+2))/3
         if(mv(i).gt.0)then
            cx1=cax(ica(i-1),ica(i))
            cy1=cay(ica(i-1),ica(i))
            cz1=caz(ica(i-1),ica(i))
         else
            cx1=ecx(i)
            cy1=ecy(i)
            cz1=ecz(i)
         endif
         if(mv(i+2).gt.0)then
            cx3=cax(ica(i+1),ica(i+2))
            cy3=cay(ica(i+1),ica(i+2))
            cz3=caz(ica(i+1),ica(i+2))
         else
            cx3=ecx(i+2)
            cy3=ecy(i+2)
            cz3=ecz(i+2)
         endif
         c13=abs(cx1*cx3+cy1*cy3+cz1*cz3)
         c13=min(0.71,c13)/0.71 !c13 is the same in [0,45]
         ESHORT5=ESHORT5-ff*c13
      enddo
c^^^^^^^^^^^^^^^^^^ bias3 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ccccccbias4a: to predicted alpha fragment --------------->
      i1=max(iiii-6,1)
      i2=min(jjjj-1,Lch-7)
      do i=i1,i2
         if(frga(i).gt.1)then
            dis=sqrt((aax(i)-aax(i+7))**2+
     $           (aay(i)-aay(i+7))**2+(aaz(i)-aaz(i+7))**2)
            ESHORT5b=ESHORT5b+abs(dis-frga(i))
         endif
      enddo
ccccccbias4b: to predicted beta fragment --------------->
      i1=max(iiii-5,1)
      i2=min(jjjj-1,Lch-6)
      do i=i1,i2
         if(frgb(i).gt.1)then
            dis=sqrt((aax(i)-aax(i+6))**2+
     $           (aay(i)-aay(i+6))**2+(aaz(i)-aaz(i+6))**2)
c            ESHORT5b=ESHORT5b+abs(dis-frgb(i))
            ESHORT5b=ESHORT5b+abs(dis-frgb(i))*2
         endif
      enddo
c^^^^^^^^^^^^^^^ predicted alpha/beta bias finished ^^^^^^^^^^^^^^^

ccc   penality for crumpling structures ---------------------------->
      i1=max(iiii-11,1)
      i2=min(jjjj-1,Lch-12)
      do i=i1,i2
         if(mv(i).gt.0)then
            ax0=x(i)
            ay0=y(i)
            az0=z(i)
         else
            ax0=ex(i)
            ay0=ey(i)
            az0=ez(i)
         endif
         if(mv(i+4).gt.0)then
            ax4=x(i+4)
            ay4=y(i+4)
            az4=z(i+4)
         else
            ax4=ex(i+4)
            ay4=ey(i+4)
            az4=ez(i+4)
         endif
         if(mv(i+8).gt.0)then
            ax8=x(i+8)
            ay8=y(i+8)
            az8=z(i+8)
         else
            ax8=ex(i+8)
            ay8=ey(i+8)
            az8=ez(i+8)
         endif
         avx1=ax4-ax0
         avy1=ay4-ay0
         avz1=az4-az0
         avx2=ax8-ax4
         avy2=ay8-ay4
         avz2=az8-az4
         aaa=avx1*avx2+avy1*avy2+avz1*avz2
         if(aaa.lt.0)then
            if(mv(i+12).gt.0)then
               ax12=x(i+12)
               ay12=y(i+12)
               az12=z(i+12)
            else
               ax12=ex(i+12)
               ay12=ey(i+12)
               az12=ez(i+12)
            endif
            avx3=ax12-ax8
            avy3=ay12-ay8
            avz3=az12-az8
            bbb=avx1*avx3+avy1*avy3+avz1*avz3
            if(bbb.gt.0)then
               ccc=avx2*avx3+avy2*avy3+avz2*avz3
               if(ccc.lt.0)then
                  ESHORT5a=ESHORT5a+1 !crumpling
               endif
            endif
         endif
      enddo
c^^^^^^^^^^^^^ penality of bizard structure finished (ESC1) ^^^^^^^^^^^^

c     Further penalize deriviation of restrain, if larger than (colim, dilim):
      IF(istat.eq.10) then      !calculate E from the beginning
cc         if(dinew.gt.dilim) ESHORT=ESHORT+er1*(dinew-dilim)
cc         if(conew.gt.colim) ESHORT=ESHORT+er3*(conew-colim)
         if(conew.gt.colim) ESHORT4a=conew-colim !colim=number of restrains.
      endif
c     codevsum and didevsum are backuped as total deviation for old 
c     conformation after the first istat=10 are finished

c     the panalty was not counted in Enew, the following is the difference
c     of panlity on new and old conformation, i.e. dE=Enew-Eold. Since
c     this part energy appear only at Eold, the sign is oppsite:
      if(istat.lt.0) then       !return to old energy
         codev=codevsum+conew-coold !total deviation for new conformation
         didev=didevsum+dinew-diold !total deviation for new conformation
c     total penalty-energy beacuse of restrain-deviation on old conformation
cc         if(codevsum.gt.colim) ESHORT=ESHORT+er3*(codevsum-colim)
cc         if(didevsum.gt.dilim) ESHORT=ESHORT+er1*(didevsum-dilim)
         if(codevsum.gt.colim) ESHORT4a=codevsum-colim
c     total penalty-energy beacuse of restrain-deviation on new conformation:
cc         if(codev.gt.colim) ESHORT=ESHORT-er3*(codev-colim)
cc         if(didev.gt.dilim) ESHORT=ESHORT-er1*(didev-dilim)
         if(codev.gt.colim) ESHORT4a=ESHORT4a-(codev-colim)
      endif

***   panelty for the deviation from template---------------->
      if(switch.eq.5)then
         if(ISTAT.eq.10)then    !total energy
            do i=1,nfr
               nn=0
               do j=nfr_i(i),nfr_f(i)
                  nn=nn+1
                  r_1(1,nn)=ex0(j)
                  r_1(2,nn)=ey0(j)
                  r_1(3,nn)=ez0(j)
                  r_2(1,nn)=ex(j)
                  r_2(2,nn)=ey(j)
                  r_2(3,nn)=ez(j)
               enddo
               call u3b(w,r_1,r_2,nn,0,rms,u,t,ier)
               armsd=dsqrt(rms/nn) !RMSD is real, rms is double precision
               ESHORT11=ESHORT11+armsd**4
            enddo
         else                   !fragment energy
            if(i_E_defo.eq.1)then !for defo_M but not trot_M
               if(iiii.eq.1)then
                  i=i_chunk(1)
               elseif(iiii.eq.Lch)then
                  i=i_chunk(Lch)
               else
                  i=i_chunk(iiii+2) !iiii+2 is first residue of the fragment
               endif
               if(i.gt.0)then
                  nn=0
                  do j=nfr_i(i),nfr_f(i)
                     nn=nn+1
                     r_1(1,nn)=ex0(j)
                     r_1(2,nn)=ey0(j)
                     r_1(3,nn)=ez0(j)
                     r_2(1,nn)=ex(j)
                     r_2(2,nn)=ey(j)
                     r_2(3,nn)=ez(j)
                  enddo
                  call u3b(w,r_1,r_2,nn,0,rms,u,t,ier)
                  armsd=dsqrt(rms/nn) !RMSD is real, rms is double precision
                  ESHORT11=ESHORT11+armsd**4
               endif
            endif
         endif
      endif
***^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

********distance restraints------------------------->
      if(nstick.eq.1)then
         do i=iiii,jjjj
            if(iq(i,n_tem(itemp)).eq.1)then
               adis=(aax(i)-ax00(i,n_tem(itemp)))**2+
     &              (aay(i)-ay00(i,n_tem(itemp)))**2+
     &              (aaz(i)-az00(i,n_tem(itemp)))**2
               adis=sqrt(adis)
               ESHORT12=ESHORT12+adis 
            endif
         enddo
      endif
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

********distance restraints for RMSD ------------------------->
      if(nrmsd.eq.1)then
         nn=0
         do i=1,Lch
            if(iq(i,n_tem(itemp)).eq.1)then
               nn=nn+1
               r_1(1,nn)=ax00(i,n_tem(itemp))
               r_1(2,nn)=ay00(i,n_tem(itemp))
               r_1(3,nn)=az00(i,n_tem(itemp))
               r_2(1,nn)=aax(i)
               r_2(2,nn)=aay(i)
               r_2(3,nn)=aaz(i)
            endif
         enddo
         call u3b(w,r_1,r_2,nn,0,rms,u,t,ier)
         armsd=dsqrt(rms/nn)    !RMSD is real, rms is double precision
         ESHORT13=armsd**2
         if(armsd.lt.armsd_min0)then
           armsd_min0=armsd
           itemp0=itemp
           icycle0=icycle
           do i=1,Lch
             bx00(i)=aax(i)
             by00(i)=aay(i)
             bz00(i)=aaz(i)
c             write(*,*)i,bx00(i),by00(i),bz00(i)
           enddo
         endif
      else
         ESHORT13=0
      endif
*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      ESHORT=
     $     +es2*ESHORT2
     $     +er1*ESHORT3
     $     +er3*ESHORT4
     $     +er4*ESHORT4a
     $     +es3*ESHORT5
     $     +es3a*ESHORT5a
     $     +es3b*ESHORT5b
     $     +es3c*ESHORT5c
     $     +es4*ESHORT6
     $     +es5*ESHORT7
     $     +es6*ESHORT8
     $     +er5*ESHORT9
     $     +er6*ESHORT10
     $     +er7*ESHORT11
     $     +astick*ESHORT12
     $     +ermsd*ESHORT13
     $     +er21*ESHORT21       !er21=1, weight tuned by contact.map, CA
     $     +er22*ESHORT22       !er22=1, weight tuned by contact.map, CB
     $     +er23*ESHORT23       !er23=1, weight tuned by contact.map, SG
      
c      write(*,*)iiii,jjjj,ESHORT,ESHORT2,ESHORT3,
c     &     ESHORT4,ESHORT4a,ESHORT5,ESHORT5a,ESHORT5b,
c     &     ESHORT5c,ESHORT6,ESHORT7,ESHORT8,ESHORT9,
c     &     ESHORT10,ESHORT11,ESHORT12,ESHORT13,
c     &     ESHORT21,ESHORT22,ESHORT23
c     write(*,*)eshort,eshort12,astick,itemp
c      write(*,*)iiii,jjjj,ESHORT21,ESHORT22,ESHORT23

c ^^^^^^^^^^ E_short finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  2-bond movement: prefabricated
c  residues of 'm, m1=m+1,m2=m+2' are involved,
c  but positions of m1 are changed.
c  old ---> new ---> old ----> new
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move2
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      double precision bNa,bNt,bNNa,bNNt
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy

c     choose the position of (m) ----------->
      i=int(aranzy(nozy)*nfl2)+1  ![1,nfl2]
      m=ras2(i)
      m1=m+1
      m2=m+2

cccc  all the pairs have at least one another pair, so nc>=2
      nc=Np2(ica(m),ica(m1))    !number of pathes from m to m2

c     choose p'th new path from (m) to (m2) -------------->
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m))) goto 201
      nn(m1)=v22(ica(m),ica(m1),p)
      if(.not.goodc(nn(m1),ica(m2))) goto 201
c^^^^^^^^^^ new conformation chosen finished ^^^^^^^^^^^^^^^^^^^^^^^^
      if(nn(m).eq.ica(m))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m2))then        ! check excluded volumn for passage of [m,m2]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m2,1)+ESHORT(m,m2,1) !icnt,nop are repeated.
         do kkk=m,m2
            afsn(kkk)=afs(kkk)  !afs(i) is useful in ESHORT
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m2,-1)+ESHORT(m,m2,-1) !repeated part of icnt, nop removed

c     calculate eprofn while dord was calculated when call EHB(m,m2,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)          
            ia=noa(pp)          
            ip=nop(pp)          !nopp(old)+1(new)-1(old)=new
            im=nom(pp)
c            if(ip.lt.0.or.ip.gt.15.or
c     &           .ia.lt.0.or.ia.gt.15.or
c     &           .im.lt.0.or.im.gt.15)then
c               write(*,*)'ia,ip,im=',ia,ip,im
c            endif
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(2)=bNa(2)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt          !backup1
            sumcto=sumct        !backup2
            do pp=1,Lch
               nopp(pp)=nop(pp) !backup3
               nomm(pp)=nom(pp) !backup4
               noaa(pp)=noa(pp) !backup5
            enddo
            eprofo=eprofn       !backup6
            codevsum=codev      !backup7
            didevsum=didev      !backup8
            do kkk=m,m2
               afs(kkk)=afsn(kkk) !backup9
            enddo

c            energ=energ+de
c     change the conformation to the new position--------->
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif                  !for id
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
 113  continue
      bNt(2)=bNt(2)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move2 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  3-bond movement: single-move
c  residues of 'i,i+1,i+2,i+3' are involved,
c  but positions bond-vectors of i+1, i+2 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move3s
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy

c     choose the position of (m) ----------->
      i=int(aranzy(nozy)*nfl3)+1
      m=ras3(i)
      m1=m+1
      m2=m+2
      m3=m+3

ccccccccccccc 1th 2-bond movement cccccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m)))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccc 2th 2-bond movement cccccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m1)=v21(t1,ica(m2),p)
      if(.not.goodc(nn(m),nn(m1)))goto 202
      nn(m2)=v22(t1,ica(m2),p)
      if(.not.goodc(nn(m2),ica(m3)))goto 202
c^^^^^^^^ three-bonds movement conformation finished ^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m3))then        ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m3,1)+ESHORT(m,m3,1) !use ifs
         do kkk=m,m3
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m3,-1)+ESHORT(m,m3,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(3)=bNa(3)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m3
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
 113  continue
      bNt(3)=bNt(3)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move3s finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  3-bond movement: double move
c  residues of 'i,i+1,i+2,i+3' are involved,
c  but positions bond-vectors of i+1, i+2 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move3d
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy

c     choose the position of (m) ----------->
      i=int(aranzy(nozy)*nfl3)+1
      m=ras3(i)
      m1=m+1
      m2=m+2
      m3=m+3

ccccccccccccc temporal 1th 2-bond movement cccccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccc temporal 2th 2-bond movement cccccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202

ccccccccccccc 1th 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(t,tt1)
 203  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 203
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,t2))goto 203
ccccccccccccc 2th 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt1,t2)
 204  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,t2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 204
      nn(m2)=v22(ttt1,t2,p)
      if(.not.goodc(nn(m2),ica(m3)))goto 204
c^^^^^^^^ three-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1))goto 113


c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m3))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m3,1)+ESHORT(m,m3,1) !use ifs
         do kkk=m,m3
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m3,-1)+ESHORT(m,m3,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(4)=bNa(4)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m3
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
113   continue
      bNt(4)=bNt(4)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move3d finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     4-bond movement: single move
c     residues of 'i,i+1,i+2,i+3,i+4' are involved,
c     but positions bond-vectors of i+1, i+2, i+3 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move4s
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy

c     choose the position of (m) ----------->
      i=int(aranzy(nozy)*nfl4)+1
      m=ras4(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4

ccccccccccccccc 1th 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m)))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m1)=v21(t1,ica(m2),p)
      if(.not.goodc(nn(m),nn(m1)))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m2)=v21(t2,ica(m3),p)
      if(.not.goodc(nn(m1),nn(m2)))goto 203
      nn(m3)=v22(t2,ica(m3),p)
      if(.not.goodc(nn(m3),ica(m4)))goto 203
c^^^^^^^^ four-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1).
     $     and.nn(m2).eq.ica(m2))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      ox(m3)=x(m3)              !memory of old path
      oy(m3)=y(m3)              !memory of old path
      oz(m3)=z(m3)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      oo(m3)=ica(m3)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path
      nx(m3)=nx(m2)+vx(nn(m2))  !memory of new path
      ny(m3)=ny(m2)+vy(nn(m2))  !memory of new path
      nz(m3)=nz(m2)+vz(nn(m2))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      x(m3)=nx(m3)
      y(m3)=ny(m3)
      z(m3)=nz(m3)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(m3)=nn(m3)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m4))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m4,1)+ESHORT(m,m4,1) !use ifs
         do kkk=m,m4
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m4,-1)+ESHORT(m,m4,-1)
         
c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(5)=bNa(5)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m4
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            x(m3)=nx(m3)
            y(m3)=ny(m3)
            z(m3)=nz(m3)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(m3)=nn(m3)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
113   continue
      bNt(5)=bNt(5)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move4s finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     4-bond movement: double move
c     residues of 'i,i+1,i+2,i+3,i+4' are involved,
c     but positions bond-vectors of i+1, i+2, i+3 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move4d
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy

c     choose the position of (m) ----------->
      i=int(aranzy(nozy)*nfl4)+1
      m=ras4(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4

ccccccccccccccc 1th temperor 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th temperor 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th temperor 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      tt2=v21(t2,ica(m3),p)
      if(.not.goodc(tt1,tt2))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203

ccccccccccccc first 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(t,tt1)
 204  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 204
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,tt2))goto 204
ccccccccccccc second 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt1,tt2)
 205  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,tt2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 205
      ttt2=v22(ttt1,tt2,p)
      if(.not.goodc(ttt2,t3))goto 205
ccccccccccccc third 2-bond movement cccccccccccccccccccccccccc
      nc=Np2(ttt2,t3)
 206  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m2)=v21(ttt2,t3,p)
      if(.not.goodc(nn(m1),nn(m2)))goto 206
      nn(m3)=v22(ttt2,t3,p)
      if(.not.goodc(nn(m3),ica(m4)))goto 206
c^^^^^^^^ four-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1).
     $     and.nn(m2).eq.ica(m2))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      ox(m3)=x(m3)              !memory of old path
      oy(m3)=y(m3)              !memory of old path
      oz(m3)=z(m3)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      oo(m3)=ica(m3)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path
      nx(m3)=nx(m2)+vx(nn(m2))  !memory of new path
      ny(m3)=ny(m2)+vy(nn(m2))  !memory of new path
      nz(m3)=nz(m2)+vz(nn(m2))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      x(m3)=nx(m3)
      y(m3)=ny(m3)
      z(m3)=nz(m3)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(m3)=nn(m3)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m4))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m4,1)+ESHORT(m,m4,1) !use ifs
         do kkk=m,m4
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m4,-1)+ESHORT(m,m4,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(6)=bNa(6)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m4
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            x(m3)=nx(m3)
            y(m3)=ny(m3)
            z(m3)=nz(m3)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(m3)=nn(m3)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
113   continue
      bNt(6)=bNt(6)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move4d finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     5-bond movement: single move
c     residues of 'i,i+1,i+2,i+3,i+4' are involved,
c     but positions bond-vectors of i+1, i+2, i+3 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move5s
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy

c     choose the position of (m) ----------->
      i=int(aranzy(nozy)*nfl5)+1
      m=ras5(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4
      m5=m+5

ccccccccccccccc 1th 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),nn(m)))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc 2th 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m1)=v21(t1,ica(m2),p)
      if(.not.goodc(nn(m),nn(m1)))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc 3th 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m2)=v21(t2,ica(m3),p)
      if(.not.goodc(nn(m1),nn(m2)))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203
ccccccccccccccc 4th 2-bond movement cccccccccccccccc
      nc=Np2(t3,ica(m4))
 204  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m3)=v21(t3,ica(m4),p)
      if(.not.goodc(nn(m2),nn(m3)))goto 204
      nn(m4)=v22(t3,ica(m4),p)
      if(.not.goodc(nn(m4),ica(m5)))goto 204
c^^^^^^^^ four-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1).and.
     $     nn(m2).eq.ica(m2).and.nn(m3).eq.ica(m3))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      ox(m3)=x(m3)              !memory of old path
      oy(m3)=y(m3)              !memory of old path
      oz(m3)=z(m3)              !memory of old path
      ox(m4)=x(m4)              !memory of old path
      oy(m4)=y(m4)              !memory of old path
      oz(m4)=z(m4)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      oo(m3)=ica(m3)            !memory of old path
      oo(m4)=ica(m4)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path
      nx(m3)=nx(m2)+vx(nn(m2))  !memory of new path
      ny(m3)=ny(m2)+vy(nn(m2))  !memory of new path
      nz(m3)=nz(m2)+vz(nn(m2))  !memory of new path
      nx(m4)=nx(m3)+vx(nn(m3))  !memory of new path
      ny(m4)=ny(m3)+vy(nn(m3))  !memory of new path
      nz(m4)=nz(m3)+vz(nn(m3))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      x(m3)=nx(m3)
      y(m3)=ny(m3)
      z(m3)=nz(m3)
      x(m4)=nx(m4)
      y(m4)=ny(m4)
      z(m4)=nz(m4)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(m3)=nn(m3)
      ica(m4)=nn(m4)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m5))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,m5,1)+ESHORT(m,m5,1) !use ifs
         do kkk=m,m5
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         x(m4)=ox(m4)
         y(m4)=oy(m4)
         z(m4)=oz(m4)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(m4)=oo(m4)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m5,-1)+ESHORT(m,m5,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(8)=bNa(8)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m5
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            x(m3)=nx(m3)
            y(m3)=ny(m3)
            z(m3)=nz(m3)
            x(m4)=nx(m4)
            y(m4)=ny(m4)
            z(m4)=nz(m4)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(m3)=nn(m3)
            ica(m4)=nn(m4)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         x(m4)=ox(m4)
         y(m4)=oy(m4)
         z(m4)=oz(m4)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(m4)=oo(m4)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
113   continue
      bNt(8)=bNt(8)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move5s finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     5-bond movement: double move
c     residues of 'i,i+1,i+2,i+3,i+4,i+5' are involved,
c     but positions bond-vectors of i+1, i+2, i+3, i+4 are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move5d
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy

c     choose the position of (m) ----------->
      i=int(aranzy(nozy)*nfl5)+1
      m=ras5(i)
      m1=m+1
      m2=m+2
      m3=m+3
      m4=m+4
      m5=m+5

ccccccccccccccc temperor 1th 2-bond movement cccccccccccccccc
      nc=Np2(ica(m),ica(m1))
 201  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t=v21(ica(m),ica(m1),p)
      if(.not.goodc(ica(m-1),t))goto 201
      t1=v22(ica(m),ica(m1),p)
      if(.not.goodc(t1,ica(m2)))goto 201
ccccccccccccccc temperor 2th 2-bond movement cccccccccccccccc
      nc=Np2(t1,ica(m2))
 202  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      tt1=v21(t1,ica(m2),p)
      if(.not.goodc(t,tt1))goto 202
      t2=v22(t1,ica(m2),p)
      if(.not.goodc(t2,ica(m3)))goto 202
ccccccccccccccc temperor 3th 2-bond movement cccccccccccccccc
      nc=Np2(t2,ica(m3))
 203  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      tt2=v21(t2,ica(m3),p)
      if(.not.goodc(tt1,tt2))goto 203
      t3=v22(t2,ica(m3),p)
      if(.not.goodc(t3,ica(m4)))goto 203
ccccccccccccccc temperor 4th 2-bond movement cccccccccccccccc
      nc=Np2(t3,ica(m4))
 204  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      tt3=v21(t3,ica(m4),p)
      if(.not.goodc(tt2,tt3))goto 204
      t4=v22(t3,ica(m4),p)
      if(.not.goodc(t4,ica(m5)))goto 204

ccccccccccccccc 1th 2-bond movement cccccccccccccccc
      nc=Np2(t,tt1)
 205  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m)=v21(t,tt1,p)
      if(.not.goodc(ica(m-1),nn(m)))goto 205
      ttt1=v22(t,tt1,p)
      if(.not.goodc(ttt1,tt2))goto 205
ccccccccccccccc 2th 2-bond movement cccccccccccccccc
      nc=Np2(ttt1,tt2)
 206  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m1)=v21(ttt1,tt2,p)
      if(.not.goodc(nn(m),nn(m1)))goto 206
      ttt2=v22(ttt1,tt2,p)
      if(.not.goodc(ttt2,tt3))goto 206
ccccccccccccccc 3th 2-bond movement cccccccccccccccc
      nc=Np2(ttt2,tt3)
 207  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m2)=v21(ttt2,tt3,p)
      if(.not.goodc(nn(m1),nn(m2)))goto 207
      ttt3=v22(ttt2,tt3,p)
      if(.not.goodc(ttt3,t4))goto 207
ccccccccccccccc 4th 2-bond movement cccccccccccccccc
      nc=Np2(ttt3,t4)
 208  p=int(aranzy(nozy)*nc)+1    ![1,nc]
      nn(m3)=v21(ttt3,t4,p)
      if(.not.goodc(nn(m2),nn(m3)))goto 208
      nn(m4)=v22(ttt3,t4,p)
      if(.not.goodc(nn(m4),ica(m5)))goto 208
c^^^^^^^^ four-bonds movement conformation finished ^^^^^^^^^^
      if(nn(m).eq.ica(m).and.nn(m1).eq.ica(m1).and.
     $     nn(m2).eq.ica(m2).and.nn(m3).eq.ica(m3))goto 113

c     prepare new path and backup old path ------------>
      ox(m1)=x(m1)              !memory of old path
      oy(m1)=y(m1)              !memory of old path
      oz(m1)=z(m1)              !memory of old path
      ox(m2)=x(m2)              !memory of old path
      oy(m2)=y(m2)              !memory of old path
      oz(m2)=z(m2)              !memory of old path
      ox(m3)=x(m3)              !memory of old path
      oy(m3)=y(m3)              !memory of old path
      oz(m3)=z(m3)              !memory of old path
      ox(m4)=x(m4)              !memory of old path
      oy(m4)=y(m4)              !memory of old path
      oz(m4)=z(m4)              !memory of old path
      oo(m)=ica(m)              !memory of old path
      oo(m1)=ica(m1)            !memory of old path
      oo(m2)=ica(m2)            !memory of old path
      oo(m3)=ica(m3)            !memory of old path
      oo(m4)=ica(m4)            !memory of old path
      nx(m1)=x(m)+vx(nn(m))     !memory of new path
      ny(m1)=y(m)+vy(nn(m))     !memory of new path
      nz(m1)=z(m)+vz(nn(m))     !memory of new path
      nx(m2)=nx(m1)+vx(nn(m1))  !memory of new path
      ny(m2)=ny(m1)+vy(nn(m1))  !memory of new path
      nz(m2)=nz(m1)+vz(nn(m1))  !memory of new path
      nx(m3)=nx(m2)+vx(nn(m2))  !memory of new path
      ny(m3)=ny(m2)+vy(nn(m2))  !memory of new path
      nz(m3)=nz(m2)+vz(nn(m2))  !memory of new path
      nx(m4)=nx(m3)+vx(nn(m3))  !memory of new path
      ny(m4)=ny(m3)+vy(nn(m3))  !memory of new path
      nz(m4)=nz(m3)+vz(nn(m3))  !memory of new path

c     change conformation to new path -------------->
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      x(m2)=nx(m2)
      y(m2)=ny(m2)
      z(m2)=nz(m2)
      x(m3)=nx(m3)
      y(m3)=ny(m3)
      z(m3)=nz(m3)
      x(m4)=nx(m4)
      y(m4)=ny(m4)
      z(m4)=nz(m4)
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(m2)=nn(m2)
      ica(m3)=nn(m3)
      ica(m4)=nn(m4)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,m5))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         
         Enew=EHB(m,m5,1)+ESHORT(m,m5,1) !use ifs
         do kkk=m,m5
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         x(m4)=ox(m4)
         y(m4)=oy(m4)
         z(m4)=oz(m4)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(m4)=oo(m4)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,m5,-1)+ESHORT(m,m5,-1)
         
c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(9)=bNa(9)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,m5
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            x(m2)=nx(m2)
            y(m2)=ny(m2)
            z(m2)=nz(m2)
            x(m3)=nx(m3)
            y(m3)=ny(m3)
            z(m3)=nz(m3)
            x(m4)=nx(m4)
            y(m4)=ny(m4)
            z(m4)=nz(m4)
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(m2)=nn(m2)
            ica(m3)=nn(m3)
            ica(m4)=nn(m4)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         x(m2)=ox(m2)
         y(m2)=oy(m2)
         z(m2)=oz(m2)
         x(m3)=ox(m3)
         y(m3)=oy(m3)
         z(m3)=oz(m3)
         x(m4)=ox(m4)
         y(m4)=oy(m4)
         z(m4)=oz(m4)
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(m2)=oo(m2)
         ica(m3)=oo(m3)
         ica(m4)=oo(m4)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
113   continue
      bNt(9)=bNt(9)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move5d finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     6-bond movement: translation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move6
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

c      nob=6+int(aranzy(nozy)*6.999) !number of moved bonds, [6,12]
      nob=6

c     choose the position (m) to be moved
      i=int(aranzy(nozy)*nfl6)+1 ![1,nfl6]
      m=ras6(i)
      m1=m+1
      mnob=m+nob
      mnob1=mnob-1
      mnob2=mnob-2

 201  mx=int(aranzy(nozy)*2.999999)-1 ![-1,0,1]
      my=int(aranzy(nozy)*2.999999)-1 ![-1,0,1]
      mz=int(aranzy(nozy)*2.999999)-1 ![-1,0,1]
      if((mx*mx+my*my+mz*mz).eq.0) goto 201

      wx=vx(ica(m))+mx
      wy=vy(ica(m))+my
      wz=vz(ica(m))+mz
      ir=wx*wx+wy*wy+wz*wz
      if(ir.lt.latt1)goto 201
      if(ir.gt.latt2)goto 201
      nn(m)=vector(wx,wy,wz)    !new vector of vertix 'm'
      if(.not.goodc(ica(m-1),nn(m))) goto 202
      if(.not.goodc(nn(m),ica(m1))) goto 202

      wx=vx(ica(mnob1))-mx
      wy=vy(ica(mnob1))-my
      wz=vz(ica(mnob1))-mz
      ir=wx*wx+wy*wy+wz*wz
      if(ir.lt.latt1) goto 202
      if(ir.gt.latt2) goto 202
      nn(mnob1)=vector(wx,wy,wz) !new vector of vertix 'm+5'
      if(.not.goodc(ica(mnob2),nn(mnob1))) goto 202
      if(.not.goodc(nn(mnob1),ica(mnob))) goto 202

c     prepare new path and backup old path ------------>
      do i=m1,mnob1
         ox(i)=x(i)             !memory of old path
         oy(i)=y(i)             !memory of old path
         oz(i)=z(i)             !memory of old path
         nx(i)=x(i)+mx          !memory of new path
         ny(i)=y(i)+my          !memory of new path
         nz(i)=z(i)+mz          !memory of new path
      enddo
      oo(m)=ica(m)              !memory of old path
      oo(mnob1)=ica(mnob1)      !memory of old path

c     change conformation to new path -------------->
      do i=m1,mnob1
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
      enddo
      ica(m)=nn(m)
      ica(mnob1)=nn(mnob1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      if(look(m,mnob))then      ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,mnob,1)+ESHORT(m,mnob,1) !use ifs
         do kkk=m,mnob
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         do i=m1,mnob1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
         enddo
         ica(m)=oo(m)
         ica(mnob1)=oo(mnob1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m,mnob,-1)+ESHORT(m,mnob,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(10)=bNa(10)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,mnob
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            do i=m1,mnob1
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
            enddo
            ica(m)=nn(m)
            ica(mnob1)=nn(mnob1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         do i=m1,mnob1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
         enddo
         ica(m)=oo(m)
         ica(mnob1)=oo(mnob1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)

 202  continue
      bNt(10)=bNt(10)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move6 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0)
c     Move the off-lattice N-fragment.
c     This movement include an off-lattice rotation of fragment
c     and a two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tran_N
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      n=nfr_f(ifr)+2            !ending point of movement
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
***   check n-terminal of fragment--------------->
      ax1=ax0+ex(n2)
      ay1=ay0+ey(n2)
      az1=az0+ez(n2)
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=1,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=1,n2
         ex_n(i)=ax0+ex(i)      !CA
         ey_n(i)=ay0+ey(i)
         ez_n(i)=az0+ez(i)
         egx_n(i)=ax0+egx(i)    !SG
         egy_n(i)=ay0+egy(i)
         egz_n(i)=az0+egz(i)
         ecx_n(i)=ecx(i)        !cc
         ecy_n(i)=ecy(i)
         ecz_n(i)=ecz(i)
         ebx_n(i)=ebx(i)        !Hb
         eby_n(i)=eby(i)
         ebz_n(i)=ebz(i)
         etx_n(i)=ax0+etx(i)    !CB
         ety_n(i)=ay0+ety(i)
         etz_n(i)=az0+etz(i)
      enddo
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202
      
c     change conformation to new path -------------->
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=1,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo

      if(look(1,n))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(1,n,1)+ESHORT(1,n,1) !use ifs
         do kkk=1,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(1,n,-1)+ESHORT(1,n,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(16)=bNa(16)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=1,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=1,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(16)=bNt(16)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ tran_N finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0) + rotation (angle0)
c     Move the off-lattice N-fragment.
c     This movement include an off-lattice rotation of fragment
c     and a two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trot_N
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      n=nfr_f(ifr)+2            !ending point of movement
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
c     rotation axis (awx,awy,awz):
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
c     rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy(nozy)*2-1)*angle00(ifr) ![-angle00,angle00]
c     rotation points:
      ax=ex(n2)
      ay=ey(n2)
      az=ez(n2)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
***   check n-terminal of fragment--------------->
      ax1=ax0+ax+(ex(n2)-ax)*a11+(ey(n2)-ay)*a12+(ez(n2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(n2)-ax)*a21+(ey(n2)-ay)*a22+(ez(n2)-az)*a23
      az1=az0+az+(ex(n2)-ax)*a31+(ey(n2)-ay)*a32+(ez(n2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=1,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=1,n2
        ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
        ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
        ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
        egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
        egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
        egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
        ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
        ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
        ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
        ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
        eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
        ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
        etx_n(i)=ax0+ax+(etx(i)-ax)*a11+(ety(i)-ay)*a12+(etz(i)-az)*a13 !CB
        ety_n(i)=ay0+ay+(etx(i)-ax)*a21+(ety(i)-ay)*a22+(etz(i)-az)*a23
        etz_n(i)=az0+az+(etx(i)-ax)*a31+(ety(i)-ay)*a32+(etz(i)-az)*a33
      enddo
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=1,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo

      if(look(1,n))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(1,n,1)+ESHORT(1,n,1) !use ifs
         do kkk=1,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=ebx_o(i)     !CB
            ety(i)=eby_o(i)
            etz(i)=ebz_o(i)
         enddo
         Eold=EHB(1,n,-1)+ESHORT(1,n,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(22)=bNa(22)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=1,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=1,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(22)=bNt(22)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ trot_N finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0) + deformation (angle0)
c     Move the off-lattice N-fragment.
c     This movement include an off-lattice rotation of fragment
c     and a two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defo_N
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)
      common/defoangle/defo_angle
      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      n=nfr_f(ifr)+2            !ending point of movement
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
ccc   d_xyz0=0.5               !range of translation movement
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
ccc   rotation axis (awx,awy,awz)->N-terminal; (bwx,bwy,bwz)->C-terminal:
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      bwx=asin_theta*cos(aphi)
      bwy=asin_theta*sin(aphi)
      bwz=acos_theta
ccc   rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy(nozy)*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
      bng=(aranzy(nozy)*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
ccc   rotation points:
      n_rot=1+int(n2*aranzy(nozy)) ![1,n2]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
***   rotation matrix:
      acos=cos(ang)             !for N-terminal of the fragments--->
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
      acos=cos(ang)             !for C-terminal of the fragments--->
      asin=sin(ang)
      b11=bwx*bwx*(1-acos)+acos
      b12=bwx*bwy*(1-acos)-bwz*asin
      b13=bwx*bwz*(1-acos)+bwy*asin
      b21=bwx*bwy*(1-acos)+bwz*asin
      b22=bwy*bwy*(1-acos)+acos
      b23=bwy*bwz*(1-acos)-bwx*asin
      b31=bwx*bwz*(1-acos)-bwy*asin
      b32=bwy*bwz*(1-acos)+bwx*asin
      b33=bwz*bwz*(1-acos)+acos
***   check n-terminal of fragment--------------->
      ax1=ax0+ax+(ex(n2)-ax)*b11+(ey(n2)-ay)*b12+(ez(n2)-az)*b13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(n2)-ax)*b21+(ey(n2)-ay)*b22+(ez(n2)-az)*b23
      az1=az0+az+(ex(n2)-ax)*b31+(ey(n2)-ay)*b32+(ez(n2)-az)*b33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=1,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=1,n_rot
         ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
         ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
         egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
         egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
         ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
         ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
         ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
         ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
         eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
         ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
         etx_n(i)=ax0+ax+(etx(i)-ax)*a11+(ety(i)-ay)*a12+(etz(i)-az)*a13 !CB
         ety_n(i)=ay0+ay+(etx(i)-ax)*a21+(ety(i)-ay)*a22+(etz(i)-az)*a23
         etz_n(i)=az0+az+(etx(i)-ax)*a31+(ety(i)-ay)*a32+(etz(i)-az)*a33
      enddo
      do i=n_rot+1,n2
         ex_n(i)=ax0+ax+(ex(i)-ax)*b11+(ey(i)-ay)*b12+(ez(i)-az)*b13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*b21+(ey(i)-ay)*b22+(ez(i)-az)*b23
         ez_n(i)=az0+az+(ex(i)-ax)*b31+(ey(i)-ay)*b32+(ez(i)-az)*b33
         egx_n(i)=ax0+ax+(egx(i)-ax)*b11+(egy(i)-ay)*b12+(egz(i)-az)*b13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*b21+(egy(i)-ay)*b22+(egz(i)-az)*b23
         egz_n(i)=az0+az+(egx(i)-ax)*b31+(egy(i)-ay)*b32+(egz(i)-az)*b33
         ecx_n(i)=ecx(i)*b11+ecy(i)*b12+ecz(i)*b13 !cc
         ecy_n(i)=ecx(i)*b21+ecy(i)*b22+ecz(i)*b23
         ecz_n(i)=ecx(i)*b31+ecy(i)*b32+ecz(i)*b33
         ebx_n(i)=ebx(i)*b11+eby(i)*b12+ebz(i)*b13 !bb
         eby_n(i)=ebx(i)*b21+eby(i)*b22+ebz(i)*b23
         ebz_n(i)=ebx(i)*b31+eby(i)*b32+ebz(i)*b33
         etx_n(i)=ax0+ax+(etx(i)-ax)*b11+(ety(i)-ay)*b12+(etz(i)-az)*b13 !CB
         ety_n(i)=ay0+ay+(etx(i)-ax)*b21+(ety(i)-ay)*b22+(etz(i)-az)*b23
         etz_n(i)=az0+az+(etx(i)-ax)*b31+(ety(i)-ay)*b32+(etz(i)-az)*b33
      enddo
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=1,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo

      if(look(1,n))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(1,n,1)+ESHORT(1,n,1) !use ifs
         do kkk=1,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(1,n,-1)+ESHORT(1,n,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(25)=bNa(25)+1
            bNNa(itemp)=bNNa(itemp)+1
            armsd_sum(itemp)=armsd_sum(itemp)+armsd
            N_rmsd(itemp)=N_rmsd(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=1,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=1,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=1,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(25)=bNt(25)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ defo_N finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0)
c     This movement include an off-lattice rotation of fragment
c     and 2 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3, n3.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tran_M
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      parameter(apai=3.1415926)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      m=nfr_i(ifr)-2            !starting point of movement
      n=nfr_f(ifr)+2            !ending point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** translate the fragment **************************
***   
c     d_xyz0=0.5                  !maximum translation coordinates
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
***   check m-terminal of fragment--------------->
      ax1=ax0+ex(m2)
      ay1=ay0+ey(m2)
      az1=az0+ez(m2)
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202
***   check n-terminal of fragment--------------->
 203  ax1=ax0+ex(n2)
      ay1=ay0+ey(n2)
      az1=az0+ez(n2)
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=m2,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=m2,n2
         ex_n(i)=ax0+ex(i)      !CA
         ey_n(i)=ay0+ey(i)
         ez_n(i)=az0+ez(i)
         egx_n(i)=ax0+egx(i)    !SG
         egy_n(i)=ay0+egy(i)
         egz_n(i)=az0+egz(i)
         ecx_n(i)=ecx(i)        !cc
         ecy_n(i)=ecy(i)
         ecz_n(i)=ecz(i)
         ebx_n(i)=ebx(i)        !Hb
         eby_n(i)=eby(i)
         ebz_n(i)=ebz(i)
         etx_n(i)=ax0+etx(i)    !CB
         ety_n(i)=ay0+ety(i)
         etz_n(i)=az0+etz(i)
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=m2,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo
      
      if(look(m,n))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,n,1)+ESHORT(m,n,1) !use ifs
         do kkk=m,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(m,n,-1)+ESHORT(m,n,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(17)=bNa(17)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=m2,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(17)=bNt(17)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ tran_M finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Rotation (angle).
c     This movement include an off-lattice rotation of fragment
c     and 2 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3, n3.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rot_M
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      parameter(apai=3.1415926)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      m=nfr_i(ifr)-2            !starting point of movement
      n=nfr_f(ifr)+2            !ending point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
c     rotation axis (awx,awy,awz):
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
c     rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy(nozy)*2-1)*angle00(ifr) ![-angle00,angle00]
c     rotation points:
      n_rot=m2+int((n2-m2+1)*aranzy(nozy)) ![m2,n2]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202
***   check n-terminal of fragment--------------->
 203  ax1=ax+(ex(n2)-ax)*a11+(ey(n2)-ay)*a12+(ez(n2)-az)*a13
      ay1=ay+(ex(n2)-ax)*a21+(ey(n2)-ay)*a22+(ez(n2)-az)*a23
      az1=az+(ex(n2)-ax)*a31+(ey(n2)-ay)*a32+(ez(n2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=m2,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=m2,n2
         ex_n(i)=ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
         ey_n(i)=ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
         ez_n(i)=az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
         egx_n(i)=ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
         egy_n(i)=ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
         egz_n(i)=az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
         ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
         ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
         ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
         ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
         eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
         ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
         etx_n(i)=ax+(etx(i)-ax)*a11+(ety(i)-ay)*a12+(etz(i)-az)*a13 !CB
         ety_n(i)=ay+(etx(i)-ax)*a21+(ety(i)-ay)*a22+(etz(i)-az)*a23
         etz_n(i)=az+(etx(i)-ax)*a31+(ety(i)-ay)*a32+(etz(i)-az)*a33
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=m2,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo

      if(look(m,n))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,n,1)+ESHORT(m,n,1) !use ifs
         do kkk=m,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(m,n,-1)+ESHORT(m,n,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(20)=bNa(20)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=m2,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=ebx_n(i)  !CB
               ety(i)=eby_n(i)
               etz(i)=ebz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(20)=bNt(20)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ rot_M finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Translation (d_xyz0) + rotation (angle).
c     This movement include an off-lattice rotation of fragment
c     and 2 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3, n3.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trot_M
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      parameter(apai=3.1415926)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      m=nfr_i(ifr)-2            !starting point of movement
      n=nfr_f(ifr)+2            !ending point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5                  !maximum translation coordinates
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
c     rotation axis (awx,awy,awz):
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
c     rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy(nozy)*2-1)*angle00(ifr) ![-angle00,angle00]
c     rotation points:
      n_rot=m2+int((n2-m2+1)*aranzy(nozy)) ![m2,n2]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax0+ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az0+az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202
***   check n-terminal of fragment--------------->
 203  ax1=ax0+ax+(ex(n2)-ax)*a11+(ey(n2)-ay)*a12+(ez(n2)-az)*a13
      ay1=ay0+ay+(ex(n2)-ax)*a21+(ey(n2)-ay)*a22+(ez(n2)-az)*a23
      az1=az0+az+(ex(n2)-ax)*a31+(ey(n2)-ay)*a32+(ez(n2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=m2,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=m2,n2
        ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
        ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
        ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
        egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
        egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
        egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
        ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
        ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
        ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
        ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
        eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
        ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
        etx_n(i)=ax0+ax+(etx(i)-ax)*a11+(ety(i)-ay)*a12+(etz(i)-az)*a13 !CB
        ety_n(i)=ay0+ay+(etx(i)-ax)*a21+(ety(i)-ay)*a22+(etz(i)-az)*a23
        etz_n(i)=az0+az+(etx(i)-ax)*a31+(ety(i)-ay)*a32+(etz(i)-az)*a33
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=m2,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo

      if(look(m,n))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,n,1)+ESHORT(m,n,1) !use ifs
         do kkk=m,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(m,n,-1)+ESHORT(m,n,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(23)=bNa(23)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=m2,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(23)=bNt(23)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ trot_M finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Translation (d_xyz0) + deformation (angle).
c     This movement include an off-lattice rotation of fragment
c     and 2 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3, n3.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defo_M
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      parameter(apai=3.1415926)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/defoangle/defo_angle
      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      m=nfr_i(ifr)-2            !starting point of movement
      n=nfr_f(ifr)+2            !ending point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
      n1=n-1
      n2=n-2                    !ending point of fragment
      n3=n-3
************** rotate the fragment **************************
***   
ccc   d_xyz0=0.5                  !maximum translation coordinates
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
ccc   rotation axis (awx,awy,awz)->N-terminal; (bwx,bwy,bwz)->C-terminal:
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      bwx=asin_theta*cos(aphi)
      bwy=asin_theta*sin(aphi)
      bwz=acos_theta
ccc   rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy(nozy)*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
      bng=(aranzy(nozy)*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
ccc   rotation points:
      n_rot=m2+int((n2-m2+1)*aranzy(nozy)) ![m2,n2]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
ccc   rotation matrix:
      acos=cos(ang)             !for N-terminal of the fragments--->
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
      acos=cos(bng)             !for N-terminal of the fragments--->
      asin=sin(bng)
      b11=bwx*bwx*(1-acos)+acos
      b12=bwx*bwy*(1-acos)-bwz*asin
      b13=bwx*bwz*(1-acos)+bwy*asin
      b21=bwx*bwy*(1-acos)+bwz*asin
      b22=bwy*bwy*(1-acos)+acos
      b23=bwy*bwz*(1-acos)-bwx*asin
      b31=bwx*bwz*(1-acos)-bwy*asin
      b32=bwy*bwz*(1-acos)+bwx*asin
      b33=bwz*bwz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax0+ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az0+az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202
***   check n-terminal of fragment--------------->
 203  ax1=ax0+ax+(ex(n2)-ax)*b11+(ey(n2)-ay)*b12+(ez(n2)-az)*b13
      ay1=ay0+ay+(ex(n2)-ax)*b21+(ey(n2)-ay)*b22+(ez(n2)-az)*b23
      az1=az0+az+(ex(n2)-ax)*b31+(ey(n2)-ay)*b32+(ez(n2)-az)*b33
***   draw distant fragment closer ----->
      r2_bond=(x(n1)-nint(ex(n2)))**2+(y(n1)-nint(ey(n2)))**2
     $     +(z(n1)-nint(ez(n2)))**2
      nbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(n1)-ex(n2))**2+(y(n1)-ey(n2))**2+(z(n1)-ez(n2))**2
 33      t4=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(t4,ica(n)))goto 33
         xn1=x(n)-vx(t4)
         yn1=y(n)-vy(t4)
         zn1=z(n)-vz(t4)
         dis2_new=(xn1-ax1)**2+(yn1-ay1)**2+(zn1-az1)**2
         if(dis2_new.lt.dis2_old)then
 34         t3=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t3,t4))goto 34
            nbig=1
            goto 204            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x(n)-x1
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y(n)-y1
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z(n)-z1
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t3=W21(ijx,ijy,ijz,p)
      t4=W22(ijx,ijy,ijz,p)
      if(.not.goodc(t4,ica(n)))goto 202

c     backup old path ------------>
 204  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      oo(n2)=ica(n2)
      oo(n1)=ica(n1)
      ox(n1)=x(n1)
      oy(n1)=y(n1)
      oz(n1)=z(n1)
      do i=m2,n2
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      nn(n2)=t3
      nn(n1)=t4
      nx(n1)=x(n)-vx(nn(n1))
      ny(n1)=y(n)-vy(nn(n1))
      nz(n1)=z(n)-vz(nn(n1))
      do i=m2,n_rot
         ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
         ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
         egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
         egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
         ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
         ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
         ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
         ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
         eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
         ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
         etx_n(i)=ax0+ax+(etx(i)-ax)*a11+(ety(i)-ay)*a12+(etz(i)-az)*a13 !CB
         ety_n(i)=ay0+ay+(etx(i)-ax)*a21+(ety(i)-ay)*a22+(etz(i)-az)*a23
         etz_n(i)=az0+az+(etx(i)-ax)*a31+(ety(i)-ay)*a32+(etz(i)-az)*a33
      enddo
      do i=n_rot+1,n2
         ex_n(i)=ax0+ax+(ex(i)-ax)*b11+(ey(i)-ay)*b12+(ez(i)-az)*b13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*b21+(ey(i)-ay)*b22+(ez(i)-az)*b23
         ez_n(i)=az0+az+(ex(i)-ax)*b31+(ey(i)-ay)*b32+(ez(i)-az)*b33
         egx_n(i)=ax0+ax+(egx(i)-ax)*b11+(egy(i)-ay)*b12+(egz(i)-az)*b13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*b21+(egy(i)-ay)*b22+(egz(i)-az)*b23
         egz_n(i)=az0+az+(egx(i)-ax)*b31+(egy(i)-ay)*b32+(egz(i)-az)*b33
         ecx_n(i)=ecx(i)*b11+ecy(i)*b12+ecz(i)*b13 !cc
         ecy_n(i)=ecx(i)*b21+ecy(i)*b22+ecz(i)*b23
         ecz_n(i)=ecx(i)*b31+ecy(i)*b32+ecz(i)*b33
         ebx_n(i)=ebx(i)*b11+eby(i)*b12+ebz(i)*b13 !bb
         eby_n(i)=ebx(i)*b21+eby(i)*b22+ebz(i)*b23
         ebz_n(i)=ebx(i)*b31+eby(i)*b32+ebz(i)*b33
         etx_n(i)=ax0+ax+(etx(i)-ax)*b11+(ety(i)-ay)*b12+(etz(i)-az)*b13 !CB
         ety_n(i)=ay0+ay+(etx(i)-ax)*b21+(ety(i)-ay)*b22+(etz(i)-az)*b23
         etz_n(i)=az0+az+(etx(i)-ax)*b31+(ety(i)-ay)*b32+(etz(i)-az)*b33
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202
      d2=(nx(n1)-ex_n(n3))**2+(ny(n1)-ey_n(n3))**2+(nz(n1)-ez_n(n3))**2
      if((d2.lt.23.or.d2.gt.78).and.nbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      ica(n2)=nn(n2)
      ica(n1)=nn(n1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(n1)=nx(n1)
      y(n1)=ny(n1)
      z(n1)=nz(n1)
      do i=m2,n2
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo

      if(look(m,n))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,n,1)+ESHORT(m,n,1) !use ifs
         do kkk=m,n
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(m,n,-1)+ESHORT(m,n,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(26)=bNa(26)+1
            bNNa(itemp)=bNNa(itemp)+1
            armsd_sum(itemp)=armsd_sum(itemp)+armsd
            N_rmsd(itemp)=N_rmsd(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,n
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            ica(n2)=nn(n2)
            ica(n1)=nn(n1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(n1)=nx(n1)
            y(n1)=ny(n1)
            z(n1)=nz(n1)
            do i=m2,n2
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         ica(n2)=oo(n2)
         ica(n1)=oo(n1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(n1)=ox(n1)
         y(n1)=oy(n1)
         z(n1)=oz(n1)
         do i=m2,n2
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(26)=bNt(26)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ defo_M finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0)
c     Move C-terminal of off-lattice fragment
c     This movement include an off-lattice rotation of fragment
c     and 1 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tran_C
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      m=nfr_i(ifr)-2            !starting point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
***   check m-terminal of fragment--------------->
      ax1=ax0+ex(m2)
      ay1=ay0+ey(m2)
      az1=az0+ez(m2)
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202

c     backup old path ------------>
 203  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      do i=m2,Lch
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      do i=m2,Lch
         ex_n(i)=ax0+ex(i)      !CA
         ey_n(i)=ay0+ey(i)
         ez_n(i)=az0+ez(i)
         egx_n(i)=ax0+egx(i)    !SG
         egy_n(i)=ay0+egy(i)
         egz_n(i)=az0+egz(i)
         ecx_n(i)=ecx(i)        !cc
         ecy_n(i)=ecy(i)
         ecz_n(i)=ecz(i)
         ebx_n(i)=ebx(i)        !Hb
         eby_n(i)=eby(i)
         ebz_n(i)=ebz(i)
         etx_n(i)=ax0+etx(i)    !CB
         ety_n(i)=ay0+ety(i)
         etz_n(i)=az0+etz(i)
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      do i=m2,Lch
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo

      if(look(m,Lch))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,Lch,1)+ESHORT(m,Lch,1) !use ifs
         do kkk=m,Lch
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(m,Lch,-1)+ESHORT(m,Lch,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(18)=bNa(18)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,Lch
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            do i=m2,Lch
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(18)=bNt(18)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ tran_C finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0) + rotation (angle0)
c     Move C-terminal of off-lattice fragment
c     This movement include an off-lattice rotation of fragment
c     and 1 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine trot_C
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      m=nfr_i(ifr)-2            !starting point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
************** rotate the fragment **************************
***   
c     d_xyz0=0.5               !range of translation movement
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
c     rotation axis (awx,awy,awz):
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
c     rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy(nozy)*2-1)*angle00(ifr) ![-angle00,angle00]
c     rotation points:
      ax=ex(m2)
      ay=ey(m2)
      az=ez(m2)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax0+ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az0+az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202

c     backup old path ------------>
 203  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      do i=m2,Lch
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      do i=m2,Lch
        ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
        ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
        ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
        egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
        egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
        egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
        ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
        ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
        ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
        ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
        eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
        ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
        etx_n(i)=ax0+ax+(etx(i)-ax)*a11+(ety(i)-ay)*a12+(etz(i)-az)*a13 !CB
        ety_n(i)=ay0+ay+(etx(i)-ax)*a21+(ety(i)-ay)*a22+(etz(i)-az)*a23
        etz_n(i)=az0+az+(etx(i)-ax)*a31+(ety(i)-ay)*a32+(etz(i)-az)*a33
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      do i=m2,Lch
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo

      if(look(m,Lch))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,Lch,1)+ESHORT(m,Lch,1) !use ifs
         do kkk=m,Lch
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(m,Lch,-1)+ESHORT(m,Lch,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(24)=bNa(24)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,Lch
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            do i=m2,Lch
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(24)=bNt(24)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ rot_C finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translation (d_xyz0) + deformation (angle0)
c     Move C-terminal of off-lattice fragment
c     This movement include an off-lattice rotation of fragment
c     and 1 two-bond movement in each end of fragment.
c     The minimum length of fragment is 2, because of m3
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine defo_C
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)       !number of replicas
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/frozen/L_cut,d_xyz0,d_xyz00(ndim),angle0,angle00(ndim)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      dimension ex_o(ndim),ey_o(ndim),ez_o(ndim) !CA
      dimension egx_o(ndim),egy_o(ndim),egz_o(ndim) !SG
      dimension ecx_o(ndim),ecy_o(ndim),ecz_o(ndim) !cc
      dimension ebx_o(ndim),eby_o(ndim),ebz_o(ndim) !Hb
      dimension etx_o(ndim),ety_o(ndim),etz_o(ndim) !CB
      dimension ex_n(ndim),ey_n(ndim),ez_n(ndim) !CA
      dimension egx_n(ndim),egy_n(ndim),egz_n(ndim) !SG
      dimension ecx_n(ndim),ecy_n(ndim),ecz_n(ndim) !cc
      dimension ebx_n(ndim),eby_n(ndim),ebz_n(ndim) !Hb
      dimension etx_n(ndim),ety_n(ndim),etz_n(ndim) !CB

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/rotat/ar0(ndim),athita0(ndim),aphi0(ndim)

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/w2a/w21(-10:10,-10:10,-10:10,80)
      common/w2b/w22(-10:10,-10:10,-10:10,80)
      common/w2c/Nw(-10:10,-10:10,-10:10)

      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/echain1/ex(ndim),ey(ndim),ez(ndim) !CA
      common/echain2/egx(ndim),egy(ndim),egz(ndim) !SG
      common/echain4/ecx(ndim),ecy(ndim),ecz(ndim) !cc
      common/echain5/ebx(ndim),eby(ndim),ebz(ndim) !Hb
      common/echain6/etx(ndim),ety(ndim),etz(ndim) !CB
      common/chainm/mv(ndim)
      common/defoangle/defo_angle
      common/rmsd_defo/armsd,armsd_sum(100),N_rmsd(100)
      common/ranzy/nozy
      common/lattice/m_latt,latt1,latt2

      m=nfr_i(ifr)-2            !starting point of movement
      m1=m+1
      m2=m+2                    !starting point of fragment
      m3=m+3
************** rotate the fragment **************************
***   
ccc   d_xyz0=0.5               !range of translation movement
      ax0=(aranzy(nozy)*2-1)*d_xyz00(ifr) !translation of coordinates, [-det,det]
      ay0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
      az0=(aranzy(nozy)*2-1)*d_xyz00(ifr)
ccc   rotation axis (awx,awy,awz)->N-terminal; (bwx,bwy,bwz)->C-terminal:
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      awx=asin_theta*cos(aphi)
      awy=asin_theta*sin(aphi)
      awz=acos_theta
      acos_theta=1-2.*aranzy(nozy)
      asin_theta=sqrt(1.0-acos_theta*acos_theta)
      aphi=2.*3.1415926*aranzy(nozy) ![0,2*pai]
      bwx=asin_theta*cos(aphi)
      bwy=asin_theta*sin(aphi)
      bwz=acos_theta
ccc   rotation angle           !in arc-angle, 1=57o, rotation angle.
      ang=(aranzy(nozy)*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
      bng=(aranzy(nozy)*2-1)*angle00(ifr)/defo_angle ![-angle00/d,angle00/d]
ccc   rotation points:
      n_rot=m2+int((Lch-m2+1)*aranzy(nozy)) ![m2,Lch]
      ax=ex(n_rot)
      ay=ey(n_rot)
      az=ez(n_rot)
***   rotation matrix:
      acos=cos(ang)
      asin=sin(ang)
      a11=awx*awx*(1-acos)+acos
      a12=awx*awy*(1-acos)-awz*asin
      a13=awx*awz*(1-acos)+awy*asin
      a21=awx*awy*(1-acos)+awz*asin
      a22=awy*awy*(1-acos)+acos
      a23=awy*awz*(1-acos)-awx*asin
      a31=awx*awz*(1-acos)-awy*asin
      a32=awy*awz*(1-acos)+awx*asin
      a33=awz*awz*(1-acos)+acos
      acos=cos(ang)
      asin=sin(ang)
      b11=bwx*bwx*(1-acos)+acos
      b12=bwx*bwy*(1-acos)-bwz*asin
      b13=bwx*bwz*(1-acos)+bwy*asin
      b21=bwx*bwy*(1-acos)+bwz*asin
      b22=bwy*bwy*(1-acos)+acos
      b23=bwy*bwz*(1-acos)-bwx*asin
      b31=bwx*bwz*(1-acos)-bwy*asin
      b32=bwy*bwz*(1-acos)+bwx*asin
      b33=bwz*bwz*(1-acos)+acos
***   check m-terminal of fragment--------------->
      ax1=ax0+ax+(ex(m2)-ax)*a11+(ey(m2)-ay)*a12+(ez(m2)-az)*a13 !v1=v0+(a)v
      ay1=ay0+ay+(ex(m2)-ax)*a21+(ey(m2)-ay)*a22+(ez(m2)-az)*a23
      az1=az0+az+(ex(m2)-ax)*a31+(ey(m2)-ay)*a32+(ez(m2)-az)*a33
***   draw distant fragment closer ----->
      r2_bond=(x(m1)-nint(ex(m2)))**2+(y(m1)-nint(ey(m2)))**2
     $     +(z(m1)-nint(ez(m2)))**2
      mbig=0
      if(r2_bond.gt.latt2)then
         dis2_old=(x(m1)-ex(m2))**2+(y(m1)-ey(m2))**2+(z(m1)-ez(m2))**2
 31      t1=int(aranzy(nozy)*anvec)+1 !random vector
         if(.not.goodc(ica(m-1),t1))goto 31
         xm1=x(m)+vx(t1)
         ym1=y(m)+vy(t1)
         zm1=z(m)+vz(t1)
         dis2_new=(xm1-ax1)**2+(ym1-ay1)**2+(zm1-az1)**2
         if(dis2_new.lt.dis2_old)then
 32         t2=int(aranzy(nozy)*anvec)+1 !random vector
            if(.not.goodc(t1,t2))goto 32
            mbig=1
            goto 203            !take t1, t2
         endif
      endif
***
      x1=nint(ax1)
      y1=nint(ay1)
      z1=nint(az1)
      ijx=x1-x(m)
      if(abs(ijx).gt.10)goto 202 !without correct movement, quit
      ijy=y1-y(m)
      if(abs(ijy).gt.10)goto 202 !without correct movement, quit
      ijz=z1-z(m)
      if(abs(ijz).gt.10)goto 202 !without correct movement, quit
      nc=Nw(ijx,ijy,ijz)
      if(nc.lt.1)goto 202       !without correct movement, quit
      p=int(aranzy(nozy)*nc)+1    ![1,nc]
      t1=W21(ijx,ijy,ijz,p)
      t2=W22(ijx,ijy,ijz,p)
      if(.not.goodc(ica(m-1),t1))goto 202

c     backup old path ------------>
 203  oo(m)=ica(m)
      oo(m1)=ica(m1)
      ox(m1)=x(m1)
      oy(m1)=y(m1)
      oz(m1)=z(m1)
      do i=m2,Lch
         ex_o(i)=ex(i)          !CA
         ey_o(i)=ey(i)
         ez_o(i)=ez(i)
         egx_o(i)=egx(i)        !SG
         egy_o(i)=egy(i)
         egz_o(i)=egz(i)
         ecx_o(i)=ecx(i)        !cc
         ecy_o(i)=ecy(i)
         ecz_o(i)=ecz(i)
         ebx_o(i)=ebx(i)        !Hb
         eby_o(i)=eby(i)
         ebz_o(i)=ebz(i)
         etx_o(i)=etx(i)        !CB
         ety_o(i)=ety(i)
         etz_o(i)=etz(i)
      enddo

c     prepare the new path------------>
      nn(m)=t1
      nn(m1)=t2
      nx(m1)=x(m)+vx(nn(m))
      ny(m1)=y(m)+vy(nn(m))
      nz(m1)=z(m)+vz(nn(m))
      do i=m2,n_rot
         ex_n(i)=ax0+ax+(ex(i)-ax)*a11+(ey(i)-ay)*a12+(ez(i)-az)*a13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*a21+(ey(i)-ay)*a22+(ez(i)-az)*a23
         ez_n(i)=az0+az+(ex(i)-ax)*a31+(ey(i)-ay)*a32+(ez(i)-az)*a33
         egx_n(i)=ax0+ax+(egx(i)-ax)*a11+(egy(i)-ay)*a12+(egz(i)-az)*a13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*a21+(egy(i)-ay)*a22+(egz(i)-az)*a23
         egz_n(i)=az0+az+(egx(i)-ax)*a31+(egy(i)-ay)*a32+(egz(i)-az)*a33
         ecx_n(i)=ecx(i)*a11+ecy(i)*a12+ecz(i)*a13 !cc
         ecy_n(i)=ecx(i)*a21+ecy(i)*a22+ecz(i)*a23
         ecz_n(i)=ecx(i)*a31+ecy(i)*a32+ecz(i)*a33
         ebx_n(i)=ebx(i)*a11+eby(i)*a12+ebz(i)*a13 !bb
         eby_n(i)=ebx(i)*a21+eby(i)*a22+ebz(i)*a23
         ebz_n(i)=ebx(i)*a31+eby(i)*a32+ebz(i)*a33
         etx_n(i)=ax0+ax+(etx(i)-ax)*a11+(ety(i)-ay)*a12+(etz(i)-az)*a13 !CB
         ety_n(i)=ay0+ay+(etx(i)-ax)*a21+(ety(i)-ay)*a22+(etz(i)-az)*a23
         etz_n(i)=az0+az+(etx(i)-ax)*a31+(ety(i)-ay)*a32+(etz(i)-az)*a33
      enddo
      do i=n_rot+1,Lch
         ex_n(i)=ax0+ax+(ex(i)-ax)*b11+(ey(i)-ay)*b12+(ez(i)-az)*b13 !CA
         ey_n(i)=ay0+ay+(ex(i)-ax)*b21+(ey(i)-ay)*b22+(ez(i)-az)*b23
         ez_n(i)=az0+az+(ex(i)-ax)*b31+(ey(i)-ay)*b32+(ez(i)-az)*b33
         egx_n(i)=ax0+ax+(egx(i)-ax)*b11+(egy(i)-ay)*b12+(egz(i)-az)*b13 !SG
         egy_n(i)=ay0+ay+(egx(i)-ax)*b21+(egy(i)-ay)*b22+(egz(i)-az)*b23
         egz_n(i)=az0+az+(egx(i)-ax)*b31+(egy(i)-ay)*b32+(egz(i)-az)*b33
         ecx_n(i)=ecx(i)*b11+ecy(i)*b12+ecz(i)*b13 !cc
         ecy_n(i)=ecx(i)*b21+ecy(i)*b22+ecz(i)*b23
         ecz_n(i)=ecx(i)*b31+ecy(i)*b32+ecz(i)*b33
         ebx_n(i)=ebx(i)*b11+eby(i)*b12+ebz(i)*b13 !bb
         eby_n(i)=ebx(i)*b21+eby(i)*b22+ebz(i)*b23
         ebz_n(i)=ebx(i)*b31+eby(i)*b32+ebz(i)*b33
         etx_n(i)=ax0+ax+(etx(i)-ax)*b11+(ety(i)-ay)*b12+(etz(i)-az)*b13 !CB
         ety_n(i)=ay0+ay+(etx(i)-ax)*b21+(ety(i)-ay)*b22+(etz(i)-az)*b23
         etz_n(i)=az0+az+(etx(i)-ax)*b31+(ety(i)-ay)*b32+(etz(i)-az)*b33
      enddo
      d2=(nx(m1)-ex_n(m3))**2+(ny(m1)-ey_n(m3))**2+(nz(m1)-ez_n(m3))**2
      if((d2.lt.23.or.d2.gt.78).and.mbig.eq.0)goto 202

c     change conformation to new path -------------->
      ica(m)=nn(m)
      ica(m1)=nn(m1)
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      x(m1)=nx(m1)
      y(m1)=ny(m1)
      z(m1)=nz(m1)
      do i=m2,Lch
         ex(i)=ex_n(i)          !CA
         ey(i)=ey_n(i)
         ez(i)=ez_n(i)
         egx(i)=egx_n(i)        !SG
         egy(i)=egy_n(i)
         egz(i)=egz_n(i)
         ecx(i)=ecx_n(i)        !cc
         ecy(i)=ecy_n(i)
         ecz(i)=ecz_n(i)
         ebx(i)=ebx_n(i)        !Hb
         eby(i)=eby_n(i)
         ebz(i)=ebz_n(i)
         etx(i)=etx_n(i)        !CB
         ety(i)=ety_n(i)
         etz(i)=etz_n(i)
      enddo

      if(look(m,Lch))then         ! check excluded volumn for passage of [m,n]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m,Lch,1)+ESHORT(m,Lch,1) !use ifs
         do kkk=m,Lch
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
         Eold=EHB(m,Lch,-1)+ESHORT(m,Lch,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(27)=bNa(27)+1
            bNNa(itemp)=bNNa(itemp)+1
            armsd_sum(itemp)=armsd_sum(itemp)+armsd
            N_rmsd(itemp)=N_rmsd(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m,Lch
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c     energ=energ+de	
            ica(m)=nn(m)
            ica(m1)=nn(m1)
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            x(m1)=nx(m1)
            y(m1)=ny(m1)
            z(m1)=nz(m1)
            do i=m2,Lch
               ex(i)=ex_n(i)    !CA
               ey(i)=ey_n(i)
               ez(i)=ez_n(i)
               egx(i)=egx_n(i)  !SG
               egy(i)=egy_n(i)
               egz(i)=egz_n(i)
               ecx(i)=ecx_n(i)  !cc
               ecy(i)=ecy_n(i)
               ecz(i)=ecz_n(i)
               ebx(i)=ebx_n(i)  !Hb
               eby(i)=eby_n(i)
               ebz(i)=ebz_n(i)
               etx(i)=etx_n(i)  !CB
               ety(i)=ety_n(i)
               etz(i)=etz_n(i)
            enddo
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(m)=oo(m)
         ica(m1)=oo(m1)
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         x(m1)=ox(m1)
         y(m1)=oy(m1)
         z(m1)=oz(m1)
         do i=m2,Lch
            ex(i)=ex_o(i)       !CA
            ey(i)=ey_o(i)
            ez(i)=ez_o(i)
            egx(i)=egx_o(i)     !SG
            egy(i)=egy_o(i)
            egz(i)=egz_o(i)
            ecx(i)=ecx_o(i)     !cc
            ecy(i)=ecy_o(i)
            ecz(i)=ecz_o(i)
            ebx(i)=ebx_o(i)     !Hb
            eby(i)=eby_o(i)
            ebz(i)=ebz_o(i)
            etx(i)=etx_o(i)     !CB
            ety(i)=ety_o(i)
            etz(i)=etz_o(i)
         enddo
      endif                     !for look(i,m)

 202  continue
      bNt(27)=bNt(27)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ rot_C finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     random walk in either N-terminal 
c     m residues of '1, 2, ..., m' are relocated by random walk.
c     Maximum number of moving point = Mend_N
c     old ---> new ---> old ----> new
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move_n_end
      implicit integer(i-z)
      parameter(ndim=1999)
                parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      double precision bNa,bNt,bNNa,bNNt
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/ranzy/nozy

      m=int(aranzy(nozy)*(Mend_N-1))+1 !m=[1,Mend_N-1], [1,m] will be moved.
      
c     back-up old path---------->
      do i=1,m
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
c     prepare new path and move to new path -------------->
      do i=m,1,-1
 111     iv=int(aranzy(nozy)*anvec)+1
         if(.not.goodc(iv,ica(i+1))) goto 111
         x(i)=x(i+1)-vx(iv)
         y(i)=y(i+1)-vy(iv)
         z(i)=z(i+1)-vz(iv)
         ica(i)=iv
         nx(i)=x(i)             !memory of new conformation
         ny(i)=y(i)             !memory of new conformation
         nz(i)=z(i)             !memory of new conformation
         nn(i)=ica(i)           !memory of new conformation
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      
      if(look(1,m+1))then       ! check excluded volumn for passage of [2,m+1]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(1,m+1,1)+ESHORT(1,m+1,1) !use ifs
         do kkk=1,m+1
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo
         
c     return back the conformation and calculate E_old --------->
         do i=1,m
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         
         Eold=EHB(1,m+1,-1)+ESHORT(1,m+1,-1) !cut off old nop by istat=-1

c     calculate eprofn while dord was calculated when call EHB(i,m,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(11)=bNa(11)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
c            energ=energ+de
            do i=1,m
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            do kkk=1,m+1
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
         do i=1,m
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
113   continue
      bNt(11)=bNt(11)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move_n_end finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     random walk in either C-terminal 
c     Maximum number of moving point = Mend_C
c     old ---> new ---> old ----> new
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move_c_end
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      double precision bNa,bNt,bNNa,bNNt
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/ranzy/nozy

      m=int(aranzy(nozy)*(Mend_C-1))+1 !m=[1,Mend_C-1], m points will be moved
      mm=Lch-m+1                !(mm,Lch] will be relocated
c     mm-1 should be moveable residue, ica(m-2) exists.

c     backup old path----------->
      do i=mm,Lch
         ox(i)=x(i)             !memory of old conformation
         oy(i)=y(i)             !memory of old conformation
         oz(i)=z(i)             !memory of old conformation
         oo(i-1)=ica(i-1)       !memory of old conformation
      enddo
c     prepare new path------------->
      do i=mm,Lch
 111     iv=int(aranzy(nozy)*anvec)+1
         if(.not.goodc(ica(i-2),iv)) goto 111
         x(i)=x(i-1)+vx(iv)
         y(i)=y(i-1)+vy(iv)
         z(i)=z(i-1)+vz(iv)
         ica(i-1)=iv
         nx(i)=x(i)             !memory of new conformation
         ny(i)=y(i)             !memory of new conformation
         nz(i)=z(i)             !memory of new conformation
         nn(i-1)=ica(i-1)       !memory of new conformation
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      
      if(look(mm-1,Lch))then    ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=1,Lch
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(mm-1,Lch,1)+ESHORT(mm-1,Lch,1) !use ifs
         do kkk=mm-1,Lch
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         do i=mm,Lch
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i-1)=oo(i-1)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(mm-1,Lch,-1)+ESHORT(mm-1,Lch,-1) !nop go to new

c     calculate eprofn while dord was calculated when call EHB(i,m,1)--->
         eprofn=0.0
         do pp=1,Lch
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(12)=bNa(12)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=1,Lch
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
c     energ=energ+de
            do i=mm,Lch
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i-1)=nn(i-1)
            enddo
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)

            do kkk=mm-1,Lch
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
         do i=mm,Lch
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i-1)=oo(i-1)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
113   continue
      bNt(12)=bNt(12)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move_c_end finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Sequence down-shift, i12<i21
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move7a
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/ranzy/nozy

      mm=12

      do i=1,lenf
         i21=4+int(aranzy(nozy)*(Lch-6)) ![4,lenf-3]
         if(u21(ica(i21),ica(i21+1)).ne.0)goto 31
      enddo
      goto 33                   !can not find a two-bond can merg into one-bond

 31   i0=max(2,i21-mm)
      do i=1,5
         i12=i0+int(aranzy(nozy)*((i21-i0)-1.0001)) ![i0,i21-2]
         nc=m12(ica(i12))       !number of ways of extending 1-->2
         if(nc.gt.0)goto 32
      enddo
      goto 33                   !cannot find a one-bond can extend to two-bond

c     ------------- 2 -> 1 ----------------------------
 32   nn(i21+1)=u21(ica(i21),ica(i21+1))
      if(.not.goodc(ica(i21-1),nn(i21+1)))goto33
      if(.not.goodc(nn(i21+1),ica(i21+2)))goto33
c     ------------- 1 -> 2 ----------------------------
      p=int(aranzy(nozy)*nc)+1
      nn(i12)=u1(ica(i12),p)
      nn(i12+1)=u2(ica(i12),p)
      if(.not.goodc(ica(i12-1),nn(i12)))goto33
      if(.not.goodc(nn(i12+1),ica(i12+1)))goto33

c     backup old path ------------>
      oo(i12)=ica(i12)
      do i=i12+1,i21+1
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
c     prepare new path ------------>
      do i=i12+2,i21
         nn(i)=ica(i-1)
      enddo
      nx(i12+1)=x(i12)+vx(nn(i12))
      ny(i12+1)=y(i12)+vy(nn(i12))
      nz(i12+1)=z(i12)+vz(nn(i12))
      do i=i12+2,i21+1
         nx(i)=nx(i-1)+vx(nn(i-1))
         ny(i)=ny(i-1)+vy(nn(i-1))
         nz(i)=nz(i-1)+vz(nn(i-1))
      enddo

c     change conformation to new path -------------->
      ica(i12)=nn(i12)
      do i=i12+1,i21+1
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
         ica(i)=nn(i)
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      m1=i12       !fixed
      m2=i21+2     !fixed
      if(look(m1,m2))then       ! check excluded volumn for passage of [m1,m2]
c     calculate E_new--------------->
         do pp=2,lenf1
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m1,m2,1)+ESHORT(m1,m2,1) !use ifs
         do kkk=m1,m2
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(i12)=oo(i12)
         do i=i12+1,i21+1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m1,m2,-1)+ESHORT(m1,m2,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=2,lenf1
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(13)=bNa(13)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m1,m2
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            ica(i12)=nn(i12)
            do i=i12+1,i21+1
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(i12)=oo(i12)
         do i=i12+1,i21+1
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)

 33   continue
      bNt(13)=bNt(13)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move7a finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Sequence up-shift, i12>i21
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move7b
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/ranzy/nozy

      mm=12

      do i=1,lenf
         i21=2+int(aranzy(nozy)*(Lch-6)) ![2,lenf-5]
         if(u21(ica(i21),ica(i21+1)).ne.0)goto 31
      enddo
      goto 33                   !can not find a two-bond can merg into one-bond

 31   i0=min(lenf2,i21+mm)
      do i=1,5
         i12=i21+3+int(aranzy(nozy)*(i0-i21-2.0001)) ![i21+3,i0]
         nc=m12(ica(i12))       !number of ways of extending 1-->2
         if(nc.gt.0)goto 32
      enddo
      goto 33                !can not find a one-bond can extend into two-bond

c     ------------- 2 -> 1 ----------------------------
 32   nn(i21)=u21(ica(i21),ica(i21+1))
      if(.not.goodc(ica(i21-1),nn(i21)))goto33
      if(.not.goodc(nn(i21),ica(i21+2)))goto33
c     ------------- 1 -> 2 ----------------------------
      p=int(aranzy(nozy)*nc)+1
      nn(i12-1)=u1(ica(i12),p)
      nn(i12)=u2(ica(i12),p)
      if(.not.goodc(ica(i12-1),nn(i12-1)))goto33
      if(.not.goodc(nn(i12),ica(i12+1)))goto33

c     backup old path ------------>
      oo(i21)=ica(i21)
      do i=i21+1,i12
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo
c     prepare new path ------------>
      do i=i21+1,i12-2
         nn(i)=ica(i+1)
      enddo
      nx(i21+1)=x(i21)+vx(nn(i21))
      ny(i21+1)=y(i21)+vy(nn(i21))
      nz(i21+1)=z(i21)+vz(nn(i21))
      do i=i21+2,i12
         nx(i)=nx(i-1)+vx(nn(i-1))
         ny(i)=ny(i-1)+vy(nn(i-1))
         nz(i)=nz(i-1)+vz(nn(i-1))
      enddo

c     change conformation to new path -------------->
      ica(i21)=nn(i21)
      do i=i21+1,i12
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
         ica(i)=nn(i)
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)

      m1=i21       !fixed
      m2=i12+1     !fixed
      if(look(m1,m2))then       ! check excluded volumn for passage of [m1,m2]
c     calculate E_new--------------->
         do pp=2,lenf1
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(m1,m2,1)+ESHORT(m1,m2,1) !use ifs
         do kkk=m1,m2
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         ica(i21)=oo(i21)
         do i=i21+1,i12
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m1,m2,-1)+ESHORT(m1,m2,-1)

c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=2,lenf1
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c     E=energ
c     Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(14)=bNa(14)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m1,m2
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            ica(i21)=nn(i21)
            do i=i21+1,i12
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         ica(i21)=oo(i21)
         do i=i21+1,i12
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)

 33   continue
      bNt(14)=bNt(14)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move7b finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  q-bond random walk movement:  (Mran1 <= q <= Mran2)
c  residues of 'i,i+1, ..., i+q,i+q+1' are involved, 
c  but positions bond-vectors of i+1, i+2, ..., i+q are changed.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move8
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/maxdis2/maxdis2(ndim)
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)
      common/shift/u21(nvec,nvec),m12(nvec),u1(nvec,100),u2(nvec,100)

      common/move31/v31(-15:15,-15:15,-15:15,6)
      common/move32/v32(-15:15,-15:15,-15:15,6)
      common/move33/v33(-15:15,-15:15,-15:15,6)
      common/move34/Np3(-15:15,-15:15,-15:15)

      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)
      common/ranzy/nozy
      common/nwmax1/nvecr

c     q bonds will be moved, q is in [Mran1,Mran2]
      Mran1=4
      M_ran2=6
      q=int(aranzy(nozy)*(Mran2-Mran1+1))+Mran1
      q3=q-3

c     choose the position (m1) to be moved
c     m1>1, m2<lenf, otherwise can not do goodc()!
c     [m1+1,m2-1] will be really moved ---------------->
 111  m1=int(aranzy(nozy)*(lenf-q-2.00001))+2 !m1=[2,lenf-q-1]
      m2=m1+q                   !m2<=lenf-1

      mx=x(m2)-x(m1)
      my=y(m2)-y(m1)
      mz=z(m2)-z(m1)
cccccccccccccccFirst q-3 bonds cccccccccccccccccccccccccc
      nn(m1-1)=ica(m1-1)
      do i=1,q3
         uu1=0
 112     uu=1+int(nvecr*aranzy(nozy)) ! uu in [1,nvec]
         uu1=uu1+1
         if(uu1.gt.nvecr)goto 111 ! path maybe wrong!
         xx=mx-vx(uu)
         yy=my-vy(uu)
         zz=mz-vz(uu)

         xx2=xx*xx
         if(xx2.ge.maxdis2(q-i))goto 112
         yy2=yy*yy
         if(yy2.ge.maxdis2(q-i))goto 112
         rr=xx2+yy2+zz*zz   ! rr=xx*xx+yy*yy+zz*zz
         if(rr.ge.maxdis2(q-i))goto 112 !too distant for remaining walks.
         if(rr.lt.12)goto 112 !too close for overlap
         if(.not.goodc(uu,nn(m1+i-2)))goto 112 !check neighbor

         nn(m1+i-1)=uu
         mx=xx
         my=yy
         mz=zz
      enddo
      if(mx.gt.12.or.mx.lt.-12) goto 111 !no defined vector in 3-bond move
      if(my.gt.12.or.my.lt.-12) goto 111 !no defined vector in 3-bond move
      if(mz.gt.12.or.mz.lt.-12) goto 111 !no defined vector in 3-bond move
      
ccccccccccccccccccc Last 3 bonds ccccccccccccccccccccccccc
      nc=Np3(mx,my,mz)          !number of pathes from i to i+m ???
      if(nc.eq.0)then
         write(*,*)'absent q-movement in move5',q,mx,my,mz,m1,m2
         goto 114
      endif
      uu1=0
 113  p=int(aranzy(nozy)*(nc-0.00001))+1    ![1,nc]
      uu1=uu1+1
      if(uu1.gt.nc)goto 111
      nn(m2-3)=v31(mx,my,mz,p)
      if(.not.goodc(nn(m2-4),nn(m2-3))) goto 113
      nn(m2-1)=v33(mx,my,mz,p)  !???
      if(.not.goodc(nn(m2-1),ica(m2))) goto 113
      nn(m2-2)=v32(mx,my,mz,p)  !???

c     backup old path ------------>
      do i=1,q
         ox(m1+i)=x(m1+i)
         oy(m1+i)=y(m1+i)
         oz(m1+i)=z(m1+i)
         oo(m1+i-1)=ica(m1+i-1)
      enddo
c     prepare new path ------------>
      nx(m1+1)=x(m1)+vx(nn(m1))
      ny(m1+1)=y(m1)+vy(nn(m1))
      nz(m1+1)=z(m1)+vz(nn(m1))
      do i=2,q
         nx(m1+i)=nx(m1+i-1)+vx(nn(m1+i-1))
         ny(m1+i)=ny(m1+i-1)+vy(nn(m1+i-1))
         nz(m1+i)=nz(m1+i-1)+vz(nn(m1+i-1))
      enddo

c     change conformation to new path -------------->
      do i=1,q
         x(m1+i)=nx(m1+i)
         y(m1+i)=ny(m1+i)
         z(m1+i)=nz(m1+i)
         ica(m1+i-1)=nn(m1+i-1)
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      
      if(look(m1,m2))then       ! check excluded volumn for passage of [i,m]
c     calculate E_new--------------->
         do pp=2,lenf1
            nop(pp)=nopp(pp)    !prepare to calculate energy
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
c            nhbn(pp)=nhbnn(pp)
         enddo
         Enew=EHB(m1,m2,1)+ESHORT(m1,m2,1) !use ifs
         do kkk=m1,m2
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         do i=1,q
            x(m1+i)=ox(m1+i)
            y(m1+i)=oy(m1+i)
            z(m1+i)=oz(m1+i)
            ica(m1+i-1)=oo(m1+i-1)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(m1,m2,-1)+ESHORT(m1,m2,-1)

c     calculate eprofn while dord was calculated when call EHB(i,m,1)--->
         eprofn=0.0
         do pp=2,lenf1
            is=seq(pp)
            ia=noa(pp)        !noa is now 
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo

c     Metropolis ------------------>
         dE=Enew-Eold+dord+eprofn-eprofo
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.3)then        !rejected
            icnt=icnto
            sumct=sumcto 		
         else                   !accepted
            bNa(7)=bNa(7)+1
            bN5a(q)=bN5a(q)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
c               nhbnn(pp)=nhbn(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=m1,m2
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            do i=1,q
               x(m1+i)=nx(m1+i)
               y(m1+i)=ny(m1+i)
               z(m1+i)=nz(m1+i)
               ica(m1+i-1)=nn(m1+i-1)
            enddo
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center       !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         do i=1,q
            x(m1+i)=ox(m1+i)
            y(m1+i)=oy(m1+i)
            z(m1+i)=oz(m1+i)
            ica(m1+i-1)=oo(m1+i-1)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif                     !for look(i,m)
 114  continue
      bNt(7)=bNt(7)+1
      bN5t(q)=bN5t(q)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move8 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     globe movement, there is biase.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine move9
      implicit integer(i-z)
      parameter(nrep=100)       !number of replicas
      parameter(ndim=1999)
      parameter(nvec=416)
      common/logica/goodc
      logical look, goodc(nvec,nvec)
      common/temperature/itemp,atemp
      common/lengths/Lch,Lch1,Lch2
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)
      common/shape/amx,amy,amz,afs(ndim),afsn(ndim)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/seqe/seq(ndim),sec(ndim)
      common/moverange/Mran1,Mran2,Mend,Mend_N,Mend_C
      common/naccept/bNa(100),bNt(100),bNNa(100),bNNt(100)
      common/naccept1/bN5a(100),bN5t(100)
      double precision bNa,bNt,bNNa,bNNt,bN5a,bN5t
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev
      common/commonuse4/al2,al4,al5,al6,al7,al8,anvec
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1

      common/mov2/v21(nvec,nvec,80),v22(nvec,nvec,80),Np2(nvec,nvec)

      dimension fax(ndim),fay(ndim),faz(ndim)
      dimension ox(ndim),oy(ndim),oz(ndim),oo(ndim)
      dimension nx(ndim),ny(ndim),nz(ndim),nn(ndim)

      common/gapmove2/ras2(ndim),nfl2
      common/gapmove3/ras3(ndim),nfl3
      common/gapmove4/ras4(ndim),nfl4
      common/gapmove5/ras5(ndim),nfl5
      common/gapmove6/ras6(ndim),nfl6
      common/ranzy/nozy
      common/nwmax1/nvecr

******backup old conformation --------------->
      do i=1,lenf
         ox(i)=x(i)
         oy(i)=y(i)
         oz(i)=z(i)
         oo(i)=ica(i)
      enddo

******prepare new confromation ---------------->
c     (x,y,z)->(fax,fay,faz) based on [x(k),y[k],z[k]]:
      k=int(lenf*aranzy(nozy))+1   ![1,lenf]
      ddxyz=1.5
 99   fdx=(aranzy(nozy)*ddxyz*2.0)-ddxyz ![-1.5, 1.5]
      fdy=(aranzy(nozy)*ddxyz*2.0)-ddxyz
      fdz=(aranzy(nozy)*ddxyz*2.0)-ddxyz
      ar=fdx*fdx+fdy*fdy+fdz*fdz
      if(ar.lt.0.75) go to 99
      do i=1,lenf
         ar=sqrt(float((x(k)-x(i))**2+(y(k)-y(i))**2+(z(k)-z(i))**2))
         fax(i)=x(i)+fdx*(1.0-ar/(1.5*acrit))
         fay(i)=y(i)+fdy*(1.0-ar/(1.5*acrit))
         faz(i)=z(i)+fdz*(1.0-ar/(1.5*acrit))
      enddo

******project the conformation onto lattices, i.e. (fax,fay,faz)->(nx,ny,nz):
      px=nint(fax(1))
      py=nint(fay(1))
      pz=nint(faz(1))
      nx(1)=px
      ny(1)=py
      nz(1)=pz
      DO 101 i=2,lenf
         armin=10000.
         kk=0
         do 1009 k=1,nvecr
            if(i.ne.2) then
               if(.not.goodc(kkkk,k))goto 1009
            endif
            kx=px+vx(k)
            ky=py+vy(k)
            kz=pz+vz(k)
            bx=fax(i)-float(kx)
            by=fay(i)-float(ky)
            bz=faz(i)-float(kz)
            ar=bx*bx+by*by+bz*bz
            if(ar.lt.armin) then
               kk=k
               mx=kx
               my=ky
               mz=kz
               armin=ar	
            endif
 1009    continue
         kkkk=kk
         if(kk.EQ.0) then
            write(20,*)' 	ERROR in the CHAIN global move' 
            GO TO  113          !do not move
         endif
         nx(i)=mx
         ny(i)=my
         nz(i)=mz
         px=mx
         py=my
         pz=mz
 101  continue
***************** new ica(i):
      ic=0
      do i=1,lenf1
         j=i+1
         wx=nx(j)-nx(i)
         wy=ny(j)-ny(i)
         wz=nz(j)-nz(i)
         nn(i)=vector(wx,wy,wz)
         if(nn(i).ne.oo(i)) ic=ic+1
      enddo
      if(ic.eq.0) go to 113     !do not move
c^^^^^^^^^^^^^^^^^^ new conformation obtained ^^^^^^^^^^^^^^^^^^^^^^^^^

c     change conformation to new path -------------->
      do i=1,lenf
         x(i)=nx(i)
         y(i)=ny(i)
         z(i)=nz(i)
         ica(i)=nn(i)
      enddo
      ica(0)=ica(2)
      ica(Lch)=ica(Lch-2)
      
      if(look(2,lenf1))then
c     calculate E_new--------------->
         do pp=2,lenf1
            nop(pp)=nopp(pp)
            nom(pp)=nomm(pp)
            noa(pp)=noaa(pp)
         enddo
         Enew=EHB(2,lenf1,1)+ESHORT(2,lenf1,1)
         do kkk=2,lenf1
            afsn(kkk)=afs(kkk)  !keep memory of ifs
         enddo

c     return back the conformation and calculate E_old --------->
         do i=1,lenf
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
         Eold=EHB(2,lenf1,-1)+ESHORT(2,lenf1,-1)
      
c     calculate eprofn while dord was calculated when call EHB(m,m3,1)--->
         eprofn=0.0
         do pp=2,lenf1
            is=seq(pp)
            ia=noa(pp)
            ip=nop(pp)
            im=nom(pp)
            eprofn=eprofn+envir(ia,im,ip,is,3)
         enddo
      
c     Metropolis ------------------>
         dE=Enew-Eold+dord+en1*(eprofn-eprofo)
c         E=energ
c         Et=E+de

         call metro(dE,id) !id=1, accepted; id=3, rejected

         if(id.eq.1)then        !accepted
            bNa(15)=bNa(15)+1
            bNNa(itemp)=bNNa(itemp)+1
            icnto=icnt
            sumcto=sumct 	
            do pp=2,lenf1
               nopp(pp)=nop(pp)
               nomm(pp)=nom(pp)
               noaa(pp)=noa(pp)
            enddo
            eprofo=eprofn

c     change the conformation to the new position--------->
            do kkk=2,lenf1
               afs(kkk)=afsn(kkk) !change ifs to the new one
            enddo
            codevsum=codev
            didevsum=didev	
c            energ=energ+de	
            do i=1,lenf
               x(i)=nx(i)
               y(i)=ny(i)
               z(i)=nz(i)
               ica(i)=nn(i)
            enddo
            ica(0)=ica(2)
            ica(Lch)=ica(Lch-2)
            call get_center     !update (amx,amy,amz)
c            call get_axis
         endif
      else                      !for look(i,m) loop
c     change back to old conformation ----------------------->
         do i=1,lenf
            x(i)=ox(i)
            y(i)=oy(i)
            z(i)=oz(i)
            ica(i)=oo(i)
         enddo
         ica(0)=ica(2)
         ica(Lch)=ica(Lch-2)
      endif
 113  continue
      bNt(15)=bNt(15)+1
      bNNt(itemp)=bNNt(itemp)+1

c ^^^^^^^^^^^^^^^^^ move9 finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

ccccccccccccccc Try to swap i1-replica and i2-replica ccccccccccccccccc
      subroutine swap(i1,i2)
      implicit integer (i-z)
      parameter(nrep=100)
      parameter(ndim=1999)	!maximum length of chain-length
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/nswap/bNSa(100),bNSt(100)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/chain0/ras(ndim),nfl
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      
      common/trackn/n_tem(100)

      if(metro_swap(i1,i2).eq.1)then
         bNSa(i1)=bNSa(i1)+1    !number of accepted swaps
         do i=1,Lch
*     x:
            tempx=xrep(i,i1)
            tempy=yrep(i,i1)
            tempz=zrep(i,i1)
            xrep(i,i1)=xrep(i,i2)
            yrep(i,i1)=yrep(i,i2)
            zrep(i,i1)=zrep(i,i2)
            xrep(i,i2)=tempx
            yrep(i,i2)=tempy
            zrep(i,i2)=tempz
         enddo
***   swap E ----------------------------->
         attt=E_rep(i1)         !exchange E_rep
         E_rep(i1)=E_rep(i2)
         E_rep(i2)=attt
***   swap n_tem ------------------------->
         temp=n_tem(i1)
         n_tem(i1)=n_tem(i2)
         n_tem(i2)=temp
***   swap moveable points---------------->
      endif
      bNSt(i1)=bNSt(i1)+1       !number of total swaps.

      return
      end

ccccccccccccccc swap: id=1, accepted; id=3, rejected ccccccccccccccccc
      integer function metro_swap(i,j)
      parameter(nrep=100)
      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/ranzy/nozy

c     w_12=wi(j)wj(i)/wi(i)wj(j) ------------->
      aaa=(1/aT_rep(i)-1/aT_rep(j))*(E_rep(i)-E_rep(j))
      aweight=exp(aaa)

      if(aranzy(nozy).le.aweight)then
         metro_swap=1		!swap accepted
      else
         metro_swap=3		!swap rejected
      endif

      return
      end

ccccccccccccccc Try to swap i1-replica and i2-replica ccccccccccccccccc
c     swap the replicas with combined chains
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine swap_RS(i1,i2)
      implicit integer (i-z)
      parameter(nrep=100)
      parameter(ndim=1999)	!maximum length of chain-length
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/nswap/bNSa(100),bNSt(100)

      common/sw1/aT_rep(nrep),E_rep(nrep)
      common/sw2/xrep(ndim,nrep),yrep(ndim,nrep),zrep(ndim,nrep)
      common/sw3/icarep(ndim,nrep)
      common/sw4/exrep(ndim,nrep),eyrep(ndim,nrep),ezrep(ndim,nrep)
      common/sw5/egxrep(ndim,nrep),egyrep(ndim,nrep),egzrep(ndim,nrep)
      common/sw7/ecxrep(ndim,nrep),ecyrep(ndim,nrep),eczrep(ndim,nrep)
      common/sw8/ebxrep(ndim,nrep),ebyrep(ndim,nrep),ebzrep(ndim,nrep)
      common/sw9/etxrep(ndim,nrep),etyrep(ndim,nrep),etzrep(ndim,nrep)
      common/chain0/ras(ndim),nfl
      common/echain0/nfr,nfr_i(ndim),nfr_f(ndim),ifr,n_fra
      common/chainm/mv(ndim)

      if(metro_swap(i1,i2).eq.1)then
         bNSa(i1)=bNSa(i1)+1    !number of accepted swaps
         do i=1,Lch
            t_ica=icarep(i,i1)             !ica(i)
            icarep(i,i1)=icarep(i,i2)
            icarep(i,i2)=t_ica
            if(mv(i).gt.0)then
               t_x=xrep(i,i1)              !x(i)
               t_y=yrep(i,i1)
               t_z=zrep(i,i1)
               xrep(i,i1)=xrep(i,i2)
               yrep(i,i1)=yrep(i,i2)
               zrep(i,i1)=zrep(i,i2)
               xrep(i,i2)=t_x
               yrep(i,i2)=t_y
               zrep(i,i2)=t_z
            else
               e_x=exrep(i,i1)           !ex(i)
               e_y=eyrep(i,i1)
               e_z=ezrep(i,i1)
               exrep(i,i1)=exrep(i,i2)
               eyrep(i,i1)=eyrep(i,i2)
               ezrep(i,i1)=ezrep(i,i2)
               exrep(i,i2)=e_x
               eyrep(i,i2)=e_y
               ezrep(i,i2)=e_z
               e_x=egxrep(i,i1)          !egx(i)
               e_y=egyrep(i,i1)
               e_z=egzrep(i,i1)
               egxrep(i,i1)=egxrep(i,i2)
               egyrep(i,i1)=egyrep(i,i2)
               egzrep(i,i1)=egzrep(i,i2)
               egxrep(i,i2)=e_x
               egyrep(i,i2)=e_y
               egzrep(i,i2)=e_z
               e_x=ecxrep(i,i1)          !ecx(i)
               e_y=ecyrep(i,i1)
               e_z=eczrep(i,i1)
               ecxrep(i,i1)=ecxrep(i,i2)
               ecyrep(i,i1)=ecyrep(i,i2)
               eczrep(i,i1)=eczrep(i,i2)
               ecxrep(i,i2)=e_x
               ecyrep(i,i2)=e_y
               eczrep(i,i2)=e_z
               e_x=ebxrep(i,i1)          !ebx(i)
               e_y=ebyrep(i,i1)
               e_z=ebzrep(i,i1)
               ebxrep(i,i1)=ebxrep(i,i2)
               ebyrep(i,i1)=ebyrep(i,i2)
               ebzrep(i,i1)=ebzrep(i,i2)
               ebxrep(i,i2)=e_x
               ebyrep(i,i2)=e_y
               ebzrep(i,i2)=e_z
               e_x=etxrep(i,i1)          !etx(i)
               e_y=etyrep(i,i1)
               e_z=etzrep(i,i1)
               etxrep(i,i1)=etxrep(i,i2)
               etyrep(i,i1)=etyrep(i,i2)
               etzrep(i,i1)=etzrep(i,i2)
               etxrep(i,i2)=e_x
               etyrep(i,i2)=e_y
               etzrep(i,i2)=e_z
            endif
         enddo
***   swap E ----------------------------->
         attt=E_rep(i1)         !exchange E_rep
         E_rep(i1)=E_rep(i2)
         E_rep(i2)=attt
***   swap moveable points---------------->
      endif
      bNSt(i1)=bNSt(i1)+1       !number of total swaps.

c^^^^^^^^^^^^^^^^^ swap_RS finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test whether all the neighboring backbones are good neighbors
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine test_neighbor
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nrep=100)
      parameter(nvec=416)
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/lengths/Lch,Lch1,Lch2
      common/logica/goodc
      logical goodc(nvec,nvec)
      common/vectors/vx(nvec),vy(nvec),vz(nvec),vector(-5:5,-5:5,-5:5)

      do j=1,Lch
         i=j-1
         if(j.ne.1.and.j.ne.Lch) then
            ii=ica(i)
            jj=ica(j)
            if(.not.goodc(ii,jj)) then
           write(20,8112)i,j,vx(ii),vy(ii),vz(ii),vx(jj),vy(jj),vz(jj)
 8112          format(5x,'warning2 -wrong input chain - vectors ',8i4)
               stop
            endif
         endif
      enddo 

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Test for overlaps of Ca:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine test_overlap
      implicit integer(i-z)
      common/lengths/Lch,Lch1,Lch2
      common/looks/exc,exc1,exc2

      do i=1,Lch
         do j=i+3,Lch
            dis2=(aax(i)-aax(j))**2+(aay(i)-aay(j))**2+
     $           (aaz(i)-aaz(j))**2
            if(dis2.lt.exc)then
               write(*,*)i,j,dis2,'  Ca-Ca overlap'
            endif
         enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate energy from the beginning
c     parameters will be modified after running this subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function energy_tot()
      implicit integer(i-z)
      parameter(ndim=1999)
      parameter(nvec=416)
      common/lengths/Lch,Lch1,Lch2
      common/envir1/nop(ndim),nom(ndim),noa(ndim),nhbn(ndim)
      common/lim/colim,dilim,coold,conew,diold,dinew,didev,codev    
      COMMON/short2/ codevsum, didevsum, csr(ndim,2)
      common/order/acorder,en2,sumct,sumcto,icnt,icnto,dord,en3
      common/backup1/NOPP(ndim),NOMM(ndim),NOAA(ndim),NHBNN(ndim)
      common/backup2/eprofo,eprofn,energ
      common/seqe/seq(ndim),sec(ndim)
      common/ehbc/envir(0:15,0:15,0:15,0:19,4),en1
      common/one/acrit,contt,eonekd(0:19),eoinp(0:19,0:100),es2,es1
      COMMON/ENERGY/EH5,ES3,ES3a,ES3b,EH1,EH1a,EH4,EHBIJ(ndim,ndim)
      common/pair1/eh2,eh1b,eh1c
      common/shortcom/eh3,es4,es5,es6,es7,es7a,es7b,es7c
      COMMON/RCN1/ER1,arca1(ndim,ndim,50),n_resa1(ndim,ndim)
      COMMON/RES/ER3,er5,er6,er7,Mcom(ndim),Kcom(ndim,100)
      common/otherenergy/E_cord,E_cnum

      common/ehbenergy/EHB1,EHB1a,EHB1b,EHB1c,EHB2,EHB3,EHB4,EHB5,EHB6
      common/ehbenergy1/EHB5a,EHB5b

      ICNTO=0
      SUMCTO=0
      DO i=1,Lch
         NOP(i)=0
         NOA(i)=0
         NOM(i)=0
      enddo

      energy_tot=EHB(1,Lch,1)+ESHORT(1,Lch,10)

      eprofo=0.0
      DO k=1,Lch
         is=seq(k)
         ia=NOA(k)
         ip=NOP(k)
         im=NOM(k)
         eprofo=eprofo+envir(ia,im,ip,is,3)
      enddo

      E_cord=abs(ICNT/(0.00001+float(SUMCT))-acorder)
      E_cnum=abs(float(SUMCT)-contt)

      energy_tot=energy_tot+en1*eprofo+en2*E_cord+en3*E_cnum

c     following is for next movement ---------->
      icnto=icnt                !backup of contact-order
      sumcto=sumct              !backup of contact-order
      do i=1,Lch
         nopp(i)=nop(i)         !backup of number of contact-pair
         noaa(i)=noa(i)         !backup of number of contact-pair
         nomm(i)=nom(i)         !backup of number of contact-pair
      enddo
      codevsum=conew            !backup of panelity for total deviation of comb
      didevsum=dinew            !backup of panelity for total deviation of dist

c^^^^^^^^^^^^^^^^^^^^ Energy_tot finished ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

cccccccccccccccc Metropolis: id=1, accepted; id=3, rejected ccccccccccccccccc
      subroutine metro(dE,id)
      common/ranzy/nozy

      id=1
      if(dE.gt.0)then
         if(aranzy(nozy).gt.weight12(dE))then
            id=3
         endif
      endif
      return
      end

ccccccccccccccccweight factor cccccccccccccccccccccccccccccccccccccccccccc
      function weight12(dE)
      parameter(nrep=100)
      common/temperature/itemp,atemp
      common/aTs/aTs1,aTs2,aTs_rep(nrep),aTs
      common/aTTs/aTTs1,aTTs2,aTTs_rep(nrep),aTTs
      common/ichos/ichos

      if(ichos.eq.1)then
         weight12=exp(-dE/atemp)
         ichos=2
      else
         weight12=exp(-square2(dE)/aTs)
         ichos=1
      endif
c     write(*,*)ichos,aTs,aTTs,dE,weight12

      return
      end

ccccccccccccccc modified energy landscape ccccccccccccccccccccccccccc
      function square2(x)
      square2=x*x
      return
      end

ccccccccccccccc flat function ccccccccccccccccccccccccccc
      function arcsh(x)
      arcsh=log(x+sqrt(1+x*x))
      return
      end

ccccccccccccccccc calculate RMSD and DRMS cccccccccccccccccccccccccccc
c x_data(i)=r_d(1,i); y_data(i)=r_d(2,i); z_data(i)=r_d(3,i);
c x_model(i)=r_m(1,i); y_model(i)=r_m(2,i); z_model(i)=r_m(3,i);
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function zyrmsd(arms)
      implicit integer (i-z)
      parameter(ndim=1999)	!maximum length of chain-length
      common/chain1/ica(0:ndim),x(ndim),y(ndim),z(ndim)
      common/rmsdrange/nca1,nca2
      common/CA/dx(ndim),dy(ndim),dz(ndim)
      double precision r_m(3,ndim),r_d(3,ndim),w(ndim)
      double precision u(3,3),t(3),rms !armsd is real
      data w /ndim*1.0/

      k=0
      do i=nca1,nca2
         k=k+1
         r_m(1,k)=x(i+1)*0.87   !model conformation
         r_m(2,k)=y(i+1)*0.87
         r_m(3,k)=z(i+1)*0.87
         r_d(1,k)=dx(i)         !native conformation
         r_d(2,k)=dy(i)
         r_d(3,k)=dz(i)
      enddo

      nn=k                      !number of data points
      call u3b(w,r_m,r_d,nn,0,rms,u,t,ier)
      arms=dsqrt(rms/nn)	!RMSD is real, rms is double precision
      zyrmsd=arms

      return
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
c.....END U3B...........................................................
c----------------------------------------------------------
c                       THE END
c----------------------------------------------------------
c 338 "rms.for"
      end

**************************************************************
*     Calculate eigenvalue and eignevectors
*     A{3,3}-matrix E{3}-eigenvalue T{3,3}-eigenvector
*     E(1)=<x^2>; E(1)=<y^2>; E(1)=<z^2> in the rotated system
**************************************************************
      subroutine eigenvalue
      common/eigen/A(3,3),E(3),T(3,3)

c     #### Eigenvalues ############
      pi=3.1415926
      p1=-1
      p2=A(1,1)+A(2,2)+A(3,3) 
      p3=-A(1,1)*A(2,2)-A(1,1)*A(3,3)-A(2,2)*A(3,3)+
     &     A(2,3)*A(2,3)+A(1,3)*A(1,3)+A(1,2)*A(1,2)
      p4=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+
     &     A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)-
     &     A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)
      p5=(-(1.0/3)*(p2/p1)**2+p3/p1)/3 
      ap5=sqrt(-p5)
      p6=((2.0/27)*(p2/p1)**3-(1.0/3)*(p2*p3/p1**2)+p4/p1)/2
      p7=acos(-p6/sqrt(-p5**3))
      p8=2*ap5*cos(p7/3.0)
      p9=-2*ap5*cos((p7+pi)/3.0)
      p10=-2*ap5*cos((p7-pi)/3.0) 
      p11=p2/(3*p1)
      E(1)=p8-p11               !eigenvalue
      E(2)=p9-p11               !eigenvalue
      E(3)=p10-p11              !eigenvalue
      
c     ##### normalized eigenvectors #########
      do i=1,3
         fnorm1=A(2,1)*A(1,2)-(A(1,1)-E(i))*(A(2,2)-E(i))
         x=((A(2,2)-E(i))*A(1,3)-A(1,2)*A(2,3))/fnorm1
         y=((A(1,1)-E(i))*A(2,3)-A(2,1)*A(1,3))/fnorm1
         T(i,3)=1/sqrt(x*x+y*y+1)
         T(i,1)=x*T(i,3)
         T(i,2)=y*T(i,3)
c     write(*,*)i,E(i),T(i,1),T(i,2),T(i,3)
      enddo

      return
      end

****************************************************
c     common/ranzy/nozy                            *
*     nozy=-1 (less than 0)                        *
*     rrr=ranzy(nozy)                              *
****************************************************
      function aranzy(nozy)
c   random numbers uniformly distributed between 0 and 1.
c   (code by Y. Zhang, ITP, Acdemia Sincica, 1999.6
c   can be replaced by any other suitable generator for random numbers)
      common/asdfa/ix1,ix2,ix3,rm1,rm2,r(99)
      data ia1,ic1,m1/1279,351762,1664557/
      data ia2,ic2,m2/2011,221592,1048583/
      data ia3,ic3,m3/15551,6150,29101/
      if(nozy.ge.0) go to 2
      ix1=mod(-nozy,m1)
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ix1,m2)
      ix1=mod(ia1*ix1+ic1,m1)
      ix3=mod(ix1,m3)
      rm1=1./float(m1)
      rm2=1./float(m2)
      do 1 j=1,99
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
 1    r(j)=(float(ix1)+float(ix2)*rm2)*rm1
 2    ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(99*ix3)/m3
      aranzy=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      nozy=ix1
      return
      end

