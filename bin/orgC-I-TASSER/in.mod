$HOUR$
$RANDOM$	$NCYCLE$	200	40      $NRUN$
0.48	0.1	0.1	0.05	0.05	0.04	0.04	0.08	0.06
80	1.7	23	18	18	6	2.5
1.0	0.2     5	3
$SWITCH$	$ITHR0$	3
 4.1025 48      0.45    0.1     0       5
 4.6538 50      0.25    0.4     0.35    20
 8.0122 4.1010  0.01    2.011
 1.8252 9.925
 0.6314 3.0197  1.3309  1.013   2.1750
 1.5078 1.7847  0.6052
 0.6143 0.6014  0.1
 1    12    1    $IFA$      chuan,nana,m_latt,ifa
2000 6
7 0.5
0	1.6	0 	0.5
$AW1$  $AW2$  $AW4$   $CR20$
$ER6$   $DWELL$  $FW$
      
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
      read(19,*) aw1,aw2,aw4,Cr20 !comb,combCA,par/pair1
      read(19,*) er6,dwell      !weight for distL
      
