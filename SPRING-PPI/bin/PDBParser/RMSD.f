*************************************************************************
*     This program is to calculate the root-mean-square deviation (RMSD)
*
*     To calculate RMSD of 'pdb1' and 'pdb2':
*       >RMSD pdb1 pdb2
*     
*     To output transferred 'pdb1_tr'
*       >RMSD pdb1 pdb2 -o pdb1_tr
*     
*     Note:
*       1. Input structure much be in PDB format
*       2. Only C-alpha RMSD is reported
*     
*     Please address your question and comments to zhng@umich.edu
*************************************************************************
      
      program calRMSD
      PARAMETER(nmax=3000)
      
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0,d0_fix
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)

      character*500 fnam,pdb(100),outname
      character*3 aa(-1:20),seqA(nmax),seqB(nmax)
      character*500 s,du
      character seq1A(nmax),seq1B(nmax),ali(nmax)
      character sequenceA(nmax),sequenceB(nmax),sequenceM(nmax)

      dimension L_ini(100),iq(nmax)
      common/scores/score,score_maxsub,score_fix,score10
      common/GDT/n_GDT05,n_GDT1,n_GDT2,n_GDT4,n_GDT8
      double precision score,score_max,score_fix,score_fix_max
      double precision score_maxsub,score10
      dimension xa(nmax),ya(nmax),za(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

*****instructions ----------------->
      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
         write(*,*)
         write(*,*)'To calculate RMSD of ''pdb1'' and ''pdb2'':'
         write(*,*)'   >RMSD pdb1 pdb2'
         write(*,*)
         write(*,*)'To output transferred ''pdb1_tr'''
         write(*,*)'   >RMSD pdb1 pdb2 -o pdb1_tr'
         write(*,*)
         write(*,*)'Note:'
         write(*,*)'   1. Input structure much be in PDB format'
         write(*,*)'   2. Only C-alpha RMSD is reported'
         write(*,*)
         goto 9999
      endif
      
******* options ----------->
      narg=iargc()
      if(narg<2)then
         write(*,*)'Error, please enter two PDB files'
         goto 9999
      endif
      i=0
      j=0
      m_out=0
 115  continue
      i=i+1
      call getarg(i,fnam)
      if(fnam.eq.'-o')then
         m_out=1
         i=i+1
         call getarg(i,outname)
      else
         j=j+1
         pdb(j)=fnam
      endif
      if(i.lt.narg)goto 115
      
ccccccccc read data from first CA file:
      open(unit=10,file=pdb(1),status='old')
      i=0
 101  read(10,104,end=102) s
      if(s(1:3).eq.'TER') goto 102
      if(s(1:4).eq.'ATOM')then
         if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &        eq.'  CA')then
         if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
            i=i+1
            read(s,103)du,seqA(i),du,nresA(i),du,xa(i),ya(i),za(i)
         endif
         endif
      endif
      goto 101
 102  continue
 103  format(A17,A3,A2,i4,A4,3F8.3)
 104  format(A100)
      close(10)
      nseqA=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
ccccccccc read data from second CA file:
      open(unit=10,file=pdb(2),status='old')
      i=0
 201  read(10,204,end=202) s
      if(s(1:3).eq.'TER') goto 202
      if(s(1:4).eq.'ATOM')then
         if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.s(13:16).
     &        eq.'  CA')then
         if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
            i=i+1
            read(s,203)du,seqB(i),du,nresB(i),du,xb(i),yb(i),zb(i)
         endif
         endif
      endif
      goto 201
 202  continue
 203  format(A17,A3,A2,i4,A4,3F8.3)
 204  format(A100)
      close(10)
      nseqB=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
******************************************************************
*     pickup the common residues:
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
      if(n_ali.lt.1)then
        write(*,*)'There is no common residues in the input structures'
        goto 9999
      endif

***   RMSD calculation ---------------------------------->
      do m=1,n_ali
         r_1(1,m)=xa(iA(m))
         r_1(2,m)=ya(iA(m))
         r_1(3,m)=za(iA(m))
         r_2(1,m)=xb(iB(m))
         r_2(2,m)=yb(iB(m))
         r_2(3,m)=zb(iB(m))
      enddo
      call u3b(w,r_1,r_2,n_ali,1,rms,u,t,ier) !u rotate r_1 to r_2
      rmsd=dsqrt(rms/n_ali)
      call zydrmsd(r_1,r_2,n_ali,drms)
      drmsd=drms
      
******************************************************************
*     Output
******************************************************************
***   output RMSD ---------------------------->
      write(*,*)
      write(*,501)pdb(1),nseqA
 501  format('Structure1: ',A10,'  Length= ',I4)
      write(*,502)pdb(2),nseqB
 502  format('Structure2: ',A10,'  Length= ',I4)
      write(*,503)n_ali
 503  format('Number of residues in common= ',I4)
      write(*,513)rmsd
 513  format('RMSD of the common residues= ',F8.3)
      write(*,514)drmsd
 514  format('Distant RMSD of the common residues= ',F8.3)
      write(*,*)
      write(*,*) '-------- rotation matrix to rotate Chain-1 to ',
     &     'Chain-2 ------'
      write(*,*)'i          t(i)         u(i,1)         u(i,2) ',
     &     '        u(i,3)'
      do i=1,3
         write(*,304)i,t(i),u(i,1),u(i,2),u(i,3)
      enddo
      write(*,*)
304   format(I2,f18.10,f15.10,f15.10,f15.10)
***   output rotated structure of pdb1:
      if(m_out.eq.1)then
         OPEN(unit=7,file=outname,status='unknown') !pdb1_tr
         do j=1,nseqA
            xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
            yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
            zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
            write(7,1237)nresA(j),seqA(j),nresA(j),
     &           xt(j),yt(j),zt(j)
         enddo
         write(7,1238)
         close(7)
      endif
 1237 format('ATOM  ',i5,'  CA  ',A3,I6,4X,3F8.3)
 1238 format('TER')
***   

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 9999 END

ccccccccccc Calculate dRMSD between 1 and CA cccccccccccccccccccccccccccccccccc
      subroutine zydrmsd(r_1,r_2,nn,drms)
      double precision r_1(3,1),r_2(3,1),drms,d,dp
      drms=0
      do i=1,nn
         do j=i+1,nn
            d=sqrt((r_1(1,i)-r_1(1,j))**2+(r_1(2,i)-r_1(2,j))**2
     $           +(r_1(3,i)-r_1(3,j))**2)
            dp=sqrt((r_2(1,i)-r_2(1,j))**2+(r_2(2,i)-r_2(2,j))**2
     $           +(r_2(3,i)-r_2(3,j))**2)
            drms=drms+(d-dp)**2
         enddo
      enddo
      drms=sqrt(drms*2/(nn*(nn-1)))
      
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
      end
