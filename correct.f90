!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!......3-D correct function......
!
subroutine correct3d(b,bd,cpt,cjp,dcjp,dwjp,nep,iInter,inf,ninf,maxconn)
  implicit none

  integer nep,iInter
  integer maxconn,ninf,inf(maxconn)
  real(8) b(*),bd(3,*),cpt(3)
  real(8) cjp(3,nep),dcjp(3,nep),dwjp(nep)

  if (iInter .eq. 1) then
     call correct3dl(b,bd,cpt,cjp,dcjp,dwjp,nep,inf,ninf,maxconn)

! Lucy commented it out, don't know what iInter=11 is.
!  elseif(iInter .eq. 11)  then
!     call correct3dtl(b,bd,cpt,cjp,dcjp,dwjp,nep)

  else
     print *, 'wrong iInter'
     stop
  endif
  return
end subroutine correct3d

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!......2-D correct function......
!
subroutine correct2d(b,bd,cpt,cjp,dcjp,dwjp,nep,iInter,inf,ninf,maxconn)
  implicit none

  integer nep,iInter
  integer maxconn,ninf,inf(maxconn)
  real(8) b(*),bd(3,*),cpt(3)
  real(8) cjp(3,nep),dcjp(3,nep),dwjp(nep)

  if (iInter .eq. 1) then
     call correct2dl(b,bd,cpt,cjp,dcjp,dwjp,nep,inf,ninf,maxconn)
  else
     print *, 'wrong iInter'
     stop
  endif
  return
end subroutine correct2d

!********************************************************************

subroutine correct3dl(b,bd,cpt,cjp,anode,dwjp,nep,inf,ninf,maxconn)
!**************************************************************
      implicit none

      integer nep
      integer ninf, maxconn, iptc
      real(8) cjp(3,nep),anode(3,nep),dwjp(nep)
      real(8) gm(4,4),gminv(4,4)
      real(8) gmdx(4,4),gmdy(4,4),gmdz(4,4)
      real(8) adx(4,4),ady(4,4),adz(4,4)
      real(8) b(4),bd(3,4),cpt(3)
      integer  inf(maxconn)

      real(8) zero,flag,det
      real(8) am000,am100,am010,am001
      real(8)       am200,am110,am101
      real(8)             am020,am011
      real(8)                   am002

      real(8) am000dx,am100dx,am010dx,am001dx
      real(8)         am200dx,am110dx,am101dx
      real(8)                 am020dx,am011dx
      real(8)                         am002dx

      real(8) am000dy,am100dy,am010dy,am001dy
      real(8)         am200dy,am110dy,am101dy
      real(8)                 am020dy,am011dy
      real(8)                         am002dy

      real(8) am000dz,am100dz,am010dz,am001dz
      real(8)         am200dz,am110dz,am101dz
      real(8)                 am020dz,am011dz
      real(8)                         am002dz

      real(8) xi,yi,zi,dsj,dx,dy,dz
      real(8) xp,yp,zp,a1,a2,a3
      real(8) xx,yy,zz
      real(8) r100,r010,r001,r200,r110,r101
      real(8) r020,r011,r002
      real(8) aw,awdx,awdy,awdz
      real(8) awdxy,awdyz,awdzx
      real(8) awdxx,awdyy,awdzz

      integer i,j,k,ipt,nn,nd
!
!.....set the initial value for moment and its derivatives:
!
      zero  = 0.00
      nn    = 4
      nd    = 4

      am000 = 0.00
      am100 = 0.00
      am010 = 0.00
      am001 = 0.00

      am200 = 0.00
      am110 = 0.00
      am101 = 0.00

      am020 = 0.00
      am011 = 0.00

      am002 = 0.00

      am000dx = 0.00
      am100dx = 0.00
      am010dx = 0.00
      am001dx = 0.00

      am200dx = 0.00
      am110dx = 0.00
      am101dx = 0.00

      am020dx = 0.00
      am011dx = 0.00

      am002dx = 0.00

      am000dy = 0.00
      am100dy = 0.00
      am010dy = 0.00
      am001dy = 0.00

      am200dy = 0.00
      am110dy = 0.00
      am101dy = 0.00

      am020dy = 0.00
      am011dy = 0.00

      am002dy = 0.00

      am000dz = 0.00
      am100dz = 0.00
      am010dz = 0.00
      am001dz = 0.00

      am200dz = 0.00
      am110dz = 0.00
      am101dz = 0.00

      am020dz = 0.00
      am011dz = 0.00

      am002dz = 0.00
!
!.....set the initial value for all array: 
!
      do i = 1, nn
        do  j = 1, nn
          gm(i,j)    = 0.00
          gminv(i,j) = 0.00
          gmdx(i,j)  = 0.00
          gmdy(i,j)  = 0.00
          gmdz(i,j)  = 0.00
        end do
        b(i)    = 0.00
        bd(1,i) = 0.00
        bd(2,i) = 0.00
        bd(3,i) = 0.00
      end do  
!
!.....input the value cpt 
!
      xp = cpt(1)
      yp = cpt(2)
      zp = cpt(3)
!
!.....main loop: calculate moment by Trapezodial rule
!
      do 30 iptc = 1, ninf
        ipt = inf(iptc)
!
!........define intermediate variable
!
!........a1, a2 and a3
!
        a1 =  anode(1,ipt)
        a2 =  anode(2,ipt)
        a3 =  anode(3,ipt)

        dsj = dwjp(ipt)

        xi  = cjp(1,ipt)
        yi  = cjp(2,ipt)
        zi  = cjp(3,ipt)
        dx  = -1.0/a1
        dy  = -1.0/a2
        dz  = -1.0/a3

        r100 = (xi - xp)/a1
        r010 = (yi - yp)/a2
        r001 = (zi - zp)/a3

        r200 = r100*r100
        r110 = r100*r010
        r101 = r100*r001

        r020 = r010*r010
        r011 = r010*r001

        r002 = r001*r001
  
        xx = abs(r100)
        yy = abs(r010)
        zz = abs(r001)

        if((xx.ge.2.0) .or. (yy.ge.2.0) .or. (zz .ge. 2.0)) go to 30

        call window3d(aw,awdx,awdy,awdz,  &
            awdxx,awdyy,awdzz,            &
            awdxy,awdyz,awdzx,            &
            xp,yp,zp,xi,yi,zi,a1,a2,a3)
  
        aw = aw*dsj
!  
!-----------------------------------
!  
        am000 = am000  +  aw
        am100 = am100  +  r100*aw
        am010 = am010  +  r010*aw
        am001 = am001  +  r001*aw
 
        am200 = am200  +  r200*aw
        am110 = am110  +  r110*aw
        am101 = am101  +  r101*aw
  
        am020 = am020  +  r020*aw
        am011 = am011  +  r011*aw
  
        am002 = am002  +  r002*aw
!  
!------------------------------------
!  
        awdx = awdx*dsj
        awdy = awdy*dsj
        awdz = awdz*dsj
!  
!  (1)
!  
        am000dx = am000dx + awdx
        am100dx = am100dx + r100*awdx + dx*aw
        am010dx = am010dx + r010*awdx
        am001dx = am001dx + r001*awdx
  
        am200dx = am200dx + r200*awdx + 2.0*dx*r100*aw
        am110dx = am110dx + r110*awdx + dx*r010*aw
        am101dx = am101dx + r101*awdx + dx*r001*aw    
  
        am020dx = am020dx + r020*awdx
        am011dx = am011dx + r011*awdx
  
        am002dx = am002dx + r002*awdx
!  
!  (2)
!  
        am000dy = am000dy + awdy
        am100dy = am100dy + r100*awdy
        am010dy = am010dy + r010*awdy + dy *aw
        am001dy = am001dy + r001*awdy
  
        am200dy = am200dy + r200*awdy
        am110dy = am110dy + r110*awdy + dy*r100*aw  
        am101dy = am101dy + r101*awdy
  
        am020dy = am020dy + r020*awdy + 2.0*dy*r010*aw 
        am011dy = am011dy + r011*awdy + dy*r001*aw    
  
        am002dy = am002dy + r002*awdy
!  
!  (3)
!  
        am000dz = am000dz + awdz
        am100dz = am100dz + r100*awdz
        am010dz = am010dz + r010*awdz
        am001dz = am001dz + r001*awdz + dz*aw
  
        am200dz = am200dz + r200*awdz
        am110dz = am110dz + r110*awdz
        am101dz = am101dz + r101*awdz + dz*r100*aw
  
        am020dz = am020dz + r020*awdz
        am011dz = am011dz + r011*awdz + dz*r010*aw   
  
        am002dz = am002dz + r002*awdz + 2.0*dz*r001*aw
  
 30   continue
!  
!  
!.....end of the main loop
!  
!.....assemble the M matrix (gm(i,j))
!  
      gm(1,1) = am000
      gm(1,2) = am100
      gm(1,3) = am010
      gm(1,4) = am001
  
      gm(2,1) = am100
      gm(2,2) = am200
      gm(2,3) = am110
      gm(2,4) = am101
  
      gm(3,1) = am010
      gm(3,2) = am110
      gm(3,3) = am020
      gm(3,4) = am011
  
      gm(4,1) = am001
      gm(4,2) = am101
      gm(4,3) = am011
      gm(4,4) = am002
!  
!.....calculate the determinat det
!  
      call gjinv(gm,gminv,nn,nd,det,flag)
!  
!.....test convergence criteria
!  
      if(det .le. zero) then
        print *, 'det =', det
        print *, 'STOP! the determinat det <= 0 '
        stop
      else
      end if
!  
!.....assemble the gminv(i,j)
!  
!.....calculate the derivative of gm
!  
      gmdx(1,1) = am000dx
      gmdx(1,2) = am100dx
      gmdx(1,3) = am010dx
      gmdx(1,4) = am001dx
  
      gmdx(2,1) = am100dx
      gmdx(2,2) = am200dx
      gmdx(2,3) = am110dx
      gmdx(2,4) = am101dx
  
      gmdx(3,1) = am010dx
      gmdx(3,2) = am110dx
      gmdx(3,3) = am020dx
      gmdx(3,4) = am011dx
  
      gmdx(4,1) = am001dx
      gmdx(4,2) = am101dx
      gmdx(4,3) = am011dx
      gmdx(4,4) = am002dx
  
  
      gmdy(1,1) = am000dy
      gmdy(1,2) = am100dy
      gmdy(1,3) = am010dy
      gmdy(1,4) = am001dy
  
      gmdy(2,1) = am100dy
      gmdy(2,2) = am200dy
      gmdy(2,3) = am110dy
      gmdy(2,4) = am101dy
  
      gmdy(3,1) = am010dy
      gmdy(3,2) = am110dy
      gmdy(3,3) = am020dy
      gmdy(3,4) = am011dy
  
      gmdy(4,1) = am001dy
      gmdy(4,2) = am101dy
      gmdy(4,3) = am011dy
      gmdy(4,4) = am002dy
  
  
      gmdz(1,1) = am000dz
      gmdz(1,2) = am100dz
      gmdz(1,3) = am010dz
      gmdz(1,4) = am001dz
  
      gmdz(2,1) = am100dz
      gmdz(2,2) = am200dz
      gmdz(2,3) = am110dz
      gmdz(2,4) = am101dz
  
      gmdz(3,1) = am010dz
      gmdz(3,2) = am110dz
      gmdz(3,3) = am020dz
      gmdz(3,4) = am011dz
  
      gmdz(4,1) = am001dz
      gmdz(4,2) = am101dz
      gmdz(4,3) = am011dz
      gmdz(4,4) = am002dz
!  
!.....assign the value for b vector
!  
      do i = 1, nn
        b(i) = gminv(1,i)
      enddo
!  
!.....find the value for bd(2,3)
!  
!.....compute M^(-1)dM/dx and M^(-1)dM/dy
!  
      do i = 1,nn
        do j = 1,nn
          adx(i,j) = 0.00    
          ady(i,j) = 0.00    
          adz(i,j) = 0.00    
          do k = 1,nn
            adx(i,j) = adx(i,j) + gminv(i,k)*gmdx(k,j)
            ady(i,j) = ady(i,j) + gminv(i,k)*gmdy(k,j)
            adz(i,j) = adz(i,j) + gminv(i,k)*gmdz(k,j)
          enddo
        enddo
      enddo
  
  
      do i = 1,nn
        do j = 1,nn
          bd(1,i) = bd(1,i) - adx(i,j)*b(j)
          bd(2,i) = bd(2,i) - ady(i,j)*b(j)
          bd(3,i) = bd(3,i) - adz(i,j)*b(j)
        end do
      end do


! 999  continue
 
      return
      end
      

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine correct3dtl(b,bd,cpt,cjp,anode,dwjp,nep)
!
!     This subroutine is to calculate the b vector and its
!     derivatives,for the moving least square reproducing kernel 
!     interpolation shape function.
!
!     This code is only offering b vector and its 1st derivatives
!     for the 3-D case (R3dtl).
!
!     The interpolation function is choosen as
!     
!     p = (1 ,x ,y ,z, xy, yz, zx)
!
!     Date: June, 1998
!
!     -arguments:
!
!  i      nep: the numbers of particles, therefore, the numbers
!              of elements too.
!
!  i      cjp(3,nep):  the particle's global coordinations;
!                      such as x(nep):= cjp(1,nep)
!                              y(nep):= cjp(2,nep)
!                              z(nep):= cjp(3,nep)
!
!  i      anode(3,nep): the array that stores all the dilation
!                     vectors of the particles .
!             The dilation parameter at X direction: anode(1,nep)
!             The dilation parameter at Y direction: anode(2,nep)
!             The dilation parameter at Z direction: anode(3,nep)
!
!  i      dwjp(nep): the integration weight at each particle point
!
!
!  i      cpt(3): the point at where the shape function is
!                     evaluated.;
!                 x(pt):= cpt(1)
!                 y(pt):= cpt(2)
!                 z(pt):= cpt(3)
!
!    
!   l     a1 : dilation parameter; scalar;
!
!   l     a2 : dilation parameter; scalar;
!
!   l     a3 : dilation parameter; scalar;
!
!               we input the adjust constants ax, ay, and we construct
!               the dilation parameter in such a way: 
!               (for non-uniform mesh )
!
!                 a1 = anode(1,j)
!                 a2 = anode(2,j)
!                 a3 = anode(3,j)
!
!    l   gm(7,7) : the basic moment matrix;
!
!    l   gminv(7,7): the inverse of the gm, i.e. gm^{-1}
!
!    l   gmdx(7,7): the derivatives of gm: gm_x;
!
!    l   gmdy(7,7): the derivatives of gm: gm_y;
!
!    l   adx(7,7):  M^(-1) M_x
!
!    l   ady(7,7):  M^(-1) M_y
!
!    l   adz(7,7):  M^(-1) M_z
!
!
!  o     b(7)   :   the b vector, also the first row of gminv
!
!  o     bd(3,7):   the 1st derivatives of the b vector
!
!                   bd(1,4) :=  b_x
!                   bd(2,4) :=  b_y
!                   bd(3,4) :=  b_z
!
!     The subroutine is only designed to generate 2-D
!     moving least square shape function based on linear polynomial,
!
!    l   p(7) := (1,x,y,z,xy,yz,zx);
!
!     The mathematical formulation is as follows
!
!     b   = (1/det)* [ a11, - a12, a13, -a14, a15, -a16, a17 ];
!
!     b_x = - M^{-1} M_x b;
!
!     b_y = - M^{-1} M_y b;
!
!     b_z = - M^{-1} M_z b;
!
!
!
!    There two different subroutines could be called by this subroutine:
!
!    cubic3d.f ----- the cubic spline window function;
!  
!    quintw3.f ----- the fifth order spline window function.
!
!**************************************************************
!
      implicit none
!
      integer nep
      real(8) cjp(3,nep),anode(3,nep),dwjp(nep)
      real(8) gm(7,7),gminv(7,7)
      real(8) gmdx(7,7),gmdy(7,7),gmdz(7,7)
      real(8) adx(7,7),ady(7,7),adz(7,7)
      real(8) b(*),bd(3,*),cpt(3)
!
      real(8) zero,flag,det
      real(8) am000,am100,am010,am001,am110,am011,am101
      real(8)       am200,            am210,am111,am201
      real(8)             am020,      am120,am021
      real(8)                   am002,      am012,am102
      real(8)                         am220,am121,am211
      real(8)                               am022,am112
      real(8)                                     am202
!
      real(8) am000dx,am100dx,am010dx,am001dx,am110dx,am011dx,am101dx
      real(8)         am200dx,                am210dx,am111dx,am201dx
      real(8)                 am020dx,        am120dx,am021dx
      real(8)                         am002dx,        am012dx,am102dx
      real(8)                                 am220dx,am121dx,am211dx
      real(8)                                         am022dx,am112dx
      real(8)                                                 am202dx
!
      real(8) am000dy,am100dy,am010dy,am001dy,am110dy,am011dy,am101dy
      real(8)         am200dy,                am210dy,am111dy,am201dy
      real(8)                 am020dy,        am120dy,am021dy
      real(8)                         am002dy,        am012dy,am102dy
      real(8)                                 am220dy,am121dy,am211dy
      real(8)                                         am022dy,am112dy
      real(8)                                                 am202dy
!
      real(8) am000dz,am100dz,am010dz,am001dz,am110dz,am011dz,am101dz
      real(8)         am200dz,                am210dz,am111dz,am201dz
      real(8)                 am020dz,        am120dz,am021dz
      real(8)                         am002dz,        am012dz,am102dz
      real(8)                                 am220dz,am121dz,am211dz
      real(8)                                         am022dz,am112dz
      real(8)                                                 am202dz
!
      real(8) xi,yi,zi,dsj,dx,dy,dz
      real(8) xp,yp,zp,a1,a2,a3
      real(8) xx,yy,zz
      real(8) r100,r010,r001,r110,r011,r101
      real(8) r200,r210,r111,r201
      real(8) r020,r120,r021
      real(8) r002,r012,r102
      real(8) r220,r121,r211
      real(8) r022,r112,r202
      real(8) aw,awdx,awdy,awdz
      real(8) awdxy,awdyz,awdzx
      real(8) awdxx,awdyy,awdzz
!
      integer i,j,k,ipt,nn,nd
!
!.....set the initial value for moment and its derivatives:
!
      zero  = 0.00
      nn    = 7
      nd    = 7
!
      am000 = 0.00
      am100 = 0.00
      am010 = 0.00
      am001 = 0.00
      am110 = 0.00
      am011 = 0.00
      am101 = 0.00
!
      am200 = 0.00
      am210 = 0.00
      am111 = 0.00
      am201 = 0.00
!
      am020 = 0.00
      am120 = 0.00
      am021 = 0.00
!
      am002 = 0.00
      am012 = 0.00
      am102 = 0.00
!
      am220 = 0.00
      am121 = 0.00
      am211 = 0.00
!
      am022 = 0.00
      am112 = 0.00
!
      am202 = 0.00
!
      am000dx = 0.00
      am100dx = 0.00
      am010dx = 0.00
      am001dx = 0.00
      am110dx = 0.00
      am011dx = 0.00
      am101dx = 0.00
!
      am200dx = 0.00
      am210dx = 0.00
      am111dx = 0.00
      am201dx = 0.00
!
      am020dx = 0.00
      am120dx = 0.00
      am021dx = 0.00
!
      am002dx = 0.00
      am012dx = 0.00
      am102dx = 0.00
!
      am220dx = 0.00
      am121dx = 0.00
      am211dx = 0.00
!
      am022dx = 0.00
      am112dx = 0.00
!
      am202dx = 0.00

      am000dy = 0.00
      am100dy = 0.00
      am010dy = 0.00
      am001dy = 0.00
      am110dy = 0.00
      am011dy = 0.00
      am101dy = 0.00

      am200dy = 0.00
      am210dy = 0.00
      am111dy = 0.00
      am201dy = 0.00

      am020dy = 0.00
      am120dy = 0.00
      am021dy = 0.00

      am002dy = 0.00
      am012dy = 0.00
      am102dy = 0.00

      am220dy = 0.00
      am121dy = 0.00
      am211dy = 0.00

      am022dy = 0.00
      am112dy = 0.00

      am202dy = 0.00

      am000dz = 0.00
      am100dz = 0.00
      am010dz = 0.00
      am001dz = 0.00
      am110dz = 0.00
      am011dz = 0.00
      am101dz = 0.00

      am200dz = 0.00
      am210dz = 0.00
      am111dz = 0.00
      am201dz = 0.00

      am020dz = 0.00
      am120dz = 0.00
      am021dz = 0.00

      am002dz = 0.00
      am012dz = 0.00
      am102dz = 0.00

      am220dz = 0.00
      am121dz = 0.00
      am211dz = 0.00

      am022dz = 0.00
      am112dz = 0.00

      am202dz = 0.00
!
!
!.....set the initial value for all array: 
!
      do i = 1, nn
         do j = 1, nn
            gm(i,j)    = 0.00
            gminv(i,j) = 0.00
            gmdx(i,j)  = 0.00
            gmdy(i,j)  = 0.00
            gmdz(i,j)  = 0.00
         end do
         b(i)    = 0.00
         bd(1,i) = 0.00
         bd(2,i) = 0.00
         bd(3,i) = 0.00
      end do  
!
!.....input the value cpt 
!
      xp = cpt(1)
      yp = cpt(2)
      zp = cpt(3)
!
!.....main loop: calculate moment by Trapezodial rule
!
      do 30 ipt = 1, nep
!
!........define intermediate variable
!
!........a1, a2 and a3
!
       a1 =  anode(1,ipt)
       a2 =  anode(2,ipt)
       a3 =  anode(3,ipt)

       dsj = dwjp(ipt)

       xi  = cjp(1,ipt)
       yi  = cjp(2,ipt)
       zi  = cjp(3,ipt)
         dx  = -1.0/a1
       dy  = -1.0/a2
       dz  = -1.0/a3

       r100 = (xi - xp)/a1
       r010 = (yi - yp)/a2
       r001 = (zi - zp)/a3
       r110 = r100*r010
       r011 = r010*r001
       r101 = r100*r001

       r200 = r100*r100
       r210 = r200*r010
       r111 = r110*r001
       r201 = r200*r001

       r020 = r010*r010
       r120 = r100*r020
       r021 = r020*r001

       r002 = r001*r001
       r012 = r010*r002
       r102 = r100*r002

       r220 = r200*r020
       r121 = r101*r020
       r211 = r200*r011

       r022 = r020*r002
       r112 = r110*r002

       r202 = r200*r002
!
!  unfinidhes ??
!
       xx = abs(r100)
       yy = abs(r010)
       zz = abs(r001)

       if((xx .ge. 2.0) .or. (yy .ge. 2.0) .or. (zz .ge. 2.0)) go to 30

       call window3d(aw,awdx,awdy,awdz,    &
                      awdxx,awdyy,awdzz,   &
                      awdxy,awdyz,awdzx,   &
             xp,yp,zp,xi,yi,zi,a1,a2,a3)

       aw    = aw*dsj
!
!-----------------------------------
!
       am000 = am000  +  aw
       am100 = am100  +  r100*aw
       am010 = am010  +  r010*aw
       am001 = am001  +  r001*aw
       am110 = am110  +  r110*aw
       am011 = am011  +  r011*aw
       am101 = am101  +  r101*aw

       am200 = am200  +  r200*aw
       am210 = am210  +  r210*aw
       am111 = am111  +  r111*aw
       am201 = am201  +  r201*aw

       am020 = am020  +  r020*aw
       am120 = am120  +  r120*aw
       am021 = am021  +  r021*aw

       am002 = am002  +  r002*aw
       am012 = am012  +  r012*aw
       am102 = am102  +  r102*aw

       am220 = am220  +  r220*aw
       am121 = am121  +  r121*aw
       am211 = am211  +  r211*aw

       am022 = am022  +  r022*aw
       am112 = am112  +  r112*aw

       am202 = am202  +  r202*aw
!
!------------------------------------
!
       awdx = awdx*dsj
       awdy = awdy*dsj
       awdz = awdz*dsj
!
! (1)
!
         am000dx = am000dx  +  awdx
         am100dx = am100dx  +  r100*awdx  + dx*aw
         am010dx = am010dx  +  r010*awdx
       am001dx = am001dx  +  r001*awdx
       am110dx = am110dx  +  r110*awdx  + dx*r010*aw
       am011dx = am011dx  +  r011*awdx
       am101dx = am101dx  +  r101*awdx  + dx*r001*aw

         am200dx = am200dx  +  r200*awdx  + 2.0*dx*r100*aw
         am210dx = am210dx  +  r210*awdx  + 2.0*dx*r110*aw
       am111dx = am111dx  +  r111*awdx  + dx*r011*aw
       am201dx = am201dx  +  r201*awdx  + 2.0*dx*r101*aw

         am020dx = am020dx  +  r020*awdx
         am120dx = am120dx  +  r120*awdx  + dx*r020*aw
       am021dx = am021dx  +  r021*awdx

         am002dx = am002dx  +  r002*awdx
       am012dx = am012dx  +  r012*awdx
       am102dx = am102dx  +  r102*awdx  + dx*r002*aw

       am220dx = am220dx  +  r220*awdx  + 2.0*dx*r120*aw
       am121dx = am121dx  +  r121*awdx  + dx*r021*aw
       am211dx = am211dx  +  r211*awdx  + 2.0*dx*r111*aw

       am022dx = am022dx  +  r022*awdx
       am112dx = am112dx  +  r112*awdx  + dx*r012*aw

       am202dx = am202dx  +  r202*awdx  + 2.0*dx*r102*aw
!
! (2)
!
         am000dy = am000dy  +  awdy
         am100dy = am100dy  +  r100*awdy
         am010dy = am010dy  +  r010*awdy  + dy*aw
       am001dy = am001dy  +  r001*awdy
       am110dy = am110dy  +  r110*awdy  + dy*r100*aw
       am011dy = am011dy  +  r011*awdy  + dy*r001*aw
       am101dy = am101dy  +  r101*awdy

         am200dy = am200dy  +  r200*awdy
         am210dy = am210dy  +  r210*awdy  + dy*r200*aw
       am111dy = am111dy  +  r111*awdy  + dy*r101*aw
       am201dy = am201dy  +  r201*awdy

         am020dy = am020dy  +  r020*awdy  + 2.0*dy*r010*aw
         am120dy = am120dy  +  r120*awdy  + 2.0*dy*r110*aw
         am021dy = am021dy  +  r021*awdy  + 2.0*dy*r011*aw

         am002dy = am002dy  +  r002*awdy
       am012dy = am012dy  +  r012*awdy  + dy*r002*aw
       am102dy = am102dy  +  r102*awdy

       am220dy = am220dy  +  r220*awdy  + 2.0*dy*r210*aw
       am121dy = am121dy  +  r121*awdy  + 2.0*dy*r111*aw
       am211dy = am211dy  +  r211*awdy  + dy*r201*aw

       am022dy = am022dy  +  r022*awdy  + 2.0*dy*r012*aw
       am112dy = am112dy  +  r112*awdy  + dy*r102*aw

       am202dy = am202dy  +  r202*awdy
!
! (3)
!
         am000dz = am000dz  +  awdz
         am100dz = am100dz  +  r100*awdz
         am010dz = am010dz  +  r010*awdz
       am001dz = am001dz  +  r001*awdz  + dz*aw
       am110dz = am110dz  +  r110*awdz
       am011dz = am011dz  +  r011*awdz  + dz*r010*aw
       am101dz = am101dz  +  r101*awdz  + dz*r100*aw

         am200dz = am200dz  +  r200*awdz
         am210dz = am210dz  +  r210*awdz 
       am111dz = am111dz  +  r111*awdz  + dz*r110*aw
       am201dz = am201dz  +  r201*awdz  + dz*r200*aw

         am020dz = am020dz  +  r020*awdz
         am120dz = am120dz  +  r120*awdz
       am021dz = am021dz  +  r021*awdz  + dz*r020*aw

         am002dz = am002dz  +  r002*awdz  + 2.0*dz*r001*aw
       am012dz = am012dz  +  r012*awdz  + 2.0*dz*r011*aw
       am102dz = am102dz  +  r102*awdz  + 2.0*dz*r101*aw

       am220dz = am220dz  +  r220*awdz
       am121dz = am121dz  +  r121*awdz  + dz*r120*aw
       am211dz = am211dz  +  r211*awdz  + dz*r210*aw

       am022dz = am022dz  +  r022*awdz  + 2.0*dz*r021*aw
       am112dz = am112dz  +  r112*awdz  + 2.0*dz*r111*aw

       am202dz = am202dz  +  r202*awdz  + 2.0*dz*r201*aw


  30  continue

!
!.....end of the main loop
!
!.....assemble the M matrix (gm(i,j))
!
      gm(1,1) = am000
      gm(1,2) = am100
      gm(1,3) = am010
      gm(1,4) = am001
      gm(1,5) = am110
      gm(1,6) = am011
      gm(1,7) = am101

      gm(2,1) = am100
      gm(2,2) = am200
      gm(2,3) = am110
      gm(2,4) = am101
      gm(2,5) = am210
      gm(2,6) = am111
      gm(2,7) = am201

      gm(3,1) = am010
      gm(3,2) = am110
      gm(3,3) = am020
      gm(3,4) = am011
      gm(3,5) = am120
      gm(3,6) = am021
      gm(3,7) = am111

      gm(4,1) = am001
      gm(4,2) = am101
      gm(4,3) = am011
      gm(4,4) = am002
      gm(4,5) = am111
      gm(4,6) = am012
      gm(4,7) = am102

      gm(5,1) = am110
      gm(5,2) = am210
      gm(5,3) = am120
      gm(5,4) = am111
      gm(5,5) = am220
      gm(5,6) = am121
      gm(5,7) = am211

      gm(6,1) = am011
      gm(6,2) = am111
      gm(6,3) = am021
      gm(6,4) = am012
      gm(6,5) = am121
      gm(6,6) = am022
      gm(6,7) = am112

      gm(7,1) = am101
      gm(7,2) = am201
      gm(7,3) = am111
      gm(7,4) = am102
      gm(7,5) = am211
      gm(7,6) = am112
      gm(7,7) = am202
!
!.....calculate the determinat det
!
      call gjinv(gm,gminv,nn,nd,det,flag)
!
!.....test convergence criteria
!
      if(det .le. zero) then
      print *, 'det =', det
      print *, 'STOP! the determinat det < 0 '
      stop
      else
      end if
!
!.....assemble the gminv(i,j)
!

!
!.....calculate the derivative of gm
!

      gmdx(1,1) = am000dx
      gmdx(1,2) = am100dx
      gmdx(1,3) = am010dx
      gmdx(1,4) = am001dx
      gmdx(1,5) = am110dx
      gmdx(1,6) = am011dx
      gmdx(1,7) = am101dx

      gmdx(2,1) = am100dx
      gmdx(2,2) = am200dx
      gmdx(2,3) = am110dx
      gmdx(2,4) = am101dx
      gmdx(2,5) = am210dx
      gmdx(2,6) = am111dx
      gmdx(2,7) = am201dx

      gmdx(3,1) = am010dx
      gmdx(3,2) = am110dx
      gmdx(3,3) = am020dx
      gmdx(3,4) = am011dx
      gmdx(3,5) = am120dx
      gmdx(3,6) = am021dx
      gmdx(3,7) = am111dx

      gmdx(4,1) = am001dx
      gmdx(4,2) = am101dx
      gmdx(4,3) = am011dx
      gmdx(4,4) = am002dx
      gmdx(4,5) = am111dx
      gmdx(4,6) = am012dx
      gmdx(4,7) = am102dx

      gmdx(5,1) = am110dx
      gmdx(5,2) = am210dx
      gmdx(5,3) = am120dx
      gmdx(5,4) = am111dx
      gmdx(5,5) = am220dx
      gmdx(5,6) = am121dx
      gmdx(5,7) = am211dx

      gmdx(6,1) = am011dx
      gmdx(6,2) = am111dx
      gmdx(6,3) = am021dx
      gmdx(6,4) = am012dx
      gmdx(6,5) = am121dx
      gmdx(6,6) = am022dx
      gmdx(6,7) = am112dx

      gmdx(7,1) = am101dx
      gmdx(7,2) = am201dx
      gmdx(7,3) = am111dx
      gmdx(7,4) = am102dx
      gmdx(7,5) = am211dx
      gmdx(7,6) = am112dx
      gmdx(7,7) = am202dx


      gmdy(1,1) = am000dy
      gmdy(1,2) = am100dy
      gmdy(1,3) = am010dy
      gmdy(1,4) = am001dy
      gmdy(1,5) = am110dy
      gmdy(1,6) = am011dy
      gmdy(1,7) = am101dy

      gmdy(2,1) = am100dy
      gmdy(2,2) = am200dy
      gmdy(2,3) = am110dy
      gmdy(2,4) = am101dy
      gmdy(2,5) = am210dy
      gmdy(2,6) = am111dy
      gmdy(2,7) = am201dy

      gmdy(3,1) = am010dy
      gmdy(3,2) = am110dy
      gmdy(3,3) = am020dy
      gmdy(3,4) = am011dy
      gmdy(3,5) = am120dy
      gmdy(3,6) = am021dy
      gmdy(3,7) = am111dy

      gmdy(4,1) = am001dy
      gmdy(4,2) = am101dy
      gmdy(4,3) = am011dy
      gmdy(4,4) = am002dy
      gmdy(4,5) = am111dy
      gmdy(4,6) = am012dy
      gmdy(4,7) = am102dy

      gmdy(5,1) = am110dy
      gmdy(5,2) = am210dy
      gmdy(5,3) = am120dy
      gmdy(5,4) = am111dy
      gmdy(5,5) = am220dy
      gmdy(5,6) = am121dy
      gmdy(5,7) = am211dy

      gmdy(6,1) = am011dy
      gmdy(6,2) = am111dy
      gmdy(6,3) = am021dy
      gmdy(6,4) = am012dy
      gmdy(6,5) = am121dy
      gmdy(6,6) = am022dy
      gmdy(6,7) = am112dy

      gmdy(7,1) = am101dy
      gmdy(7,2) = am201dy
      gmdy(7,3) = am111dy
      gmdy(7,4) = am102dy
      gmdy(7,5) = am211dy
      gmdy(7,6) = am112dy
      gmdy(7,7) = am202dy


      gmdz(1,1) = am000dz
      gmdz(1,2) = am100dz
      gmdz(1,3) = am010dz
      gmdz(1,4) = am001dz
      gmdz(1,5) = am110dz
      gmdz(1,6) = am011dz
      gmdz(1,7) = am101dz

      gmdz(2,1) = am100dz
      gmdz(2,2) = am200dz
      gmdz(2,3) = am110dz
      gmdz(2,4) = am101dz
      gmdz(2,5) = am210dz
      gmdz(2,6) = am111dz
      gmdz(2,7) = am201dz

      gmdz(3,1) = am010dz
      gmdz(3,2) = am110dz
      gmdz(3,3) = am020dz
      gmdz(3,4) = am011dz
      gmdz(3,5) = am120dz
      gmdz(3,6) = am021dz
      gmdz(3,7) = am111dz

      gmdz(4,1) = am001dz
      gmdz(4,2) = am101dz
      gmdz(4,3) = am011dz
      gmdz(4,4) = am002dz
      gmdz(4,5) = am111dz
      gmdz(4,6) = am012dz
      gmdz(4,7) = am102dz

      gmdz(5,1) = am110dz
      gmdz(5,2) = am210dz
      gmdz(5,3) = am120dz
      gmdz(5,4) = am111dz
      gmdz(5,5) = am220dz
      gmdz(5,6) = am121dz
      gmdz(5,7) = am211dz

      gmdz(6,1) = am011dz
      gmdz(6,2) = am111dz
      gmdz(6,3) = am021dz
      gmdz(6,4) = am012dz
      gmdz(6,5) = am121dz
      gmdz(6,6) = am022dz
      gmdz(6,7) = am112dz

      gmdz(7,1) = am101dz
      gmdz(7,2) = am201dz
      gmdz(7,3) = am111dz
      gmdz(7,4) = am102dz
      gmdz(7,5) = am211dz
      gmdz(7,6) = am112dz
      gmdz(7,7) = am202dz

!
!.....assign the value for b vector
!
      do i = 1, nn
         b(i) = gminv(1,i)
      enddo
!
!.....find the value for bd(2,3)
!
!.....compute M^(-1)dM/dx and M^(-1)dM/dy
!
      do i = 1,nn
       do j = 1,nn
          adx(i,j) = 0.00    
          ady(i,j) = 0.00    
          adz(i,j) = 0.00    
          do k = 1,nn
             adx(i,j) = adx(i,j) + gminv(i,k)*gmdx(k,j)
             ady(i,j) = ady(i,j) + gminv(i,k)*gmdy(k,j)
             adz(i,j) = adz(i,j) + gminv(i,k)*gmdz(k,j)
            enddo
         enddo
      enddo


      do i = 1,nn
         do j = 1,nn
          bd(1,i) = bd(1,i) - adx(i,j)*b(j)
          bd(2,i) = bd(2,i) - ady(i,j)*b(j)
          bd(3,i) = bd(3,i) - adz(i,j)*b(j)
         end do
      end do

      return
      end

!************************************************************************
      subroutine correct2dl(b,bd,cpt,cjp,dcjp,dwjp,nep,inf,ninf,maxconn)
!************************************************************************
!
!     This subroutine is to calculate the b vector and its
!     derivatives,for the moving least square reproducing kernel 
!     interpolation shape function.
!
!     This code is only offering b vector and its 1st derivatives
!     for the 2-D case (RGT1a).
!
!     The interpolation function is choosen as
!     
!     p = (1 , x , y )
!
!     Date: June, 1994
!
!     -arguments:
!
!  i      nep: the numbers of particles, therefore, the numbers
!              of elements too.
!
!  i      cjp(2,nep):  the particle's global coordinations;
!                      such as x(nep):= cjp(1,nep)
!                              y(nep):= cjp(2,nep)
!
!  i      dcjp(2,nep): the array that stores all the dilation
!                     vectors of the particles .
!                     The dilation parameter at X direction: dcjp(1,nep)
!                     The dilation parameter at Y direction: dcjp(2,nep)
!
!  i      dwjp(nep): the integration weight at each particle point
!
!
!  i      cpt(2): the point at where the shape function is
!                     evaluated.;
!                 xp:= cpt(1)
!                 yp:= cpt(2)
!
!    
!   l     ha1 : dilation parameter; scalar;
!
!   l     ha2 : dilation parameter; scalar;
!
!               we input the adjust constants ax, ay, and we construct
!               the dilation parameter in such a way: 
!               (for non-uniform mesh )
!
!                 ha1 = dcjp(1,j)
!                 ha2 = dcjp(2,j)
!
!    l   gm(3,3) : the basi! moment matrix;
!
!    l   gminv(3,3): the inverse of the gm, i.e. gm^{-1}
!
!    l   gmdx(3,3): the derivatives of gm: gm_x;
!
!    l   gmdy(3,3): the derivatives of gm: gm_y;
!
!    l   adx(3,3):  M^(-1) M_x
!
!    l   ady(3,3):  M^(-1) M_y
!
!
!  o     b(3)   :   the b vector, also the first row of gminv
!
!  o     bd(2,3):   the 1st derivatives of the b vector
!
!                   bd(1,3) :=  b_x
!                   bd(2,3) :=  b_y
!
!     The subroutine is only designed to generate 2-D
!     moving least square shape function based on linear polynomial,
!
!    l   p(3) := (1,x,y);
!
!     The mathematical formulation is as follows
!
!     b   = (1/det)* [ a11, - a12, a13 ];
!
!     b_x = - M^{-1} M_x b;
!
!     b_y = - M^{-1} M_y b;
!
!
!    There two different subroutines could be called by this subroutine:
!
!    cubic2d.f ----- the cubi! spline window function;
!  
!    quintw2.f ----- the fifth order spline window function.
!
!**************************************************************
!      implicit double precision (a-h,o-z)
!      include 'parameter.h'
!
!      integer n_support
!      dimension Lmap(mnsch)


	  real* 8 am00,am10,am01,am20,am11,am02
	  real* 8 am00dx,am10dx,am01dx,am20dx,am11dx,am02dx
	  real* 8 am00dy,am10dy,am01dy,am20dy,am11dy,am02dy
	  real* 8 xp,yp,ha1,ha2,dsj,xj,yj,dx,dy
	  real* 8 r10,r01,r20,r11,r02
	  real* 8 xx,yy,aw
	  real* 8 awdx,awdy
	  real* 8 a11,a12,a21,a13,a31,a22,a23,a32,a33,det,cdet

      real* 8 b(3),bd(2,3),cpt(2)
      real* 8 cjp(2,nep),dcjp(2,nep),dwjp(nep)
      real* 8 gminv(3,3),gmdx(3,3),gmdy(3,3)
      real* 8 adx(3,3),ady(3,3)
	  integer jpt
	  integer ninf, maxconn, inf(maxconn)
!
!.....set the initial value for moment and its derivatives:
!
      zero = 0.00

      am00 = 0.00
      am10 = 0.00
      am01 = 0.00
      am20 = 0.00
      am11 = 0.00
      am02 = 0.00

      am00dx = 0.00
      am10dx = 0.00
      am01dx = 0.00
      am20dx = 0.00
      am11dx = 0.00
      am02dx = 0.00

      am00dy = 0.00
      am10dy = 0.00
      am01dy = 0.00
      am20dy = 0.00
      am11dy = 0.00
      am02dy = 0.00

!.....set the initial value for all array: 

      do i = 1, 3
		 do  j = 1, 3
		    gminv(i,j)  = 0.00
			gmdx(i,j)   = 0.00
			gmdy(i,j)   = 0.00
         end do
		b(i)    = 0.00
		bd(1,i) = 0.00
		bd(2,i) = 0.00
      end do  

!.....input the value cpt 

      xp = cpt(1)
      yp = cpt(2)

!.....main loop: calculate moment by Trapezodial rule

      do 30 j = 1, ninf
	  jpt = inf(j)

!........define intermediate variable

!........ha1 and ha2

	 ha1 =  dcjp(1,jpt)
	 ha2 =  dcjp(2,jpt)

	 dsj = dwjp(jpt)

	 xj  = cjp(1,jpt)
	 yj  = cjp(2,jpt)
     dx  = -1.0/ha1
	 dy  = -1.0/ha2

	 r10 = (xj - xp)/ha1
	 r01 = (yj - yp)/ha2
	 r20 = r10*r10
     r11 = r10*r01
     r02 = r01*r01

     xx = dabs(r10)
	 yy = dabs(r01)
	 if((xx.ge.2.0) .or. (yy.ge.2.0)) go to 30

     call window(aw,awdx,awdy,awdxx,awdxy,awdyy,xp,yp,xj,yj,ha1,ha2)

     aw   = aw*dsj

     am00 = am00  +  aw
         am10 = am10  +  r10*aw
         am01 = am01  +  r01*aw
         am20 = am20  +  r20*aw
         am11 = am11  +  r11*aw
         am02 = am02  +  r02*aw

	 awdx = awdx*dsj
	 awdy = awdy*dsj

     am00dx  = am00dx  + awdx
	 am10dx  = am10dx  + dx*aw + r10*awdx
	 am01dx  = am01dx  + r01*awdx
	 am20dx  = am20dx  + 2.0*dx*r10*aw + r20*awdx
	 am11dx  = am11dx  + dx*r01*aw + r11*awdx
	 am02dx  = am02dx  + r02*awdx

     am00dy  = am00dy  + awdy
	 am10dy  = am10dy  + r10*awdy
	 am01dy  = am01dy  + dy*aw + r01*awdy
	 am20dy  = am20dy  + r20*awdy
	 am11dy  = am11dy  + r10*dy*aw + r11*awdy
	 am02dy  = am02dy  + 2.0*dy*r01*aw + r02*awdy

  30  continue


!.....end of the main loop

!.....assemble the cofactor matrices ( a(i,j))

      a11 =  am20*am02  - am11*am11
      a12 =  am10*am02  - am11*am01
      a21 =  a12
      a13 =  am10*am11  - am20*am01
      a31 =  a13
      a22 =  am00*am02  - am01*am01
      a23 =  am00*am11  - am10*am01
      a32 =  a23
      a33 =  am00*am20  - am10*am10

!.....calculate the determinat det

      det = am00*a11 - am10*a12 + am01*a13 

!.....test convergence criteria

      zero = 0.0d0
      if(det .le. zero) then
		print *, 'det =', det
		print *, 'STOP! the determinat det < 0 '
		print *, am00, am10, am01 
		print *, am10, am20, am11 
		print *, am01, am11, am02 
		stop
      end if

!.....assemble the gminv(i,j)

      cdet = 1.0/det
      gminv(1,1) =  a11*cdet
      gminv(1,2) = -a12*cdet
      gminv(2,1) = -a12*cdet
      gminv(1,3) =  a13*cdet
      gminv(3,1) =  a13*cdet
      gminv(2,2) =  a22*cdet
      gminv(2,3) = -a23*cdet
      gminv(3,2) = -a23*cdet
      gminv(3,3) =  a33*cdet

!.....calculate the derivative of gm

      gmdx(1,1)  = am00dx
      gmdx(1,2)  = am10dx
      gmdx(2,1)  = am10dx
      gmdx(1,3)  = am01dx
      gmdx(3,1)  = am01dx
      gmdx(2,2)  = am20dx
      gmdx(2,3)  = am11dx
      gmdx(3,2)  = am11dx
      gmdx(3,3)  = am02dx

      gmdy(1,1)  = am00dy
      gmdy(1,2)  = am10dy
      gmdy(2,1)  = am10dy
      gmdy(1,3)  = am01dy
      gmdy(3,1)  = am01dy
      gmdy(2,2)  = am20dy
      gmdy(2,3)  = am11dy
      gmdy(3,2)  = am11dy
      gmdy(3,3)  = am02dy

!.....assign the value for b vector

      b(1) = gminv(1,1)
      b(2) = gminv(1,2)
      b(3) = gminv(1,3)

!.....Find the value for bd(2,3)

!.....compute M^(-1)dM/dx and M^(-1)dM/dy

      do i = 1,3
		 do j = 1,3
			adx(i,j) = 0.0d0    
			ady(i,j) = 0.0d0    
			do k = 1,3
			 adx(i,j) = adx(i,j) + gminv(i,k)*gmdx(k,j)
			 ady(i,j) = ady(i,j) + gminv(i,k)*gmdy(k,j)
            enddo
         enddo
      enddo
!
      do i = 1,3
         do j = 1,3
			bd(1,i) = bd(1,i) - adx(i,j)*b(j)
			bd(2,i) = bd(2,i) - ady(i,j)*b(j)
         enddo
      enddo

  999 continue

      return
      end
