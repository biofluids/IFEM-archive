      subroutine correct3dl(b,bd,cpt,cjp,anode,dwjp,nep,inf,ninf,maxconn)
c**************************************************************
c
c
      implicit none
c
      integer nep
      integer ninf, maxconn, iptc
      real*8 cjp(3,nep),anode(3,nep),dwjp(nep)
      real*8 gm(4,4),gminv(4,4)
      real*8 gmdx(4,4),gmdy(4,4),gmdz(4,4)
      real*8 adx(4,4),ady(4,4),adz(4,4)
      real*8 b(4),bd(3,4),cpt(3)
      integer  inf(maxconn)
c
      real*8 zero,flag,det
      real*8 am000,am100,am010,am001
      real*8       am200,am110,am101
      real*8             am020,am011
      real*8                   am002
c
      real*8 am000dx,am100dx,am010dx,am001dx
      real*8         am200dx,am110dx,am101dx
      real*8                 am020dx,am011dx
      real*8                         am002dx
c
      real*8 am000dy,am100dy,am010dy,am001dy
      real*8         am200dy,am110dy,am101dy
      real*8                 am020dy,am011dy
      real*8                         am002dy
c
      real*8 am000dz,am100dz,am010dz,am001dz
      real*8         am200dz,am110dz,am101dz
      real*8                 am020dz,am011dz
      real*8                         am002dz
c
      real*8 xi,yi,zi,dsj,dx,dy,dz
      real*8 xp,yp,zp,a1,a2,a3
      real*8 xx,yy,zz
      real*8 r100,r010,r001,r200,r110,r101
      real*8 r020,r011,r002
      real*8 aw,awdx,awdy,awdz
      real*8 awdxy,awdyz,awdzx
      real*8 awdxx,awdyy,awdzz
c
      integer i,j,k,ipt,nn,nd
c
c.....set the initial value for moment and its derivatives:
c
      zero  = 0.00
      nn    = 4
      nd    = 4
c
      am000 = 0.00
      am100 = 0.00
      am010 = 0.00
      am001 = 0.00
c
      am200 = 0.00
      am110 = 0.00
      am101 = 0.00
c
      am020 = 0.00
      am011 = 0.00
c
      am002 = 0.00
c
      am000dx = 0.00
      am100dx = 0.00
      am010dx = 0.00
      am001dx = 0.00
c
      am200dx = 0.00
      am110dx = 0.00
      am101dx = 0.00
c
      am020dx = 0.00
      am011dx = 0.00
c
      am002dx = 0.00
c
      am000dy = 0.00
      am100dy = 0.00
      am010dy = 0.00
      am001dy = 0.00
c
      am200dy = 0.00
      am110dy = 0.00
      am101dy = 0.00
c
      am020dy = 0.00
      am011dy = 0.00
c
      am002dy = 0.00
c
      am000dz = 0.00
      am100dz = 0.00
      am010dz = 0.00
      am001dz = 0.00
c
      am200dz = 0.00
      am110dz = 0.00
      am101dz = 0.00
c
      am020dz = 0.00
      am011dz = 0.00
c
      am002dz = 0.00
c
c.....set the initial value for all array: 
c
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
c
c.....input the value cpt 
c
      xp = cpt(1)
      yp = cpt(2)
      zp = cpt(3)
c
c.....main loop: calculate moment by Trapezodial rule
c
      do 30 iptc = 1, ninf
        ipt = inf(iptc)
c
c........define intermediate variable
c
c........a1, a2 and a3
c
        a1 =  anode(1,ipt)
        a2 =  anode(2,ipt)
        a3 =  anode(3,ipt)
c
        dsj = dwjp(ipt)
c
        xi  = cjp(1,ipt)
        yi  = cjp(2,ipt)
        zi  = cjp(3,ipt)
        dx  = -1.0/a1
        dy  = -1.0/a2
        dz  = -1.0/a3
c
        r100 = (xi - xp)/a1
        r010 = (yi - yp)/a2
        r001 = (zi - zp)/a3
c
        r200 = r100*r100
        r110 = r100*r010
        r101 = r100*r001
c
        r020 = r010*r010
        r011 = r010*r001
c
        r002 = r001*r001
c  
        xx = abs(r100)
        yy = abs(r010)
        zz = abs(r001)

        if((xx.ge.2.0) .or. (yy.ge.2.0) 
     &       .or. (zz .ge. 2.0)) go to 30
c
        call window3d(aw,awdx,awdy,awdz,
     &       awdxx,awdyy,awdzz,
     &       awdxy,awdyz,awdzx,
     &       xp,yp,zp,xi,yi,zi,a1,a2,a3)
c  
        aw    = aw*dsj
c  
c-----------------------------------
c  
        am000 = am000  +  aw
        am100 = am100  +  r100*aw
        am010 = am010  +  r010*aw
        am001 = am001  +  r001*aw
c  
        am200 = am200  +  r200*aw
        am110 = am110  +  r110*aw
        am101 = am101  +  r101*aw
c  
        am020 = am020  +  r020*aw
        am011 = am011  +  r011*aw
c  
        am002 = am002  +  r002*aw
c  
c------------------------------------
c  
        awdx = awdx*dsj
        awdy = awdy*dsj
        awdz = awdz*dsj
c  
c  (1)
c  
        am000dx = am000dx + awdx
        am100dx = am100dx + r100*awdx + dx*aw
        am010dx = am010dx + r010*awdx
        am001dx = am001dx + r001*awdx
c  
        am200dx = am200dx + r200*awdx + 2.0*dx*r100*aw
        am110dx = am110dx + r110*awdx + dx*r010*aw
        am101dx = am101dx + r101*awdx + dx*r001*aw    
c  
        am020dx = am020dx + r020*awdx
        am011dx = am011dx + r011*awdx
c  
        am002dx = am002dx + r002*awdx
c  
c  (2)
c  
        am000dy = am000dy + awdy
        am100dy = am100dy + r100*awdy
        am010dy = am010dy + r010*awdy + dy *aw
        am001dy = am001dy + r001*awdy
c  
        am200dy = am200dy + r200*awdy
        am110dy = am110dy + r110*awdy + dy*r100*aw  
        am101dy = am101dy + r101*awdy
c  
        am020dy = am020dy + r020*awdy + 2.0*dy*r010*aw 
        am011dy = am011dy + r011*awdy + dy*r001*aw    
c  
        am002dy = am002dy + r002*awdy
c  
c  (3)
c  
        am000dz = am000dz + awdz
        am100dz = am100dz + r100*awdz
        am010dz = am010dz + r010*awdz
        am001dz = am001dz + r001*awdz + dz*aw
c  
        am200dz = am200dz + r200*awdz
        am110dz = am110dz + r110*awdz
        am101dz = am101dz + r101*awdz + dz*r100*aw
c  
        am020dz = am020dz + r020*awdz
        am011dz = am011dz + r011*awdz + dz*r010*aw   
c  
        am002dz = am002dz + r002*awdz + 2.0*dz*r001*aw
c  
 30   continue
c  
c  
c.....end of the main loop
c  
c.....assemble the M matrix (gm(i,j))
c  
      gm(1,1) = am000
      gm(1,2) = am100
      gm(1,3) = am010
      gm(1,4) = am001
c  
      gm(2,1) = am100
      gm(2,2) = am200
      gm(2,3) = am110
      gm(2,4) = am101
c  
      gm(3,1) = am010
      gm(3,2) = am110
      gm(3,3) = am020
      gm(3,4) = am011
c  
      gm(4,1) = am001
      gm(4,2) = am101
      gm(4,3) = am011
      gm(4,4) = am002
c  
c.....calculate the determinat det
c  
      call gjinv(gm,gminv,nn,nd,det,flag)
c  
c.....test convergence criteria
c  
      if(det .le. zero) then
        print *, 'det =', det
        print *, 'STOP! the determinat det < 0 '
        stop
      else
      end if
c  
c.....assemble the gminv(i,j)
c  
c.....calculate the derivative of gm
c  
      gmdx(1,1) = am000dx
      gmdx(1,2) = am100dx
      gmdx(1,3) = am010dx
      gmdx(1,4) = am001dx
c  
      gmdx(2,1) = am100dx
      gmdx(2,2) = am200dx
      gmdx(2,3) = am110dx
      gmdx(2,4) = am101dx
c  
      gmdx(3,1) = am010dx
      gmdx(3,2) = am110dx
      gmdx(3,3) = am020dx
      gmdx(3,4) = am011dx
c  
      gmdx(4,1) = am001dx
      gmdx(4,2) = am101dx
      gmdx(4,3) = am011dx
      gmdx(4,4) = am002dx
c  
c  
      gmdy(1,1) = am000dy
      gmdy(1,2) = am100dy
      gmdy(1,3) = am010dy
      gmdy(1,4) = am001dy
c  
      gmdy(2,1) = am100dy
      gmdy(2,2) = am200dy
      gmdy(2,3) = am110dy
      gmdy(2,4) = am101dy
c  
      gmdy(3,1) = am010dy
      gmdy(3,2) = am110dy
      gmdy(3,3) = am020dy
      gmdy(3,4) = am011dy
c  
      gmdy(4,1) = am001dy
      gmdy(4,2) = am101dy
      gmdy(4,3) = am011dy
      gmdy(4,4) = am002dy
c  
c  
      gmdz(1,1) = am000dz
      gmdz(1,2) = am100dz
      gmdz(1,3) = am010dz
      gmdz(1,4) = am001dz
c  
      gmdz(2,1) = am100dz
      gmdz(2,2) = am200dz
      gmdz(2,3) = am110dz
      gmdz(2,4) = am101dz
c  
      gmdz(3,1) = am010dz
      gmdz(3,2) = am110dz
      gmdz(3,3) = am020dz
      gmdz(3,4) = am011dz
c  
      gmdz(4,1) = am001dz
      gmdz(4,2) = am101dz
      gmdz(4,3) = am011dz
      gmdz(4,4) = am002dz
c  
c.....assign the value for b vector
c  
      do i = 1, nn
        b(i) = gminv(1,i)
      enddo
c  
c.....find the value for bd(2,3)
c  
c.....compute M^(-1)dM/dx and M^(-1)dM/dy
c  
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
c  
c  
      do i = 1,nn
        do j = 1,nn
          bd(1,i) = bd(1,i) - adx(i,j)*b(j)
          bd(2,i) = bd(2,i) - ady(i,j)*b(j)
          bd(3,i) = bd(3,i) - adz(i,j)*b(j)
        end do
      end do
c  
c  
c  
 999  continue
c  
      return
      end
c  
      
      
      
      
      
      
      
      
      
