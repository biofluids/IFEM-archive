c     
c     i1,i2,i3,j1 j2 j3 and derivative
c     
      subroutine r_smaterj(wto,toc,xmi,xmj,dxmj,ddxmj)
      implicit real*8 (a-h,o-z) 
ccccccccccccccccccccccccccccccccccccccccccccccc
c     dxmj(i,j)      ------ i is xmj(i)
c     ddxmj(i,j,k)   ------ i is xmj(i)
c     1-- 11
c     2-- 22
c     3-- 12
c     4-- 33
cccccccccccccccccccccccccccccccccccccccccccccccc
      include 'r_common'
      dimension xmj(3),xmi(3),toc(3,3),
     $     dli(3,6),ddli(3,6,6),dxmj(3,6),ddxmj(3,6,6)
c     dimension otc(3,3),tb(3,1)
      cc=0.0d0
cccccccccccccccccccccccccccccccccccccccccccc

      do 20 i=1,2
         do 21 j=1,2
            cc=cc+toc(i,j)*toc(i,j)
c            otc(i,j)=toc(i,j)
   21    continue
   20 continue
c      otc(1,3)=0.0d0
c      otc(2,3)=0.0d0
c      otc(3,1)=0.0d0
c      otc(3,2)=0.0d0
c      otc(3,3)=1.0d0
c      tb(3,1)=0.0d0
c      tb(3,2)=0.0d0
c      tb(3,3)=0.0d0
c      call gaussj(otc,3,3,tb,1,1)
cccccccccccccccccccccccccccccccccccccccccccc
c     note the formulation is for 3-d
cccccccccccccccccccccccccccccccccccccccccccc
      cc=cc+1.0d0
      xmi(1)=toc(1,1)+toc(2,2)+1.0d0
      xmi(2)=0.5d0*(xmi(1)**2-cc)
      xmi(3)=toc(1,1)*toc(2,2)-toc(1,2)*toc(2,1)
      xmi1=dexp(-x13*dlog(xmi(3)))
      xmi2=dexp(-x23*dlog(xmi(3)))
      xmi4=dexp(-x43*dlog(xmi(3)))
      xmi5=dexp(-x53*dlog(xmi(3)))
      xmi7=dexp(-x73*dlog(xmi(3)))
      xmi8=dexp(-x83*dlog(xmi(3)))
c
      xmj(1)=xmi(1)*xmi1
      xmj(2)=xmi(2)*xmi2
      xmj(3)=dsqrt(xmi(3))
cccccccccccccccccccccccccccccccc
c     strain energy
cccccccccccccccccccccccccccccccc
      wto=rc1*(xmj(1)-3.0d0)+rc2*(xmj(2)-3.0d0)+
     1     0.5d0*rk*(xmj(3)-1.0d0)**2
cccccccccccccccccccccccccccccccc
      do 11 i=1,3
         do 12 j=1,4
            dli(i,j)=0.0d0
            do 15 m=1,4
               ddli(i,m,j)=0.0d0
   15       continue
   12    continue
   11 continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     derivative of I
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dli(3,1)=2.0d0*toc(2,2)
      dli(3,2)=2.0d0*toc(1,1)
      dli(3,3)=-2.0d0*toc(2,1)
      dli(3,4)=2.0d0*(toc(1,1)*toc(2,2)-toc(2,1)*toc(1,2))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      dli(3,4)=(otc(3,3)+otc(3,3))*xmi(3)
      dli(1,1)=2.0d0
      dli(1,2)=2.0d0
      dli(1,4)=2.0d0
      dli(2,1)=2.0d0+2.0d0*toc(2,2)
      dli(2,2)=2.0d0+2.0d0*toc(1,1)
      dli(2,3)=-2.0d0*toc(1,2)
      dli(2,4)=2.0d0*(toc(2,2)+toc(1,1))
      ddli(3,1,2)=4.0d0
      ddli(3,2,1)=4.0d0
c      ddli(3,3,3)=-4.0d0
      ddli(3,3,3)=-2.0d0
      ddli(2,1,2)=4.0d0
      ddli(2,2,1)=4.0d0
      ddli(2,3,3)=-2.0d0
c      ddli(2,3,3)=-4.0d0
      do 13 i=1,3
c
c     derivative of J
c
         dxmj(1,i)=dli(1,i)*xmi1-xmi(1)*xmi4*dli(3,i)*x13
         dxmj(2,i)=dli(2,i)*xmi2-xmi(2)*xmi5*dli(3,i)*x23
         dxmj(3,i)=0.5d0/dsqrt(xmi(3))*dli(3,i)
         do 14 j=1,3
            ddxmj(1,i,j)=-x13*xmi4*(dli(1,i)*dli(3,j)+
     $           dli(1,j)*dli(3,i)+ddli(3,i,j)*xmi(1))+
     $           x49*dli(3,j)*dli(3,i)*xmi7*xmi(1)
            ddxmj(2,i,j)=ddli(2,i,j)*xmi2-xmi5*x23*
     $           (dli(2,i)*dli(3,j)+dli(2,j)*dli(3,i)+
     $           ddli(3,i,j)*xmi(2))+dli(3,j)*
     $           dli(3,i)*xmi8*xmi(2)*x109         
            ddxmj(3,i,j)=-0.25d0*dli(3,i)*
     $           dli(3,j)/xmi(3)/dsqrt(xmi(3))+
     $           0.5d0*ddli(3,i,j)/dsqrt(xmi(3))
 14      continue
 13   continue
      dxmj(1,4)=dli(1,4)*xmi1-xmi(1)*xmi4*dli(3,4)*x13
      dxmj(2,4)=dli(2,4)*xmi2-xmi(2)*xmi5*dli(3,4)*x23
      dxmj(3,4)=0.5d0/dsqrt(xmi(3))*dli(3,4)
      return
      end





