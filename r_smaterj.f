c     
c     i1,i2,i3,j1 j2 j3 and derivative
c     
      subroutine r_smaterj(wto,toc,xmi,xmj,dxmj,ddxmj)
      implicit real*8 (a-h,o-z) 
ccccccccccccccccccccccccccccccccccccccccccccccc
c     dxmj(i,j)      ------ i is xmj(i)=>invariance
c     ddxmj(i,j,k)   ------ i is xmj(i)
c-----------------------------------------------
c	2-D				3-D
c     1-- 11			1 -> 11
c     2-- 22			2 -> 22
c     3-- 12			3 -> 33
c     4-- 33			4 -> 23
c					5 -> 13
c					6 -> 12
cccccccccccccccccccccccccccccccccccccccccccccccc
      include 'r_common'
      dimension xmj(3),xmi(3),toc(3,3),
     1 dli(3,6),ddli(3,6,6),dxmj(3,6),ddxmj(3,6,6)

      cc=0.0d0
cccccccccccccccccccccccccccccccccccccccccccc
      do i=1,3
         do j=1,3
            cc=cc+toc(i,j)*toc(i,j)
	   enddo
	enddo

cccccccccccccccccccccccccccccccccccccccccccc
c     note the formulation is for 3-d
cccccccccccccccccccccccccccccccccccccccccccc

	xmi(1)=toc(1,1)+toc(2,2)+toc(3,3)
	xmi(2)=0.5d0*(xmi(1)**2-cc)
	xmi(3)=toc(1,1)*(toc(2,2)*toc(3,3)-toc(2,3)*toc(3,2))
     +   -toc(1,2)*(toc(2,1)*toc(3,3)-toc(3,1)*toc(2,3))
     +   +toc(1,3)*(toc(2,1)*toc(3,2)-toc(3,1)*toc(2,2))
c  -----
      xmi1=dexp(-x13*dlog(xmi(3)))
      xmi2=dexp(-x23*dlog(xmi(3)))
      xmi4=dexp(-x43*dlog(xmi(3)))
      xmi5=dexp(-x53*dlog(xmi(3)))
      xmi7=dexp(-x73*dlog(xmi(3)))
      xmi8=dexp(-x83*dlog(xmi(3)))
      xmj(1)=xmi(1)*xmi1
      xmj(2)=xmi(2)*xmi2
      xmj(3)=dsqrt(xmi(3))
	
cccccccccccccccccccccccccccccccc
c     strain energy
cccccccccccccccccccccccccccccccc
      wto=rc1*(xmj(1)-3.0d0)+rc2*(xmj(2)-3.0d0)+
     1     0.5d0*rk*(xmj(3)-1.0d0)**2
cccccccccccccccccccccccccccccccc
      do i=1,3
         do j=1,6
            dli(i,j)=0.0d0
            do m=1,6
               ddli(i,m,j)=0.0d0
		  enddo
	   enddo
	enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     first derivative of I => dI/dA
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  3-D case only
	dli(1,1)=2.0d0
	dli(1,2)=2.0d0
	dli(1,3)=2.0d0

	dli(2,1)=2.0d0*(toc(2,2)+toc(3,3))
	dli(2,2)=2.0d0*(toc(1,1)+toc(3,3))
	dli(2,3)=2.0d0*(toc(2,2)+toc(1,1))
	dli(2,4)=-(toc(2,3)+toc(3,2))
	dli(2,5)=-(toc(1,3)+toc(3,1))
	dli(2,6)=-(toc(1,2)+toc(2,1))

	dli(3,1)=2.0d0*(toc(2,2)*toc(3,3)-toc(2,3)*toc(3,2))
	dli(3,2)=2.0d0*(toc(1,1)*toc(3,3)-toc(1,3)*toc(3,1))
	dli(3,3)=2.0d0*(toc(1,1)*toc(2,2)-toc(1,2)*toc(2,1))
	dli(3,4)=2.0d0*(toc(1,2)*toc(3,1)-toc(1,1)*toc(3,2))
	dli(3,5)=2.0d0*(toc(2,1)*toc(3,2)-toc(3,1)*toc(2,2))
	dli(3,6)=2.0d0*(toc(3,1)*toc(2,3)-toc(2,1)*toc(3,3))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     second derivative of I => dI/dA
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  3-D case only
	ddli(2,1,2)=4.0d0
	ddli(2,1,3)=4.0d0
	ddli(2,2,1)=4.0d0
	ddli(2,2,3)=4.0d0
	ddli(2,3,1)=4.0d0
	ddli(2,3,2)=4.0d0
	ddli(2,4,4)=-2.0d0
	ddli(2,5,5)=-2.0d0
	ddli(2,6,6)=-2.0d0

	ddli(3,1,2)= 4.0d0*toc(3,3)
	ddli(3,1,3)= 4.0d0*toc(2,2)
	ddli(3,1,4)=-2.0d0*(toc(3,2)+toc(2,3))

	ddli(3,2,1)= 4.0d0*toc(3,3)
	ddli(3,2,3)= 4.0d0*toc(1,1)
	ddli(3,2,5)=-2.0d0*(toc(3,1)+toc(1,3))

	ddli(3,3,1)= 4.0d0*toc(2,2)
	ddli(3,3,2)= 4.0d0*toc(1,1)
	ddli(3,3,6)=-2.0d0*(toc(2,1)+toc(1,2))

	ddli(3,4,1)=-2.0d0*(toc(3,2)+toc(2,3))
	ddli(3,4,4)=-2.0d0*toc(1,1)
	ddli(3,4,5)= 2.0d0*(toc(1,2)+toc(2,1))
	ddli(3,4,6)= 2.0d0*(toc(3,1)+toc(1,3))

	ddli(3,5,2)=-2.0d0*(toc(3,1)+toc(1,3))
	ddli(3,5,4)= 2.0d0*(toc(2,1)+toc(1,2))
	ddli(3,5,5)=-2.0d0*toc(2,2)
	ddli(3,5,6)= 2.0d0*(toc(2,3)+toc(3,2))

	ddli(3,6,3)=-2.0d0*(toc(2,1)+toc(1,2))
	ddli(3,6,4)= 2.0d0*(toc(3,1)+toc(1,3))
	ddli(3,6,5)= 2.0d0*(toc(2,3)+toc(3,2))
	ddli(3,6,6)=-2.0d0*toc(3,3)

cccccccccccccccccccccccccccccccccccccccc
c     first derivative of J
cccccccccccccccccccccccccccccccccccccccc
      do i=1,6
         dxmj(1,i)=dli(1,i)*xmi1-xmi(1)*xmi4*dli(3,i)*x13
         dxmj(2,i)=dli(2,i)*xmi2-xmi(2)*xmi5*dli(3,i)*x23
         dxmj(3,i)=0.5d0/dsqrt(xmi(3))*dli(3,i)
	enddo

cccccccccccccccccccccccccccccccccccccccc
c     second derivative of J
cccccccccccccccccccccccccccccccccccccccc
	do i=1,6
         do j=1,6
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
	   enddo
	enddo
      return
      end