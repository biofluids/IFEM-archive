


      subroutine evl_FEM_shpdev4(lnods,xm,xsi,eta,ie,shape,cartd)
      
c*** this subroutine is to evl the shape function and derivative of 
c    globe coords (x and y ) for 4-node element

c*** output parameters:
c    shape(1:4)   : N
c    cartd(2,1:4) : dN/dx,dN/dy where the x,y is the reference coord ( big X)


c					y	
c					^
c					|
c					|
c					|
c					|
c                               4---------------3
c                               |		|
c                               |  	    	|
c                               | 		|
c				|		| -------------------> x
c				|		|
c				|  	    	|
c				|		|
c				1---------------2
c
c
      implicit none 
      include 'parameter.h'

      integer ie
      integer lnods(maxNode,maxElem) 
      real*8 xm(2,maxNumnp)
      real*8 xsi,eta      
      real*8 shape(maxNode),cartd(2,maxNode)
      
      real*8 xjacm(2,2),deriv(2,maxNode),xjaci(2,2)
      real*8 djacb
      integer inode,jnode,idim,jdim

      integer nnode
      real*8 s,t
      data nnode/4/
c
c***  shape function for 4 noded element
c
      s=xsi
      t=eta

      shape(1)=(1.-s)*(1.-t)*0.25
      shape(2)=(1.+s)*(1.-t)*0.25
      shape(3)=(1.+s)*(1.+t)*0.25
      shape(4)=(1.-s)*(1.+t)*0.25
c
c***  shape function derivatives
c
      deriv(1,1)=-(1.-t)*0.25
      deriv(1,2)=(1.-t)*0.25
      deriv(1,3)=(1.+t)*0.25
      deriv(1,4)=-(1.+t)*0.25
      
      deriv(2,1)=-(1.-s)*0.25
      deriv(2,2)=-(1.+s)*0.25
      deriv(2,3)=(1.+s)*0.25
      deriv(2,4)=(1.-s)*0.25
c

      do idim=1,2
        do jdim=1,2
           xjacm(idim,jdim)=0.0
           do jnode=1,nnode
              inode=lnods(jnode,ie)
              xjacm(idim,jdim)=xjacm(idim,jdim)
     &                        +deriv(idim,jnode)*xm(jdim,inode)
           enddo
        enddo
      enddo  
c
c***  calculate determinant and inverse of jacobian matrix
c
        djacb=xjacm(1,1)*xjacm(2,2)-xjacm(1,2)*xjacm(2,1)
        
        if ( djacb .le. 0. ) then
          write(*,901) ie
  901     format(//,36h program halted in subroutine jacob2,/,11x,
     .    22h zero or negative area,/,10x,16h element number ,i5)
          stop 3200
        endif 
        
        xjaci(1,1)=xjacm(2,2)/djacb
        xjaci(2,2)=xjacm(1,1)/djacb
        xjaci(1,2)=-xjacm(1,2)/djacb
        xjaci(2,1)=-xjacm(2,1)/djacb
        
c
c***  calculate cartesian detivatives
c
        do idim=1,2
           do jnode=1,nnode
              inode=lnods(jnode,ie)
              cartd(idim,jnode)=0.0
              do jdim=1,2
                 cartd(idim,jnode)=cartd(idim,jnode)
     &                     +xjaci(idim,jdim)*deriv(jdim,jnode)
              enddo
           enddo
        enddo
      
      
      end	! ends evl_FEM_shpdev4
