

      subroutine evl_FEM_shpdev3(lnods,xm,xsi,eta,ie,shape,cartd)
      
c*** this subroutine is to evl the shape function and derivative of 
c    globe coords (x and y ) for 3-node element

c*** output parameters:
c    shape(1:maxNode)   : N
c    cartd(2,1:maxNode) : dN/dx,dN/dy where the x,y is the reference coord ( big X)


c					y	
c					^
c					|
c					|
c					|
c					|
c                                       3 (0,1)
c                                       |  \
c                                       |   \		
c                                       |    \	
c				        |_____\        ---------> x
c			               1       2
c                                              (1,0)	
c

      implicit none 
      include 'parameter.h'

      integer nnode

      integer ie
      integer lnods(maxNode,maxElem) 
      real*8 xm(2,maxNumnp)
      real*8 xsi,eta      
      real*8 shape(maxNode),cartd(2,maxNode)
      
      real*8 xjacm(2,2),deriv(2,maxNode),xjaci(2,2)
      real*8 djacb
      integer inode,jnode,idim,jdim

      real*8 s,t
      data nnode/3/

c
c***  shape function and their derivatives for 3 noded element
c

      s=xsi
      t=eta
      
      shape(1)=s
      shape(2)=t
      shape(3)=1.-s-t
      deriv(1,1)=1.
      deriv(1,2)=0.
      deriv(1,3)=-1.
      deriv(2,1)=0.
      deriv(2,2)=1.
      deriv(2,3)=-1.

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
  901     format(//,'program halted in subroutine evl_FEM_shpdev3',
     .    /,11x,
     .    22h zero or negative area,/,10x,16h element number ,i5)
          stop 3200
        endif 
c        
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
      
      
      end	! ends evl_FEM_shpdev3
