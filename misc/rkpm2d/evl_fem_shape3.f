c
      subroutine evl_FEM_shape3(shape,xsi,eta)
c      
c*** this subroutine is to evl the shape function values for 4-node element
c
c			y	
c			^
c			|
c			|
c			|
c			|
c               4---------------3
c               |		|
c               |  	    	|
c               | 		|
c	        |		| -------------------> x
c	        |		|
c		|  	    	|
c	        |		|
c		1---------------2
c
c---------------------------------------------------------------------
c
      implicit none      
      include 'parameter.h'
      real*8 shape(maxNode)
      real*8 xsi,eta
c      
      shape(1) = xsi
      shape(2) = eta
      shape(3) = 1.0 - xsi - eta
c      
      return
      end	! ends evl_FEM_shape3
c
