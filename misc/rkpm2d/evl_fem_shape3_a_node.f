c
      subroutine evl_FEM_shape3_a_node(xsi,eta,inode,shapei)
c
c*** this subroutine cal the shape function of the 3-node elem only
c
      implicit none
      integer inode
      real*8 xsi,eta		! the first two area coordinate
      real*8 shapei
c
      if (inode .eq. 1) then
         shapei = xsi
      elseif (inode .eq. 2) then
         shapei = eta
      elseif (inode .eq. 3) then
         shapei = 1.0 - xsi - eta
      else
         write(*,*) 'wrong value of inode in subroutine evl_FEM_shape3'
      endif
c
      return      
      end	! ends evl_FEM_shape3_a_node
c
