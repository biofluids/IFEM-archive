c
      subroutine evl_FEM_shape4_a_node(xsi,eta,inode,shapei)
c
c*** this subroutine cal the shape function of the 4-node elem only
c
      implicit none
      integer inode
      real*8 xsi,eta
      real*8 shapei
c
      if (inode.eq.1) then
         shapei = 0.25*(1. - xsi)*(1. - eta)
      elseif (inode.eq.2) then
         shapei = 0.25*(1. + xsi)*(1. - eta)
      elseif (inode.eq.3) then
         shapei = 0.25*(1. + xsi)*(1. + eta)
      elseif (inode.eq.4) then
         shapei = 0.25*(1. - xsi)*(1. + eta)
      else
         write(*,*) 'wrong value of inode in subroutine evl_FEM_shape4'
      endif
c
      return
      end	! ends evl_FEM_shape4_a_node
c
