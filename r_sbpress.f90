!     
!     bar pressure and derivative 
!
subroutine r_sbpress(dxmj,ddxmj,xmj)
  use r_common
  implicit none 

  real*8 :: xmj(3),dxmj(3,6),ddxmj(3,6,6)

  integer :: i,j

  bpre=-rk*(xmj(3)-1.0d0)
  do i=1,6
     dbpre(i)=-rk*dxmj(3,i)
     do j=1,6
        ddbpre(i,j)=-rk*ddxmj(3,i,j)
     enddo
  enddo
  
  return
end subroutine r_sbpress
