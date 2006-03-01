!     
!     bar pressure and derivative 
!
subroutine r_sbpress(dxmj,ddxmj,xmj)
  use r_common
  use solid_variables, only: nsd_solid
  
  implicit none 

  real(8) :: xmj(3),dxmj(3,2*nsd_solid)
  real(8) :: ddxmj(3,2*nsd_solid,2*nsd_solid)

  integer :: i,j

  bpre=-rk*(xmj(3)-1.0d0)
  do i=1,2*nsd_solid
     dbpre(i)=-rk*dxmj(3,i)
     do j=1,2*nsd_solid
        ddbpre(i,j)=-rk*ddxmj(3,i,j)
     enddo
  enddo
  
  return
end subroutine r_sbpress
