!     
!     jacobian,bd,pd calculation
!     
subroutine r_bdpd_curr(xji)
  use solid_variables, only: nsd_solid,nen_solid
  use r_common
  implicit none

  real(8) :: xji(nsd_solid,nsd_solid)

  integer :: i,j,k
  real(8) :: dumcd

 !...compute the bd matrix-> strain displacement matrix
  do i=1,nen_solid
     do k=1,nsd_solid
        dumcd=0.0d0
        do j=1,nsd_solid
           dumcd=dumcd+xji(j,k)*r_p(j,i)
        enddo
        bd_curr(k,i)=dumcd
    enddo
  enddo

  return
end subroutine r_bdpd_curr
