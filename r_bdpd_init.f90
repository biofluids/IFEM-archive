!     
!     jacobian,bd,pd calculation
!     
subroutine r_bdpd_init(xji)
  use solid_variables, only: nen_solid
  use r_common
  implicit none

  real(8) :: xji(3,3)

  integer :: i,j,k
  real(8) :: dumcd

 !...compute the bd matrix-> strain displacement matrix
  do i=1,nen_solid
     do k=1,3
        dumcd=0.0d0
        do j=1,3
           dumcd=dumcd+xji(j,k)*r_p(j,i)
        enddo
        bd(k,i)=dumcd
    enddo
  enddo

  return
end subroutine r_bdpd_init
