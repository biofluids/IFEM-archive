!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   updated.fcm                                                          c
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine updatePML(dinc, d, dg, mn)
  use fluid_variables
  use pml_variables
  implicit none

  integer,intent(in) :: mn          ! ndf -> mn
  real(8) :: d(mn,nn),dg(mn,nn+nn_PML)
  real(8),intent(out) :: dinc(mn,nn+nn_PML)

  integer inpml

  dinc(1:mn,1:nn+nn_PML) = dg(1:mn,1:nn+nn_PML)

  d(1:mn,1:nn) = d(1:mn,1:nn) + dinc(1:mn,1:nn)

  if (nn_PML > 0) then
      do inpml=1,nn_PML
          qv(1:mn,nodePML(inpml))=qv(1:mn,nodePML(inpml))+dinc(1:mn,nn+inpml)
      enddo
  endif

return
end subroutine updatePML