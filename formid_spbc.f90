

subroutine formid_spbc(id,nn_spbc,spbcnode)

  use fluid_variables, only:ndf,nn,lambda,nsd,f_slip
  use mpi_variables

  integer id(ndf,nn)
  integer nn_spbc, spbcnode(nn_spbc)

  integer i,j
if(f_slip==0) goto 2000

  do i=1,nn_spbc
     id(1,spbcnode(i))=0
!     id(2,spbcnode(i))=1
  end do
! used to treat corner points
  id(1:nsd,1)=0
  id(1:nsd,401)=0

if(myid==0)write(*,*)'finish formid_spbc'
2000 continue
end subroutine formid_spbc
