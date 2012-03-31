

subroutine formid_spbc(id,nn_spbc,spbcnode)

  use fluid_variables, only:ndf,nn
  use mpi_variables


  integer id(ndf,nn)
  integer nn_spbc, spbcnode(nn_spbc)

  integer i,j

  do i=1,nn_spbc
     id(1,spbcnode(i))=0
!     id(2,spbcnode(i))=1
  end do

if(myid==0)write(*,*)'finish formid_spbc'
end subroutine formid_spbc
