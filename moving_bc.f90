subroutine moving_bc(d,id,node_alebc,nn_alebc,meshvel)
! Set the nodes on moving boundary 
! v=meshvel
! id(1:nsd)=0
use fluid_variables, only: ndf,nsd,nn
real(8) d(ndf,nn)
integer id(ndf,nn)
integer node_alebc(nn_alebc)
integer nn_alebc
real(8) meshvel(nsd,nn)

integer i

do i=1,nn_alebc
	d(1:nsd,node_alebc(i))=meshvel(1:nsd,node_alebc(i))
	id(1:nsd,node_alebc(i))=0
end do



end subroutine moving_bc
