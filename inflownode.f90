subroutine inflownode(x,nn,nsd,x0,nn_inflow,ginflow,linflow)
! This subroutine is to set the set the node index for the nodes on the inflow boundary
! With the index, it is able to solve the linear equation for the inflow boundary
! ginflow is the node index from global to boudary equation
! linflow is the node index from boundary equation to global
! Note the # of nodes on the inflow boundary should be given to simpliy the problem
! linflow is also easy to be given in the pre-processing part

real(8) x(nsd,nn) ! x is the coordinates of fluid
integer nn ! # of nodes in fluid
integer nsd ! # of dimension
real(8) x0 ! boundary position, very easy model just assume node with x(1)=x0 are at the inflow boundary
integer nn_inflow ! # of nodes on the inflow boundary
integer ginflow(nn)
integer linflow(nn_inflow)
integer inn
integer i

i=0

do inn=1,nn
	if (x(1,inn)==x0) then
	i=i+1
	ginflow(inn)=i
	linflow(i)=inn
	end if
end do

end

