subroutine equal_pa(x,y,ndf,nn,node_local,nn_local)
! let y=x for the nodes on its own proc
! y has to be clear to zero first

real(8) x(ndf,nn)
real(8) y(ndf,nn)
integer nn
integer ndf
integer node_local(nn_local)
integer nn_local

integer i
integer j
integer node

do i=1,nn_local
	node=node_local(i)
	do j=1,ndf
		y(j,node)=x(j,node)
	end do
end do

end subroutine equal_pa
