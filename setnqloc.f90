subroutine setnqloc(id,mn,nn_local,node_local,nq_local)
use fluid_variables, only: nn
implicit none

integer :: mn, id(mn,nn)
integer  nn_local
integer node_local(nn_local)
integer nq_local
integer inc_local
integer idf

nq_local=0

do inc_local=1,nn_local
	do idf=1,mn
		if(id(idf,node_local(inc_local)) .ne. 0) then
			nq_local=nq_local+1
		end if
	end do
end do

end subroutine setnqloc

