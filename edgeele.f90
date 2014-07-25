subroutine edgeele(edge,mrng,neface,ne,bcel,n_bcel)
! Find the element index which containing part of edge and return as bcel
! edge is from 1 to 4 fro 2-D case

integer edge ! edge index what to find
integer mrng(neface,ne) ! boundary information
integer ne ! number of element
integer neface ! edge per element
integer bcel(n_bcel) ! returning results
integer n_bcel ! number of elements having part of edge

integer i
integer j
integer k

    k=0
    do i=1,ne
    	do j=1,neface
    	    if(mrng(j,i) == edge) then
    	        k=k+1
    	        bcel(k)=i
    	    endif
    	enddo
    enddo

end
