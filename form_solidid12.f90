subroutine form_solidid12(kid,nsd,nn,ien_sbc,ne_sbc,nen,ne,ien)
! Base on the last column of ien_sbc to decide the type of BC
! 999 --> 1st type; -999 --> 2nd type
! set kid=0 at solid boundary nodes to put 1st type BC
implicit none

integer kid(nsd,nn)
integer nsd
integer nn
integer ien_sbc(ne_sbc,nen+2)
integer ne_sbc
integer nen
integer ne
integer ien(ne,nen)
integer i
integer j
integer node
integer ie

do i=1,ne_sbc
	if (ien_sbc(i,nen+2) == 999) then ! find element has 1st BC
		ie = ien_sbc(i,1) ! get the element index
		do j=1,nen
			if (ien_sbc(i,j+1) == 1) then ! decide the node index on the 1st BC
				node = ien(ie,j)
				kid(:,node) = 0 ! set kid = 0 at those points
			end if
		end do
	end if
end do

end
