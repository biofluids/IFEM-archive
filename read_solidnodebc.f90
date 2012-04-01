subroutine read_solidnodebc(node_sbc)
! Read in solid interface edge information
use solid_variables, only : nn_sbc_1
implicit none

integer node_sbc(nn_sbc_1)
! ien_sbc[elment index, flag*nen (if node is on the edge)]
integer i
integer file

file=40
  open(file, FILE="sbcnode_solid.in", STATUS="old",action="read")

	do i=1,nn_sbc_1
		read(file,100) node_sbc(i)
	end do

100 format(I8)
close(file)

return

end subroutine read_solidnodebc
