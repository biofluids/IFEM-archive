subroutine read_solidbc(ien_sbc)
! Read in solid interface edge information
use solid_variables, only : ne_sbc, nen_solid
implicit none

integer ien_sbc(ne_sbc,nen_solid+1)
! ien_sbc[elment index, flag*nen (if node is on the edge)]
integer i
integer file

file=40
  open(file, FILE="sbc_solid.in", STATUS="old",action="read")

	do i=1,ne_sbc
		if (nen_solid == 3) then
		read(file,100) ien_sbc(i,1:nen_solid+1)
		end if

                if (nen_solid == 4) then
                read(file,101) ien_sbc(i,1:nen_solid+1)
                end if

                if (nen_solid == 8) then
                read(file,102) ien_sbc(i,1:nen_solid+1)
                end if



	end do

100 format(4I8)
101 format(5I8)
102 format(9I8)

close(file)

return

end subroutine read_solidbc
