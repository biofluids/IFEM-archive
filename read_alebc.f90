subroutine read_alebc(node_alebc,nn_alebc)
integer nn_alebc
integer node_alebc(nn_alebc)
integer file
integer i

file = 21
open(file, FILE='ale_bc.in',STATUS='old')

do i=1,nn_alebc
	read(file,'(I8)') node_alebc(i)
end do
close(file)

return

end subroutine read_alebc
