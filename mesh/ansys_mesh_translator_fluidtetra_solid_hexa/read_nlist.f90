subroutine read_nlist(maxnod,n_loc)
implicit none    

integer,intent(in) :: maxnod
integer :: i_block,n_block,nn,i,nn_test

real,intent(out),dimension(1:maxnod,1:3) :: n_loc

character(len = 1) :: dummy

open(unit=10,file="nlist.lis",status="old",action="read")
write(*,*) " opened nlist.lis..."

nn = 0
read(10,*) dummy
read(10,*) dummy
read(10,*) dummy

n_block = (maxnod - mod(maxnod,20))/20

if (n_block.gt.0) then
	do i_block = 1,n_block
		do i = 1,20
			nn = nn + 1
			read(10,*) nn_test,n_loc(nn,1:3)
		end do
		if (i_block.ne.n_block) then
		read(10,*) dummy
		end if
	end do	 
endif

if (mod(maxnod,20).ne.0) then
	if (n_block.gt.0) then
		read(10,*) dummy
	end if
	do i_block = 1,mod(maxnod,20)
		nn = nn + 1
		read(10,*) nn_test,n_loc(nn,1:3)
	end do
end if
close(10)

return
end subroutine read_nlist