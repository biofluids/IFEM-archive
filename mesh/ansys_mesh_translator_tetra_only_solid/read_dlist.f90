subroutine read_dlist(maxconstrain,n_disp_constrain,maxline_d)
implicit none

integer,intent(in) :: maxconstrain,maxline_d
integer :: nb,b_block,i_block,i,j,idummy

integer :: nl
integer :: nn

integer,dimension(1:maxconstrain),intent(out) :: n_disp_constrain


real :: real_dummy

character(len = 1) :: dummy


nl = 0
nn = 0

open(unit=12,file="dlist.lis",status="old",action="read")
write(*,*) " opened dlist.lis..."
read(12,*) dummy
read(12,*) dummy
read(12,*) dummy
b_block = ((maxline_d-3)-mod(maxline_d-3,22))/22
if (b_block.gt.0) then
	do i_block = 1,b_block
		do i = 1,20
			nl = nl + 1 !... number of readline
			if (mod(nl,3) == 1) then 
				nn = nn + 1
			endif
			read(12,*) n_disp_constrain(nn),dummy,real_dummy,real_dummy
			write(*,*) n_disp_constrain(nn),dummy,real_dummy,real_dummy,nn,nl
		end do
		if (i_block.ne.b_block) then
			read(12,*) dummy
		end if
	end do	 
endif

if (mod(maxline_d-3,22).ne.0) then
	if (b_block.gt.0) then
		read(12,*) dummy
	end if
	do i = 1,mod(maxline_d-3,22)-2
		nl = nl + 1 !... number of readline
		if (mod(nl,3) == 1) then 
			nn = nn + 1
		endif
		read(12,*) n_disp_constrain(nn),dummy,real_dummy,real_dummy
		write(*,*) n_disp_constrain(nn),dummy,real_dummy,real_dummy,nn,nl
	end do
end if


close(12)


return
end subroutine read_dlist