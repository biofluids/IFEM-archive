subroutine read_elist(maxele,ele_con,cube,tet)
implicit none 

integer,intent(in) :: cube,tet,maxele
integer :: e_block,i_block,i,ne,id,ne_test


integer,dimension(1:maxele,1:8),intent(out) :: ele_con



character(len = 1) :: dummy


open(unit=11,file="elist.lis",status="old")
	write(*,*) " opened elist.lis..."
	ne = 0
	ele_con = 0
	read(11,*) dummy
	read(11,*) dummy

	e_block = (maxele - mod(maxele,20))/20
	!write(*,*) e_block
	if (e_block.gt.0) then
		do i_block = 1,e_block
			do i = 1,20
				ne = ne + 1
				read(11,*) ne_test,id,id,id,id,id,ele_con(ne,1:8)  !... !!!!! 4 or 5 dummy readings    CHECK in ELIST.lis
			end do
			if (i_block.ne.e_block) then
				read(11,*) dummy
			end if
		end do	 
	endif
	if (mod(maxele,20).ne.0) then
		if (e_block.gt.0) then
			read(11,*) dummy
		end if
		do i_block = 1,mod(maxele,20)
			ne = ne + 1
			read(11,*) ne_test,id,id,id,id,id,ele_con(ne,1:8)  !... !!!!! 4 or 5 dummy readings    CHECK in ELIST.lis
		end do
	end if

	close(11)
	
	if (tet == 1) then
		ele_con(:,4) = ele_con(:,5)
	endif
return
end