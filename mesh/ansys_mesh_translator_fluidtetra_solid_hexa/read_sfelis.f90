subroutine read_sfelis(maxele,maxlin,maxbnd,ele_con,ele_face,cube,tet)
implicit none

integer,intent(in) :: maxele,maxlin,maxbnd,cube,tet
integer :: nb,ne,b_block,i_block,i,j,idummy,n1,n2,n3,n4,nface

integer,dimension(1:maxele,1:8),intent(in) :: ele_con
integer,dimension(1:maxele,1:6) :: ele_face
integer,dimension(1:maxbnd)     :: bound_type,bound_ne
integer,dimension(1:maxbnd,1:4) :: bound_ele_face
integer,dimension(1:8)          :: test_ele

real :: rbnd_type,real_dummy,test_bnd

character(len = 1) :: dummy


open(unit=12,file="sfelis.lis",status="old",action="read")
write(*,*) " opened sfelis.lis..."
read(12,*) dummy
read(12,*) dummy
b_block = ((maxlin-4)-mod(maxlin-4,23))/23
!write(*,*) b_block
if (b_block.gt.0) then
	do i_block = 1,b_block
		do i = 1,5
			nb = nb + 1 !... number of actual boundary
			do j = 1,4 !... loop over four nodes of an element face
				if (j.eq.1) then
					!write(*,*) nb,maxbnd
			        read(12,*) bound_ne(nb),idummy,bound_ele_face(nb,j),rbnd_type,real_dummy
					bound_type(nb) = int(rbnd_type)
				else
					read(12,*) bound_ele_face(nb,j),real_dummy,real_dummy
				end if
			    !write(*,*) ne,ne_test,ele_con(ne,1:8)
			end do
		end do
		if (i_block.ne.b_block) then
			read(12,*) dummy
		end if
	end do	 
endif


if (mod(maxlin-4,23).ne.0) then
	if (b_block.gt.0) then
		read(12,*) dummy
	end if
	do i_block = 1,((mod(maxlin-4,23)-1)/4)
		nb = nb + 1
		!write(*,*) nb
		do j = 1,4
			if (j.eq.1) then
			    read(12,*) bound_ne(nb),idummy,bound_ele_face(nb,j),rbnd_type,real_dummy
				bound_type(nb) = int(rbnd_type)
			else
				read(12,*) bound_ele_face(nb,j),real_dummy,real_dummy
			end if
		end do
	end do
end if
close(12)


if (cube == 1) then

	ele_face = 0
	do nb = 1,maxbnd
		!write(*,*) " nb = ",nb
		test_bnd = 1000000000*bound_ele_face(nb,1) + 1000000*bound_ele_face(nb,2) + 1000*bound_ele_face(nb,3) + bound_ele_face(nb,4)
		do ne = 1,maxele
			!write(*,*) " ne = ",ne
			if (bound_ne(nb).eq.ne) then
				do nface = 1,6
					!write(*,*) " nface = ",nface
						!... because of node exchange in elecon faces have to be changed, too  2--5   4--6
					if (nface.eq.1) then
						n1 = 1 ; n2 = 2 ; n3 = 3 ; n4 = 4
					else if (nface.eq.5) then
						n1 = 2
						n2 = 6
						n3 = 7
						n4 = 3
					else if (nface.eq.3) then
						n1 = 5
						n2 = 6
						n3 = 7
						n4 = 8
					else if (nface.eq.6) then
						n1 = 1
						n2 = 5
						n3 = 8
						n4 = 4
					else if (nface.eq.2) then
						n1 = 1
						n2 = 2
						n3 = 6
						n4 = 5
					else if (nface.eq.4) then
						n1 = 4
						n2 = 3
						n3 = 7
						n4 = 8
					end if
					test_ele(1) = 1000000000*ele_con(ne,n1) + 1000000*ele_con(ne,n2) + 1000*ele_con(ne,n3) + ele_con(ne,n4)
					test_ele(2) = 1000000000*ele_con(ne,n2) + 1000000*ele_con(ne,n3) + 1000*ele_con(ne,n4) + ele_con(ne,n1)
					test_ele(3) = 1000000000*ele_con(ne,n3) + 1000000*ele_con(ne,n4) + 1000*ele_con(ne,n1) + ele_con(ne,n2)
					test_ele(4) = 1000000000*ele_con(ne,n4) + 1000000*ele_con(ne,n1) + 1000*ele_con(ne,n2) + ele_con(ne,n3)
					test_ele(5) = 1000000000*ele_con(ne,n4) + 1000000*ele_con(ne,n3) + 1000*ele_con(ne,n2) + ele_con(ne,n1)
					test_ele(6) = 1000000000*ele_con(ne,n3) + 1000000*ele_con(ne,n2) + 1000*ele_con(ne,n1) + ele_con(ne,n4)
					test_ele(7) = 1000000000*ele_con(ne,n2) + 1000000*ele_con(ne,n1) + 1000*ele_con(ne,n4) + ele_con(ne,n3)
					test_ele(8) = 1000000000*ele_con(ne,n1) + 1000000*ele_con(ne,n4) + 1000*ele_con(ne,n3) + ele_con(ne,n2)
					do i =1,8
						if (test_bnd.eq.test_ele(i)) then
							ele_face(ne,nface) = bound_type(nb)
						end if
					end do
				end do
			end if
		end do
	end do
elseif (tet == 1) then
	ele_face = 0
	do nb = 1,maxbnd
		!write(*,*) " nb = ",nb
		test_bnd = 100000000*bound_ele_face(nb,1) + 10000*bound_ele_face(nb,2) + bound_ele_face(nb,3)
		do ne = 1,maxele
			!write(*,*) " ne = ",ne
			if (bound_ne(nb).eq.ne) then
				do nface = 1,4
					!write(*,*) " nface = ",nface
						!... because of node exchange in elecon faces have to be changed, too
					if (nface.eq.1) then
						n1 = 3 ; n2 = 2 ; n3 = 1
					else if (nface.eq.2) then
						n1 = 1 ; n2 = 2 ; n3 = 4
					else if (nface.eq.3) then
						n1 = 2 ; n2 = 3 ; n3 = 4
					else if (nface.eq.4) then
						n1 = 3 ; n2 = 1 ; n3 = 4
					end if
					test_ele(1) = 100000000*ele_con(ne,n1) + 10000*ele_con(ne,n2) + ele_con(ne,n3)
					test_ele(2) = 100000000*ele_con(ne,n2) + 10000*ele_con(ne,n3) + ele_con(ne,n1)
					test_ele(3) = 100000000*ele_con(ne,n3) + 10000*ele_con(ne,n1) + ele_con(ne,n2)
					test_ele(4) = 100000000*ele_con(ne,n3) + 10000*ele_con(ne,n2) + ele_con(ne,n1)
					test_ele(5) = 100000000*ele_con(ne,n2) + 10000*ele_con(ne,n1) + ele_con(ne,n3)
					test_ele(6) = 100000000*ele_con(ne,n1) + 10000*ele_con(ne,n3) + ele_con(ne,n2)
					do i =1,6
						if (test_bnd.eq.test_ele(i)) then
							ele_face(ne,nface) = bound_type(nb)
						end if
					end do
				end do
			end if
		end do
	end do
endif


return
end subroutine read_sfelis