!  ansys_mesh_translator.f90 
!
!  FUNCTIONS:
!	ansys_mesh_translator      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Ansys Mesh Translator
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program ansys_mesh_translator
	implicit none

	integer,parameter :: solid_output = 0 
	integer,parameter :: fluid_output = 1 

	integer,parameter :: cube = 0
	integer,parameter :: tet  = 1

	integer,parameter :: iread_nlist = 1 !... fluid & solid
	integer,parameter :: maxnod = 4987 

	integer,parameter :: iread_dlist = 0 !... solid only
	integer,parameter :: maxline_d = 8 !... number of the last written line in dlist.lis
	integer,parameter :: maxconstrain = ((maxline_d-3)-mod(maxline_d-3,22))/22*20/3 + (mod(maxline_d-3,22) - 1)/3  !... number of nodes (each 3 constraines) in dlist.lis

	integer,parameter :: iread_elist = 1 !... fluid & solid
	integer,parameter :: maxele = 26030
	integer,parameter :: iread_sfelis = 1 !... fluid
	integer,parameter :: maxlin = 8899  !... number of the last written line in SFELIS.lis  --> number of boundaries calculated
	integer,parameter :: maxbnd = ((maxlin-4)-mod(maxlin-4,23))/23*5 + (mod(maxlin-4,23) - 1)/4  !... number of boundary information in SFELIS.lis

	real,dimension(1:maxnod,1:3)    :: n_loc = 0
	integer,dimension(1:maxnod)     :: n_disp_constrain = 0
	integer,dimension(1:maxele,1:8) :: ele_con = 0
	integer,dimension(1:maxele,1:6) :: ele_face = 0
 


	character(len = 1) :: dummy

	integer :: i,j
	integer,parameter :: zero = 0 
	real,parameter    :: one = 1.0d0
		
	real :: force
	real,parameter :: pi = 3.1415926



	write(*,*) "Translate ANSYS List output to input for fem mesh"
	write(*,*) "-------------------------------------------------"
	write(*,*) " "
 
    open(unit=23, file='tout.dat',status="replace")
	!open(unit=24, file='mesh.info',status="replace")
	!write(*,*) " prepared output files..."




	!... read node information
	if (iread_nlist.eq.1) then
		call read_nlist(maxnod,n_loc)
	end if

	!... read displacement constraineds information
	if (iread_dlist.eq.1) then
		call read_dlist(maxnod,n_disp_constrain,maxline_d)
	end if

	!... read element connectivity information
	if (iread_elist.eq.1) then
		call read_elist(maxele,ele_con,cube,tet)
	end if

	!... read boundary information
	if (iread_sfelis.eq.1) then
		call read_sfelis(maxele,maxlin,maxbnd,ele_con,ele_face,cube,tet)
	end if




	!... write fluid output
	if (fluid_output.eq.1) then
		if ( iread_nlist.eq.1 ) then
			open(unit=20, file='mxyz.dat',status="replace")
			do i = 1,maxnod
				write(20,'(3F12.5)') n_loc(i,1:3)
			end do
			close(20)
			write(*,*) "       mxyz.dat written!"
		end if
		
		if (cube == 1) then
			if ( iread_elist.eq.1 ) then
				open(unit=20, file='mien.dat',status="replace")
				do i = 1,maxele
					write(20,'(8I6)') ele_con(i,1),ele_con(i,4),ele_con(i,3),ele_con(i,2),ele_con(i,5),ele_con(i,8),ele_con(i,7),ele_con(i,6)
				end do
				close(20)
				write(*,*) "       mien.dat written!"
			end if

			if ( iread_sfelis.eq.1 ) then
				open(unit=20, file='mrng.dat',status="replace")
				do i = 1,maxele
					write(20,'(6I4)') ele_face(i,1:6)
				end do
				close(20)
				write(*,*) "       mrng.dat written!"
			end if
		elseif (tet == 1) then
			if ( iread_elist.eq.1 ) then
				open(unit=20, file='mien.dat',status="replace")
				do i = 1,maxele
					write(20,'(4I6)') ele_con(i,1),ele_con(i,2),ele_con(i,3),ele_con(i,4)
				end do
				close(20)
				write(*,*) "       mien.dat written!"
			end if

			if ( iread_sfelis.eq.1 ) then
				open(unit=20, file='mrng.dat',status="replace")
				do i = 1,maxele
					write(20,'(4I6)') ele_face(i,1),ele_face(i,2),ele_face(i,3),ele_face(i,4)
				end do
				close(20)
				write(*,*) "       mrng.dat written!"
			end if
		endif
		write(*,*) " "
	end if
	
	!... write solid output
	if (solid_output.eq.1) then

		open(unit=20, file='coortable_coord_only.dat',status="replace")
		
		write(20,'(4I6,"  /nnd,numel,nnda,numela")') maxnod,maxele,maxnod,maxele
		write(20,*) " "
		do i = 1,maxele
			write(20,'(10I6)') i,ele_con(i,1),ele_con(i,4),ele_con(i,3),ele_con(i,2),ele_con(i,5),ele_con(i,8),ele_con(i,7),ele_con(i,6),zero
		end do
		write(20,*) " "
		do i = 1,maxnod
			write(20,'(I5,3F10.5,I3)') i,n_loc(i,1:3),zero
		end do
		write(20,*) " "
		!write(20,*) "0 0 0 0 /numskew,numgb,intnum,numct"

        if (iread_dlist == 1) then

		write(20,*) "0 1 0 0 /numskew,numgb,intnum,numct"
		write(20,*) " "
		write(20,*) " 111111 /ndirgb of (numgb): 111111:all sides,110111:fix x-dir,101111:fix y-dir"
		write(20,*) " "
		write(20,'(I5,"			/numdir of (numgb)")') maxconstrain
		write(20,*) " "

		   do i = 1,maxconstrain
			  write(20,'(I5,3f8.3)') n_disp_constrain(i),one,zero,zero
		   enddo
        endif

		write(*,*) " "
		write(*,*) " coortable_coord_only.dat written!"
		write(*,*) " "
	end if


end program ansys_mesh_translator
