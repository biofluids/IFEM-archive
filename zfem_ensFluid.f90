!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! modified form io11.f file to generate 
! ensight fluid field file
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine zfem_ensFluid(d,klok)
  use r_common
  use solid_variables
  use fluid_variables
  implicit none

  real*8 :: d(ndf,nn)
  integer :: klok
  
  real*8 :: stress_f(6,nn),strain_f(6,nn)

  integer :: i,j,k,m,n
  integer :: ntem
  character*5  fileroot
  character*12 name_file1
  character*12 name_file2
  character*15 name_file3
  character*15 name_file4

  real*8 :: pave(nn_solid)

  integer,parameter :: ifileunit = 15

 !...set the fluids stress and strain = 0
  stress_f(1:6,1:nn)=0.0
  strain_f(1:6,1:nn)=0.0
  write(*,*) "calculate stresses for output"

 !...calculate pressure for the structure
  do i=1,nn_solid
     
	 tstress(1:6,i)=0.0d0
     tstrain(1:6,i)=0.0d0
     
     ntem=0
     pave(i)=0
     if (nsd_solid .ne. 0) then  !do not calculate if it is a point
        do j=1,ne_solid
           do m=1,nis
              if (nea(j,m) .eq. i) then
                 ntem=ntem+1
                 pave(i)=pave(i)+pre(1,j)
                 do n=1,6
                    tstress(n,i)=tstress(n,i)+cstr(n,j,1,1,1)
                    tstrain(n,i)=tstrain(n,i)+ge(n,j,1,1,1)
				 enddo
                 go to 541
              endif
		   enddo
541	    enddo
        do k=1,6
           tstress(k,i)=tstress(k,i)/ntem
           tstrain(k,i)=tstrain(k,i)/ntem
	    enddo
	    pave(i)=pave(i)/ntem
	 endif
  enddo

  if (klok .eq. 0) then
     write(fileroot, '(a5)') '00000'
  elseif (klok .lt. 10) then
     write(fileroot, '(a4,i1)') '0000',klok
  elseif (klok .lt. 100) then
     write(fileroot, '(a3,i2)') '000' ,klok
  elseif (klok .lt. 1000) then
     write(fileroot, '(a2,i3)') '00'  ,klok
  elseif (klok .lt. 10000) then
     write(fileroot, '(a1,i4)') '0'   ,klok
  elseif (klok .lt. 100000) then    
     write(fileroot, '(i5)')   ''  ,klok
  else
     write(0,*) 'klok .ge. 100000: modify subroutine createfileroot'
     call exit(1)
  endif

  write(name_file1,'(A7, A5)')  'fem.vel', fileroot
  write(name_file2,'(A7, A5)')  'fem.pre', fileroot
  write(name_file3,'(A10, A5)')  'fem.stress', fileroot
  write(name_file4,'(A10, A5)')  'fem.strain', fileroot

  write(*,*) 'writing... ', name_file1
  write(*,*) 'writing... ', name_file2
  write(*,*) 'writing... ', name_file3
  write(*,*) 'writing... ', name_file4

 !...Write velocity output in ens_movie.vel*
  open(ifileunit, file=name_file1, form='formatted')
  write(ifileunit, '(A)') 	'structure and fluid field: velocity vector'
  write(ifileunit,110) (solid_vel(1,i),solid_vel(2,i),solid_vel(3,i),i=1,nn_solid),(d(1,i),d(2,i),d(3,i),i=1,nn)
  close(ifileunit)

 !...Write pressure output in ens_movie.pre*
  open(ifileunit, file=name_file2, form='formatted')
  write(ifileunit, '(A)') 'structure and fluid field: pressure'  
  write(ifileunit,110) (pave(i),i=1,nn_solid),(d(4,i),i=1,nn)
  close(ifileunit)

 !...Write stress output in ens_movie.stress*
  open(ifileunit, file=name_file3, form='formatted')
  write(ifileunit, '(A)') 'structure field: stress'  
  write(ifileunit,110) (tstress(1,i),tstress(2,i),tstress(3,i),            &
      	                tstress(4,i),tstress(5,i),tstress(6,i),i=1,nn_solid),   &
                       (stress_f(1,i),stress_f(2,i),stress_f(3,i),         &
                        stress_f(4,i),stress_f(5,i),stress_f(6,i),i=1,nn)  
  close(ifileunit)

 !...Write strain output in ens_movie.strain*

  open(ifileunit, file=name_file4, form='formatted')
  write(*, '(A)') 'structure field: strain'  
  write(ifileunit,110) (tstrain(1,i),tstrain(2,i),tstrain(3,i),            &
                        tstrain(4,i),tstrain(5,i),tstrain(6,i),i=1,nn_solid),   &
                       (strain_f(1,i),strain_f(2,i),strain_f(3,i),         &
                        strain_f(4,i),strain_f(5,i),strain_f(6,i),i=1,nn)  
  close(ifileunit)

110	format(6e12.5)

  return
end subroutine zfem_ensFluid
