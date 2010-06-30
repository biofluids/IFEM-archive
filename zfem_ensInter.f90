!!!!!!!write interface output!!!!!!!!!!!!!!!!

subroutine zfem_ensInter(I_fluid,norm_fluid,norm_inter,Ic_inter,klok)
  use solid_variables
  use fluid_variables,only:nn,nsd
  use interface_varibales
  use run_variables, only: its

  real(8) I_fluid(nn)
  real(8) norm_fluid(nsd,nn)
  real(8) norm_inter(maxmatrix)
  real(8) tmp_solid(nn_solid)
  real(8) tmp_inter(nn_inter)
  real(8) Ic_inter
  integer klok
  integer in,i
  character(len=7)  :: fileroot
  character(len=14) :: name_file1
  character(len=14) :: name_file2
  integer ifileunit = 17

  tmp_solid(:) = 0.0
  tmp_inter(:) = Ic_inter

  if (klok .eq. 0) then
     write(fileroot, '(a6)') '000000'
  elseif (klok .lt. 10) then
     write(fileroot, '(a5,i1)') '00000',klok
  elseif (klok .lt. 100) then
     write(fileroot, '(a4,i2)') '0000' ,klok
  elseif (klok .lt. 1000) then
     write(fileroot, '(a3,i3)') '000'  ,klok
  elseif (klok .lt. 10000) then
     write(fileroot, '(a2,i4)') '00'   ,klok
  elseif (klok .lt. 100000) then
     write(fileroot, '(a1,i5)') '0'   ,klok
  elseif (klok .lt. 1000000) then
     write(fileroot, '(i6)')   ''  ,klok
  else
     write(0,*) 'klok >= 1000000: modify subroutine createfileroot'
     call exit(1)
  endif

  write(name_file1,'(A8,  A6)')  'fem_ind', fileroot
  write(name_file2,'(A8,  A6)')  'fem_nor', fileroot

  if (nsd==3) then
     write(*,*)'writing...',name_file2
     open(ifileunit,file=name_file2,form='formatted')
     write(ifileunit, '(A)')  'stru,fluid and interface: normal vector'
     write(ifleunit,110) (tmp_solid(in), 0.0, 0.0, in=1,nn_solid), &
	   (norm_fluid(1,in),norm_fluid(2,in),norm_fluid(3,in),in=1,nn), &
	   (norm_inter(1,in),norm_inter(2,in),norm_inter(3,in),in=1,nn_inter)
     close(ifileunit)
  else if (nsd == 2) then
     write(*,*)'writing...',name_file2
     open(ifileunit,file=name_file2,form='formatted')
     write(ifileunit, '(A)')  'stru,fluid and interface: normal vector'
     write(ifleunit,110) (tmp_solid(in), 0.0, 0.0, in=1,nn_solid), &
           (norm_fluid(1,in),norm_fluid(2,in),0.0,in=1,nn), &
           (norm_inter(1,in),norm_inter(2,in),0.0,in=1,nn_inter)
     close(ifileunit)
  end if


     write(*,*) 'writing...',name_file1
     open(ifileunit,file=name_file1,form='formatted')
     write(ifileunit, '(A)') 'stru,fluid and interface: indicator'
     write(ifileunit, 110) (tmp_solid(in),in=1,nn_solid), &
		(I_fluid(in),in=1,nn),(tmp_inter(in),in=1,nn_inter)
     close(ifileunit)
110  format(6e12.5)

end subroutine zfem_ensInter
























  


















