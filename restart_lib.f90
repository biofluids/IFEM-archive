module restart_lib
  implicit none
  save

  integer :: restart_klok
  integer :: restart_freq
  integer :: restart

  integer,parameter,private :: restart_u1 = 21, restart_u2 = 22

  integer,private :: nts_test1,nts_test2

  integer,private :: restart_unit
  character(len = 12),private :: restart_file

  character(len = 12),private,parameter :: restart_f1 = "restart1.bin"
  character(len = 12),private,parameter :: restart_f2 = "restart2.bin"


contains

subroutine restart_file_check
  implicit none

  integer :: ios1,ios2
  integer :: nts_test1,nts_test2

  write(*,*) "check restart files..."
  open(unit=restart_u1, file=restart_f1, status='old', form='unformatted', iostat = ios1)
  open(unit=restart_u2, file=restart_f2, status='old', form='unformatted', iostat = ios2)
  if ((ios1 /= 0).and.(ios2 /= 0)) then
     write(*,*) " no restart files available"
     stop
  endif

  nts_test1 = 0
  nts_test2 = 0

  if (ios1 == 0) then
     read(unit=restart_u1) nts_test1
     close(restart_u1)
  endif
  if (ios2 == 0) then
     read(unit=restart_u2) nts_test2
     close(restart_u2)
  endif

  if ((nts_test1 /= restart).and.(nts_test2 /= restart)) then
     write(*,*) " "
     write(*,*) " Choose different restart timestep!!!"
     write(*,*) " restart information choosen:                ",restart
     write(*,*) " restart information available for timestep: ",nts_test1
     if (nts_test2 > 0) write(*,*) "                                        and: ",nts_test2
     stop
  endif

  if (nts_test1 == restart) then
     restart_unit = restart_u1
     restart_file = restart_f1
  else
     restart_unit = restart_u2
     restart_file = restart_f2
  endif

end subroutine restart_file_check




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine restart_read_header(nnloc,neloc,nenloc,ne_surfloc,nen_surfloc,etypeloc,nnfaceloc,nefaceloc)
  implicit none

  integer,intent(out) :: nnloc,neloc,nenloc,ne_surfloc,nen_surfloc,etypeloc,nnfaceloc,nefaceloc

  integer :: idummy

  write(*,*) " read header from <",restart_file,">..."

  open(unit=restart_unit, file=restart_file, status='old', form='unformatted',readonly)

  read(restart_unit) idummy

  read(restart_unit) nnloc,neloc,nenloc,ne_surfloc,nen_surfloc,etypeloc,nnfaceloc,nefaceloc

  close(restart_unit)

end subroutine restart_read_header




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine restart_read_data(nts_start,nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface,   &
                              nn_s,ne_s,nen_s,nn_s1,ne_s1,n_s,nsurf,                     &
                              nsd,ndf,nsd_s,                                             &
                              tt,klok,                                                   &
                              ien,ien_surf,rng,x,xref,hg,d,bid,                          &
                              solid_fem_con,solid_surface,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel)
  implicit none

  integer,intent(out)   :: nts_start
  integer,intent(inout) :: nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface

  integer,intent(inout) :: nn_s,ne_s,nen_s,nn_s1,ne_s1,n_s,nsurf

  integer,intent(in)  :: nsd,ndf,nsd_s

  real(8),intent(out) :: tt
  integer,intent(out) :: klok

  integer,intent(out) :: ien(nen,ne)
  integer,intent(out) :: ien_surf(nen_surf,ne_surf)
  integer,intent(out) :: rng(neface,ne)
  real(8),intent(out) :: x(nsd,nn),xref(nsd,nn)
  real(8),intent(out) :: hg(ne)
  real(8),intent(out) :: d(ndf,nn)
  integer,intent(out) :: bid(nn)


  integer,intent(out),dimension(1:ne_s,1:nen_s) :: solid_fem_con
  integer,intent(out),dimension(1:ne_s,1:nsurf) :: solid_surface
  real(8),intent(out),dimension(1:nsd_s,1:nn_s) :: solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel


  open(unit=restart_unit, file=restart_file, status='old', form='unformatted')

 !...simulation data
  read(restart_unit) nts_start
 !...fluid header data
  read(restart_unit) nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface
 !...solid header data
  read(restart_unit) nn_s,ne_s,nen_s,nn_s1,ne_s1,n_s,nsurf
 !...time data
  read(restart_unit) tt,klok,restart_klok
 !...fluid data
  read(restart_unit) ien,ien_surf,rng,x,xref,hg,d,bid
 !...solid data
  read(restart_unit) solid_fem_con,solid_surface,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel

  close(restart_unit)


  write(*,*) " "
  write(*,'("restart information succesfully read for timestep: ",I5)') nts_start

end subroutine restart_read_data


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine restart_write_data(its,nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface,   &
                              nn_s,ne_s,nen_s,nn_s1,ne_s1,n_s,nsurf,                &
                              nsd,ndf,nsd_s,                                        &
                              tt,klok,                                              &
                              ien,ien_surf,rng,x,xref,hg,d,bid,                     &
                              solid_fem_con,solid_surface,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel)
  implicit none

  integer,intent(in) :: its
  integer,intent(in) :: nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface

  integer,intent(in) :: nn_s,ne_s,nen_s,nn_s1,ne_s1,n_s,nsurf

  integer,intent(in) :: nsd,ndf,nsd_s

  real(8),intent(in) :: tt
  integer,intent(in) :: klok

  integer,intent(in) :: ien(nen,ne)
  integer,intent(in) :: ien_surf(nen_surf,ne_surf)
  integer,intent(in) :: rng(neface,ne)
  real(8),intent(in) :: x(nsd,nn),xref(nsd,nn)
  real(8),intent(in) :: hg(ne)
  real(8),intent(in) :: d(ndf,nn)
  integer,intent(in) :: bid(nn)

  integer,intent(in),dimension(1:ne_s,1:nen_s) :: solid_fem_con
  integer,intent(in),dimension(1:ne_s,1:nsurf) :: solid_surface
  real(8),intent(in),dimension(1:nsd_s,1:nn_s) :: solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel


  if (mod(restart_klok,2) == 0) then
     restart_unit = restart_u1
     restart_file = restart_f1
  else
     restart_unit = restart_u2
     restart_file = restart_f2
  endif


  open(unit=restart_unit, file=restart_file, status='replace', form='unformatted')

 !...simulation data
  write(restart_unit) its
 !...fluid header data
  write(restart_unit) nn,ne,nen,ne_surf,nen_surf,etype,nnface,neface
 !...solid header data
  write(restart_unit) nn_s,ne_s,nen_s,nn_s1,ne_s1,n_s,nsurf
 !...time data
  write(restart_unit) tt,klok,restart_klok
 !...fluid data
  write(restart_unit) ien,ien_surf,rng,x,xref,hg,d,bid
 !...solid data
  write(restart_unit) solid_fem_con,solid_surface,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel

  close(restart_unit)


  write(*,*) " restart file <",restart_file,"> written..."

end subroutine restart_write_data




end module restart_lib