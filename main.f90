!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   main.fcm                                                             c
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   initialization for the fluid and general variables is done in main.f
!   initialization for the solid variables (and input file reading) in hypo.f

program main
  use parseinput
  use adaptive_meshing_2
  use restart_lib
  implicit none

  integer,parameter :: new_mesh = 0

 !...set standard values 
  call initialize  ! for fluids

 !...read configuration files
  call parseinput_fluid  !...reading fluid.in
  call parseinput_solid  !...reading coortable.in

  call nondimension

 !...echos input data
  call echoinput

 !...if a remeshing step is executed, transfer data from old to new mesh
  if (new_mesh == 1) then
     call exchange_mesh
     stop
  endif

 !...switch to main routine    
  call hypo

end program main