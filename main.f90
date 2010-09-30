!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   main.fcm                                                             c
!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   initialization for the fluid and general variables is done in main.f
!   initialization for the solid variables (and input file reading) in hypo.f

program main
  use parseinput
  implicit none

 !...set standard values 
  call initialize  ! for fluids

 !...read configuration files
  call parseinput_fluid  !...reading input_fluid.in
  call parseinput_solid  !...reading input_solid.in
  call parseinput_interface !...reading input_inter.in
!  call parseinput_meshcenter !...reading input_center_mesh.in

  call nondimension
 !...echos input data
  call echoinput

 !...switch to main routine    
  call hypo

end program main
