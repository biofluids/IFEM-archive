!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!Variables for interface!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module interface_variables
  implicit none
  save


  integer :: nn_inter     !# of interface points
  real(8) :: sur_tension  !surface tension force
  integer :: maxmatrix  
  real(8) :: scale_inter(3) ! 
  real(8) :: shift_inter(3) 


end module interface_variables
