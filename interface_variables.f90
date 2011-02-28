!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!Variables for interface!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module interface_variables
  implicit none
  save


  integer :: nn_inter     !# of interface points
  real(8) :: sur_tension  !surface tension force
  real(8) :: den_inter !density of bubble
  real(8) :: vis_inter !viscousity of bubble
  integer :: maxmatrix  
  real(8) :: scale_inter(3) ! 
  real(8) :: shift_inter(3) 
  real(8) :: curv_bound   !max curvature


end module interface_variables
