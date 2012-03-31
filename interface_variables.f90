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
  real(8) :: maxdcurv
  real(8) :: vol_corr  !volume corretion
  real(8) :: total_length
  real(8) :: hsp
  real(8) :: rkpm_scale
  real(8) :: max_dcurv
  real(8) :: max_hg
  integer :: nbc
  real(8) :: mass0
end module interface_variables

