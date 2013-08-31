!==========================================

!calculate curv numrically

!==========================================

subroutine get_curv_num_3D(x,xp,x_center,I_fluid_center,corr_Ip,&
			curv_n,curv_inter,norm_inter)

  use interface_variables
  use fluid_variables, only:nsd,ne,nn

  real(8) x(nsd),xp(nsd,maxmatrix),x_center(nsd,nn_center)
  real(8) I_fluid_center(nn_center),corr_Ip(maxmatrix)
  real(8) curv_n,curv_inter,norm_inter(nsd)

  real(8) eps
  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) x_c(nsd)
  real(8) curv_x, curv_y,curv_z
  real(8) dcurv_x,dcurv_y,dcurv_z
!  real(8) dcurvA,dcurvB,dcurvC,dcurvD,dcurv_x,dcurv_y
  real(8) curv_a,norm_a(nsd)

  eps=1.0e-5
  dI(:)=0.0
  x_c(1:nsd)=x(1:nsd)
  call get_indicator_derivative_3D(x_c,xp,x_center,I_fluid_center,corr_Ip,&
				II,dI,ddI,norm_a,curv_a)
curv_inter=curv_a
norm_inter(1:nsd)=norm_a(1:nsd)
!===================================
  x_c(1)=x_c(1)+eps
  call get_indicator_derivative_3D(x_c,xp,x_center,I_fluid_center,corr_Ip,&
                                II,dI,ddI,norm_a,curv_x)

  dcurv_x=(curv_x-curv_a)/eps
!write(*,*)'dcurv_x=',dcurv_x
!====================================
  x_c(1:nsd)=x(1:nsd)
  x_c(2)=x_c(2)+eps
  call get_indicator_derivative_3D(x_c,xp,x_center,I_fluid_center,corr_Ip,&
                                II,dI,ddI,norm_a,curv_y)

  dcurv_y=(curv_y-curv_a)/eps
!====================================
  x_c(1:nsd)=x(1:nsd)
  x_c(3)=x_c(3)+eps
  call get_indicator_derivative_3D(x_c,xp,x_center,I_fluid_center,corr_Ip,&
                                II,dI,ddI,norm_a,curv_z)

  dcurv_z=(curv_z-curv_a)/eps

  curv_n=sqrt(dcurv_x**2+dcurv_y**2+dcurv_z**2)*max_hg**2
end subroutine get_curv_num_3D  


