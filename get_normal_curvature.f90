!=====================================
!get normal and curvature considering the summation of weight
!=====================================

subroutine get_normal_curvature(x_inter,x_center,I_fluid_center,corr_Ip,infdomain_inter,norm_inter,curv_inter,hg)

  use interface_variables
  use fluid_variables, only:nsd,ne,nn
!  use mpi_variables

  real(8) x_inter(nsd,maxmatrix),x_center(nsd,ne)
  real(8) I_fluid_center(ne), corr_Ip(maxmatrix)
  integer infdomain_inter(maxmatrix)
  real(8) norm_inter(nsd,maxmatrix)
  real(8) curv_inter(maxmatrix)
  real(8) hg(ne)

  integer i,j
  real(8) Sp_sum_0d(2), Sp_sum_1d(2,nsd),Sp_sum_2d(2,3*(nsd-1))
  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) x(nsd), temp, A, B

  do i=1,nn_inter
     x(1:nsd)=x_inter(1:nsd,i)
     call get_weight_derivative(x,x_inter,x_center,infdomain_inter,hg,&
				Sp_sum_0d,Sp_sum_1d,Sp_sum_2d)
     call get_indicator_derivative(x,x_inter,x_center,infdomain_inter,hg, &
			I_fluid_center,corr_Ip,&
			Sp_sum_0d, Sp_sum_1d,Sp_sum_2d,II,dI,ddI)
     temp=sqrt(dI(1)**2+dI(2)**2)
     norm_inter(1:nsd,i)=-dI(1:nsd)/temp
     A=ddI(1)/temp+dI(1)*(-0.5)/(temp**3)*(2*dI(1)*ddI(1)+2*dI(2)*ddI(3))
     B=ddI(2)/temp+dI(2)*(-0.5)/(temp**3)*(2*dI(1)*ddI(3)+2*dI(2)*ddI(2))
!if(myid==0) then
!write(*,*)Sp_Sum_0d
!write(*,*)'xxxxxxxxxxxxxxxxxxxxxx'
!write(*,*)Sp_sum_1d(1:2)
!write(*,*)'xxxxxxxxxxxxxxxxxxxxx'
!write(*,*)Sp_sum_2d(1:3)
!end if
     curv_inter(i)=A+B
!     write(*,*)'i=',i,'curv=',curv_inter(i)
  end do
write(*,*)'max curvature=',maxval(abs(curv_inter(1:nn_inter)))

end subroutine get_normal_curvature
