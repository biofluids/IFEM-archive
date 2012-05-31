

subroutine cal_vel_norm(x,x_inter,x_center,vel_fluid,vol_nn,hg,I_fluid_center,corr_Ip,nn_inter)

  use fluid_variables, only:nsd,nn,ne,nen
  use interface_variables, only:maxmatrix
  use mpi_variables
  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,ne)
  real(8) vel_fluid(nsd,nn),vol_nn(nn),hg(ne)
  real(8) I_fluid_center(ne),corr_Ip(maxmatrix)
  integer nn_inter

  integer i,j,icount,jcount
  real(8) temp
  real(8) II,dI(nsd),ddI(3*(nsd-1)),norm_p(nsd),curv_p
  real(8) xloc(nsd),vel_loc(nsd)
  temp=0.0
  do i=1,nn_inter
     xloc(:)=x_inter(1:nsd,i)
     call get_inter_vel(x,xloc,vel_fluid,vel_loc,vol_nn)
     call get_indicator_derivative_2D_1st(xloc,x_inter,x_center,hg,I_fluid_center,corr_Ip, &
			II,dI,ddI,norm_p,curv_p)
     temp=temp+(vel_loc(1)*norm_p(1)+vel_loc(2)*norm_p(2))**2
  end do

  temp=sqrt(temp)/real(nn_inter)

if(myid==0)write(*,*)'norm for normal inter velocity=',temp

end subroutine cal_vel_norm
