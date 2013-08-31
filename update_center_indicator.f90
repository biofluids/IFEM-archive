

subroutine update_center_indicator(x,x_inter,x_center,vel_fluid,vol_nn,dt,I_fluid_center,corr_Ip,I_solid)

  use fluid_variables, only:nsd,ne,nn,nen
  use interface_variables
  use allocate_variables, only:den_domain,ne_den_domain
  use mpi_variables
  include 'mpif.h'

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,nn_center)
  real(8) vel_fluid(nsd,nn),vol_nn(nn),dt
  real(8) I_fluid_center(nn_center),corr_Ip(maxmatrix),I_solid(nn)

  real(8) xloc(nsd),vel_loc(nsd)
  real(8) II,dI(nsd),ddI(3*(nsd-1)),norm_p(nsd),curv_p
  real(8) delta_I(ne_den_domain),delta_I_temp(ne_den_domain)
  
  integer i,j,k,icount,jcount,ie,inl,node,isd

  integer nn_loc,base,top,loc_index

  delta_I(:)=0.0
  delta_I_temp(:)=0.0

  if(ne_den_domain.le.ncpus) then
    if(myid+1.le.ne_den_domain) then
      nn_loc=1
    else
      nn_loc=0
    end if
  else
    base=floor(real(ne_den_domain)/real(ncpus))
    top=ne_den_domain-base*ncpus
    if(myid+1.le.top) then
      nn_loc=base+1
    else
      nn_loc=base
    end if
  end if

  do loc_index=1,nn_loc
     i=myid+1+(loc_index-1)*ncpus
     ie=den_domain(myid+1+(loc_index-1)*ncpus)
     xloc(:)=x_center(:,ie)
     if(nsd==2) then
	call get_indicator_derivative_2D_1st(xloc,x_inter,x_center, &
			I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
     else if(nsd==3) then
	call get_indicator_derivative_3D_1st(xloc,x_inter,x_center, &
			I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
     end if
     call get_inter_vel(x,xloc,vel_fluid,vel_loc,vol_nn,I_solid)
     do isd=1,nsd
        delta_I_temp(i)=I_fluid_center(ie)-vel_loc(isd)*dI(isd)*dt
     end do
!     write(*,*)'iedelta=',delta_I_temp(i)-I_fluid_center(ie)
     if(delta_I_temp(i).gt.1.0) delta_I_temp(i)=1.0
     if(delta_I_temp(i).lt.0.0) delta_I_temp(i)=0.0
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(delta_I_temp(1),delta_I(1),ne_den_domain,mpi_double_precision,mpi_sum, &
			mpi_comm_world,ierror)

  do ie=1,ne_den_domain
     I_fluid_center(den_domain(ie))=delta_I(ie)
  end do

end subroutine update_center_indicator



     
