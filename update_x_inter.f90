

subroutine update_x_inter(x,x_inter,vel_inter,vel_fluid,vol_nn,dt)

  use fluid_variables, only:nsd,nn,ne,nen
  use interface_variables
  use mpi_variables
  include 'mpif.h'

  real(8) x_inter(nsd,maxmatrix),x(nsd,nn)
  real(8) x_inter_temp(nsd,maxmatrix)
  real(8) vel_inter(nsd,maxmatrix),vel_fluid(nsd,nn)
  real(8) vel_inter_temp(nsd,maxmatrix)
  real(8) vol_nn(nn),dt

  real(8) xloc(nsd),xloc_ini(nsd),xloc_temp(nsd),vel_loc(nsd)
  real(8) R_K(nsd)

  integer i,j,icount,jcount,isd

  integer nn_inter_loc,base,top,loc_index

  if(nn_inter.le.ncpus) then
    if(myid+1.le.nn_inter) then
      nn_inter_loc=1
    else
      nn_inter_loc=0
    end if
  else
     base=floor(real(nn_inter)/real(ncpus))
     top=nn_inter-base*ncpus
     if(myid+1.le.top) then
        nn_inter_loc=base+1
     else
       nn_inter_loc=base
     end if
  end if
  
  x_inter_temp(:,:)=0.0
  vel_inter_temp(:,:)=0.0

  do loc_index=1,nn_inter_loc
     i=myid+1+(loc_index-1)*ncpus
     xloc_ini(:)=x_inter(:,i)
     call get_inter_vel(x,xloc_ini,vel_fluid,vel_loc,vol_nn)

     vel_inter_temp(:,i)=vel_loc(:)
     R_K(:)=vel_loc(:)*dt
     xloc(:)=xloc_ini(:)+1.0/6.0*R_K(:)

     xloc_temp(:)=xloc_ini(:)+0.5*R_K(:)
     call get_inter_vel(x,xloc_temp,vel_fluid,vel_loc,vol_nn)
     R_K(:)=vel_loc(:)*dt
     xloc(:)=xloc(:)+1.0/3.0*R_K(:)

     xloc_temp(:)=xloc_ini(:)+0.5*R_K(:)
     call get_inter_vel(x,xloc_temp,vel_fluid,vel_loc,vol_nn)
     R_K(:)=vel_loc(:)*dt
     xloc(:)=xloc(:)+1.0/3.0*R_K(:)

     xloc_temp(:)=xloc_ini(:)+R_K(:)
     call get_inter_vel(x,xloc_temp,vel_fluid,vel_loc,vol_nn)
     R_K(:)=vel_loc(:)*dt
     xloc(:)=xloc(:)+1.0/6.0*R_K(:)

     x_inter_temp(:,i)=xloc(:)

  end do
call mpi_barrier(mpi_comm_world,ierror)
     x_inter(:,:)=0.0
     vel_inter(:,:)=0.0
call mpi_allreduce(x_inter_temp(1,1),x_inter(1,1),nsd*maxmatrix,mpi_double_precision,mpi_sum, &
			mpi_comm_world,ierror)
call mpi_allreduce(vel_inter_temp(1,1),vel_inter(1,1),nsd*maxmatrix,mpi_double_precision,mpi_sum,&
			mpi_comm_world,ierror)
call mpi_barrier(mpi_comm_world,ierror)



end subroutine update_x_inter




