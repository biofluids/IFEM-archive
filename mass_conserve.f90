

subroutine mass_conserve(x,x_inter,x_center,I_fluid_center,I_fluid,ien,corr_Ip,its)

  use fluid_variables
  use interface_variables
  use mpi_variables
  use run_variables,only:nts_start
  include 'mpif.h'
  

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,nn_center)
  real(8) I_fluid_center(nn_center),I_fluid(nn)
  integer ien(nen,ne)
  real(8) corr_Ip(nsd,maxmatrix)
  real(8) mass,Length

  integer nn_loc,base,top,loc_index
  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) curv_a,norm_a(nsd)
  real(8) x_inter_temp(nsd,maxmatrix)
  integer i,j,its
if(myid==0)write(*,*)'begin mass conserve'

  call get_correction_mf(x_inter,x_center,corr_Ip,I_fluid_center)

  call get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,I_fluid)

  call cal_mass(x,I_fluid,ien,nn,ne,nen,mass)

  if((its.le.5).or.(its==nts_start)) then
     mass0=mass
     if(myid==0)write(*,*)'initial mass = ',mass0
      goto 100
  end if
if(myid==0)write(*,*)'mass0=',mass0,'mass=',mass
if(abs(mass0-mass)/mass0.gt.1.0e-4) then
  mass=mass0+(mass-mass0)/abs(mass-mass0)*1.0e-4*mass0
if(myid==0)write(*,*)'new mass=',mass
end if
  call cal_Length(x_inter,mass,x,I_fluid,ien,Length)

  x_inter_temp(:,:)=0.0
  if(nn_inter.le.ncpus) then
    if(myid+1.le.nn_inter) then
      nn_loc=1
    else
      nn_loc=0
    end if
  else
    base=floor(real(nn_inter)/real(ncpus))
    top=nn_inter-base*ncpus
    if(myid+1.le.top) then
      nn_loc=base+1
    else
      nn_loc=base
    end if
  end if

  do loc_index=1,nn_loc
     i=myid+1+(loc_index-1)*ncpus
     call get_indicator_derivative_2D_1st(x_inter(:,i),x_inter,x_center,I_fluid_center,corr_Ip, &
     					II,dI,ddI,norm_a,curv_a)
     
     x_inter_temp(:,i)=x_inter(:,i)+(mass0-mass)/Length*norm_a(:)

  end do
  x_inter(:,:)=0.0
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(x_inter_temp(1,1),x_inter(1,1),nsd*maxmatrix,mpi_double_precision, &
                    mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

100 continue 
end subroutine mass_conserve



