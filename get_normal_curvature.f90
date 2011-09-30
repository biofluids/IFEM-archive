!=============================

subroutine get_normal_curvature(x_inter,x_center,I_fluid_center,corr_Ip,&
				norm_inter,curv_inter,hg,dcurv)

  use interface_variables
  use fluid_variables, only:nsd,ne,nn
  use mpi_variables
  include 'mpif.h'
  real(8) x_inter(nsd,maxmatrix),x_center(nsd,ne)
  real(8) I_fluid_center(ne), corr_Ip(maxmatrix)
  real(8) norm_inter(nsd,maxmatrix),curv_inter(maxmatrix)
  real(8) hg(ne),dcurv,curv_n(nn_inter)

  integer i
  real(8) x(nsd)
  integer nn_inter_loc,base,top
  real(8) norm_inter_temp(nsd,maxmatrix)
  real(8) curv_inter_temp(maxmatrix)
  real(8) curv_n_temp(nn_inter)
  integer loc_index

  norm_inter_temp(:,:)=0.0
  curv_inter_temp(:)=0.0
  curv_n_temp(:)=0.0
  norm_inter(:,:)=0.0
  curv_inter(:)=0.0
  curv_n(:)=0.0

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

  do i=1,nn_inter_loc 
     loc_index=myid+1+(i-1)*ncpus
     x(1:nsd)=x_inter(1:nsd,loc_index)
     if(nsd==3) then
     call get_curv_num_3D(x,x_inter,x_center,hg,I_fluid_center,corr_Ip,&
			curv_n_temp(loc_index),curv_inter_temp(loc_index),&
			norm_inter_temp(1:nsd,loc_index)) 
     else
     call get_curv_num_2D(x,x_inter,x_center,hg,I_fluid_center,corr_Ip,&
                             curv_n_temp(loc_index),curv_inter_temp(loc_index),&
			                             norm_inter_temp(1:nsd,loc_index))
     end if
  end do
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(curv_n_temp(1),curv_n(1),nn_inter,mpi_double_precision, &
		mpi_sum,0,mpi_comm_world,ierror)
  call mpi_reduce(norm_inter_temp(1,1),norm_inter(1,1),nsd*maxmatrix,mpi_double_precision, &
		mpi_sum,0,mpi_comm_world,ierror)
  call mpi_reduce(curv_inter_temp(1),curv_inter(1),maxmatrix,mpi_double_precision, &
		mpi_sum,0,mpi_comm_world,ierror)
  call mpi_bcast(curv_n(1),nn_inter,mpi_double_precision,0,mpi_comm_world,ierror)
  call mpi_bcast(norm_inter(1,1),nsd*maxmatrix,mpi_double_precision,0,mpi_comm_world,ierror)
  call mpi_bcast(curv_inter(1),maxmatrix,mpi_double_precision,0,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

!  do i=1,nn_inter
!     x(1:nsd)=x_inter(1:nsd,i)
!     call get_curv_num(x,x_inter,x_center,hg,infdomain,I_fluid_center,corr_Ip,&
!			curv_n(i),curv_inter(i),norm_inter(1:nsd,i))
!  end do

  dcurv=maxval(curv_n(1:nn_inter))
  if(myid==0) then
    write(*,*)'max curv=',maxval(abs(curv_inter(1:nn_inter)))
    write(*,*)'max dcurv=',dcurv
!    write(*,*)'curv=',curv_inter(1:nn_inter)
  end if

end subroutine get_normal_curvature
