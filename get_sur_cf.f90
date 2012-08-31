subroutine get_sur_cf(x,x_inter,x_center,I_fluid,corr_Ip,I_fluid_center,sur_fluid,hg,I_solid)

  use fluid_variables, only:nsd,ne,nn,nen,den_liq
  use interface_variables
  use allocate_variables, only:fluid_domain,nn_fluid_domain
  use mpi_variables
  include 'mpif.h'

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,ne)
  real(8) I_fluid(nn),corr_Ip(maxmatrix),I_fluid_center(ne)
  real(8) sur_fluid(nsd,nn),sur_fluid_temp(nsd,nn)
  real(8) I_solid(nn)

  real(8) den_p,den_f
  real(8) eps,hg(ne)

  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) curv_a,norm_a(nsd)

  integer nn_loc,base,top,loc_index

  if(nn_fluid_domain.le.ncpus) then
    if(myid+1.le.nn_fluid_domain) then
      nn_loc=1
    else
      nn_loc=0
    end if
  else
    base=floor(real(nn_fluid_domain)/real(ncpus))
    top=nn_fluid_domain-base*ncpus
    if(myid+1.le.top) then
      nn_loc=base+1
    else
      nn_loc=base
    end if
  end if

  sur_fluid(:,:)=0.0
  sur_fluid_temp(:,:)=0.0
  den_p=0.5*(den_inter+den_liq)
  eps=0.05

  do loc_index=1,nn_loc
     i=fluid_domain(myid+1+(loc_index-1)*ncpus)
!if(I_solid(i).ge.0.4) goto 200
     if( (I_fluid(i).lt.(1.0-eps)) .and. (I_fluid(i).gt.eps)) then
      if(nsd==2) then
       call get_indicator_derivative_2D(x(:,i),x_inter,x_center,hg,I_fluid_center,corr_Ip, &
                                        II,dI,ddI,norm_a,curv_a)
      elseif(nsd==3) then
              call get_indicator_derivative_3D(x(:,i),x_inter,x_center,hg,I_fluid_center,corr_Ip, &
                                        II,dI,ddI,norm_a,curv_a)
      end if
       den_f=den_liq+(den_inter-den_liq)*I_fluid(i)
       sur_fluid_temp(1:nsd,i)=sur_tension*(-curv_a)*dI(1:nsd)*den_f/den_p
!write(*,*)sur_fluid_temp(1:nsd,i),i
     end if
200 continue
  end do
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(sur_fluid_temp(1,1),sur_fluid(1,1),nsd*nn,mpi_double_precision, &
                mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

end subroutine get_sur_cf

