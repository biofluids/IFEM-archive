!==========================================!
!find out arclength ,surface tension, indicator for fluid field
!==========================================!

subroutine get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip,hg, &
		I_fluid)

  use interface_variables
  use fluid_variables, only:nn,ne,nsd,den_liq
  use mpi_variables
  include 'mpif.h'
  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,ne)
  real(8) I_fluid_center(ne),corr_Ip(maxmatrix),hg(ne)
  real(8) I_fluid(nn),I_fluid_temp(nn)

  integer i,j,icount,jcount,ie,inl,node,isd
  real(8) dx(nsd),hsg,Sp
  real(8) dcurv,temp,Bsum

!====used for non-uniform mesh
  real(8) M(nsd+1,nsd+1),B(nsd+1),P(nsd+1)
  real(8) vec(nsd+1)
  integer IP(nsd+1)
!===used for mpi
  integer nn_loc,base,top,loc_index

  if(nn.le.ncpus) then
    if(myid+1.le.nn) then
      nn_loc=1
    else 
      nn_loc=0
    end if
  else
     base=floor(real(nn)/real(ncpus))
     top=nn-base*ncpus
     if(myid+1.le.top) then
        nn_loc=base+1
     else
        nn_loc=base
     end if
  end if

! get I_fluid
  I_fluid(1:nn)=0.0
  I_fluid_temp(1:nn)=0.0
!  do i=1,nn
!     do j=1,nn_inter
!        dx(:)=abs(x(:,i)-x_inter(:,j))
!        call B_Spline(dx,hsp,nsd,Sp)
!        I_fluid(i)=I_fluid(i)+corr_Ip(j)*Sp
!     end do
!  end do

  do loc_index=1,nn_loc
     i=myid+1+(loc_index-1)*ncpus
     M(:,:)=0.0
     B(:)=0.0
     P(:)=0.0
     P(1)=1.0
     vec(1)=1.0
     do j=1,ne
	hsg=hg(j)
	dx(:)=abs(x(:,i)-x_center(:,j))
	call B_Spline(dx,hsp,nsd,Sp)
	vec(2:nsd+1)=x(:,i)-x_center(:,j)
	do icount=1,nsd+1
	   do jcount=1,nsd+1
	      M(icount,jcount)=M(icount,jcount)+vec(icount)*vec(jcount)*Sp/(hsp**nsd)*(hsg**nsd)
	   end do
	end do
     end do
     call DGESV(nsd+1,1,M,nsd+1,IP,P,nsd+1,INFO)
     B(:)=P(:)
     do j=1,ne
	hsg=hg(j)
        dx(:)=abs(x(:,i)-x_center(:,j))
        call B_Spline(dx,hsp,nsd,Sp)

	vec(2:nsd+1)=x(:,i)-x_center(:,j)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount)
	end do
        I_fluid_temp(i)=I_fluid_temp(i)+I_fluid_center(j)*temp*Sp/(hsp**nsd)*(hsg**nsd)
     end do
     
     do j=1,nn_inter
	dx(:)=abs(x(:,i)-x_inter(:,j))
	call B_Spline(dx,hsp,nsd,Sp)
	vec(2:nsd+1)=x(:,i)-x_inter(:,j)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount)
	end do
	I_fluid_temp(i)=I_fluid_temp(i)+corr_Ip(j)*temp*Sp
     end do

     if(I_fluid_temp(i).gt.1.0) then
	I_fluid_temp(i)=1.0
     else if(I_fluid_temp(i).lt.0.0) then
	I_fluid_temp(i)=0.0
     end if
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(I_fluid_temp(1),I_fluid(1),nn,mpi_double_precision, &
		mpi_sum,0,mpi_comm_world,ierror)
  call mpi_bcast(I_fluid(1),nn,mpi_double_precision,0,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)


end subroutine get_fluid_property











