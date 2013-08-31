!==========================================!
!find out arclength ,surface tension, indicator for fluid field
!==========================================!

subroutine get_fluid_property(x,x_inter,x_center,I_fluid_center,corr_Ip, &
		I_fluid)

  use interface_variables
  use fluid_variables, only:nn,ne,nsd,den_liq,nen
  use mpi_variables
  use allocate_variables,only:nn_fluid_domain,fluid_domain
  include 'mpif.h'
  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,nn_center)
  real(8) I_fluid_center(nn_center),corr_Ip(maxmatrix)
  real(8) I_fluid(nn),II_fluid_temp(nn_fluid_domain),II_fluid(nn_fluid_domain)

  integer i,j,icount,jcount,ie,inl,node,isd,ii
  real(8) dx(nsd),hsg,Sp
  real(8) dcurv,temp,Bsum

!====used for non-uniform mesh
  real(8) M(nsd+1,nsd+1),B(nsd+1),P(nsd+1)
  real(8) vec(nsd+1)
  integer IP(nsd+1),INFO
!===used for mpi
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

!  I_fluid(fluid_domain(1:nn_fluid_domain))=0.0
!  I_fluid_temp(fluid_domain(1:nn_fluid_domain))=0.0
  II_fluid_temp(:)=0.0
  II_fluid(:)=0.0

  do loc_index=1,nn_loc
     ii=myid+1+(loc_index-1)*ncpus
     i=fluid_domain(myid+1+(loc_index-1)*ncpus)
     M(:,:)=0.0
     B(:)=0.0
     P(:)=0.0
     P(1)=1.0
     vec(1)=1.0
     do j=1,nn_center
	hsg=c_w(j)
	dx(:)=abs(x(:,i)-x_center(:,j))
	call B_Spline1(dx,hsp,nsd,Sp,INFO)
	if(INFO==1) then
	vec(2:nsd+1)=x(:,i)-x_center(:,j)
	do icount=1,nsd+1
	   do jcount=1,nsd+1
	      M(icount,jcount)=M(icount,jcount)+vec(icount)*vec(jcount)*Sp/(hsp**nsd)*(hsg**nsd)
	   end do
	end do
	end if
     end do
     call DGESV(nsd+1,1,M,nsd+1,IP,P,nsd+1,INFO)
     B(:)=P(:)
     do j=1,nn_center
	hsg=c_w(j)
        dx(:)=abs(x(:,i)-x_center(:,j))
        call B_Spline1(dx,hsp,nsd,Sp,INFO)
	if(INFO==1) then
	vec(2:nsd+1)=x(:,i)-x_center(:,j)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount)
	end do
        II_fluid_temp(ii)=II_fluid_temp(ii)+I_fluid_center(j)*temp*Sp/(hsp**nsd)*(hsg**nsd)
	end if
     end do
     
     do j=1,nn_inter
	dx(:)=abs(x(:,i)-x_inter(:,j))
	call B_Spline1(dx,hsp,nsd,Sp,INFO)
	if(INFO==1) then
	vec(2:nsd+1)=x(:,i)-x_inter(:,j)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount)
	end do
	II_fluid_temp(ii)=II_fluid_temp(ii)+corr_Ip(j)*temp*Sp
	end if
     end do

     if(II_fluid_temp(ii).gt.1.0) then
	II_fluid_temp(ii)=1.0
     else if(II_fluid_temp(ii).lt.0.0) then
	II_fluid_temp(ii)=0.0
     end if
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(II_fluid_temp(1),II_fluid(1),nn_fluid_domain,mpi_double_precision, &
		mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)
do icount=1,nn_fluid_domain
   node=fluid_domain(icount)
   I_fluid(node)=II_fluid(icount)
end do


end subroutine get_fluid_property











