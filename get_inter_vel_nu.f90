!===============================
!get interface velocity
!=================================

subroutine get_inter_vel(x,x_inter,vel_fluid,vel_inter,hg,vol_nn)

  use fluid_variables, only:nsd,nn,ne,nen
  use interface_variables
  use mpi_variables
  include 'mpif.h'

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix)
  real(8) vel_fluid(nsd,nn),vel_inter(nsd,maxmatrix),vel_inter_temp(nsd,maxmatrix)
  real(8) hg(ne)
  real(8) vol_nn(nn)

  integer i,j,icount,jcount
  real(8) dx(nsd),Sp,hs,temp

  real(8) M(nsd+1,nsd+1),B(nsd+1),P(nsd+1)
  real(8) vec(nsd+1)
  integer IP(nsd+1),INFO  

  real(8) hsp_temp

!=====used for mpi===
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

  
  vel_inter(:,:)=0.0
  vel_inter_temp(:,:)=0.0

  hsp_temp=hsp
!  hsp=2*hsp_temp
! increase influence domain for vel intepolation
!  do i=1,nn_inter
  do loc_index=1,nn_inter_loc
     i=myid+1+(loc_index-1)*ncpus
     M(:,:)=0.0
     B(:)=0.0
     P(:)=0.0
     P(1)=1.0
     vec(1)=1.0

     do j=1,nn
	dx(:)=abs(x_inter(:,i)-x(:,j))
	call B_Spline1(dx,hsp,nsd,Sp,INFO)
	if(INFO==1) then
	vec(2:nsd+1)=x_inter(:,i)-x(:,j)
        do icount=1,nsd+1
           do jcount=1,nsd+1
              M(icount,jcount)=M(icount,jcount)+vec(icount)*vec(jcount)*Sp/(hsp**nsd)*vol_nn(j)
           end do
        end do
	end if
     end do
     call DGESV(nsd+1,1,M,nsd+1,IP,P,nsd+1,INFO)
     B(:)=P(:)
     do j=1,nn
	dx(:)=abs(x_inter(:,i)-x(:,j))
	call B_Spline1(dx,hsp,nsd,Sp,INFO)
	if(INFO==1) then
	vec(2:nsd+1)=x_inter(:,i)-x(:,j)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount)
	end do
	vel_inter_temp(:,i)=vel_inter_temp(:,i)+vel_fluid(:,j)*temp*Sp/(hsp**nsd)*vol_nn(j)

!	vel_inter(:,i)=vel_inter(:,i)+vel_fluid(:,j)*Sp
	end if
     end do
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(vel_inter_temp(1,1),vel_inter(1,1),nsd*maxmatrix,mpi_double_precision, &
                mpi_sum,0,mpi_comm_world,ierror)
  call mpi_bcast(vel_inter(1,1),nsd*maxmatrix,mpi_double_precision,0,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

  hsp=hsp_temp
end subroutine get_inter_vel
