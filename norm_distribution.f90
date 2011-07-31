!===============================
!get interface velocity
!=================================

subroutine norm_distribution(x_den,x_inter,vol_nn, &
			arc_inter,norm_inter,norm_fluid)

  use fluid_variables, only:nsd
  use denmesh_variables
  use interface_variables
  use mpi_variables
  include 'mpif.h'
  real(8) x_den(nsd,nn_den),x_inter(nsd,maxmatrix)
  real(8) vol_nn(nn_den)
  real(8) arc_inter(maxmatrix)
  real(8) norm_inter(nsd,maxmatrix)
  real(8) norm_fluid(nsd,nn_den),norm_fluid_temp(nsd,nn_den)

  integer i,j,icount,jcount
  real(8) dx(nsd),Sp,hs,temp

  real(8) M(nsd+1,nsd+1),B(nsd+1),P(nsd+1)
  real(8) vec(nsd+1)
  integer IP(nsd+1)  

!=====used for mpi
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

  norm_fluid(:,:)=0.0
  norm_fluid_temp(:,:)=0.0

!  do i=1,nn_inter
  do loc_index=1,nn_inter_loc
     i=myid+1+(loc_index-1)*ncpus
     M(:,:)=0.0
     B(:)=0.0
     P(:)=0.0
     P(1)=1.0
     vec(1)=1.0

     do j=1,nn_den
        
	dx(:)=abs(x_inter(:,i)-x_den(:,j))
	call B_Spline(dx,hsp,nsd,Sp)
	vec(2:nsd+1)=x_inter(:,i)-x_den(:,j)
        do icount=1,nsd+1
           do jcount=1,nsd+1
              M(icount,jcount)=M(icount,jcount)+vec(icount)*vec(jcount)*Sp/(hsp**nsd)*vol_nn(j)
           end do
        end do
     end do
     call DGESV(nsd+1,1,M,nsd+1,IP,P,nsd+1,INFO)
     B(:)=P(:)
     do j=1,nn_den
	dx(:)=abs(x_inter(:,i)-x_den(:,j))
	call B_Spline(dx,hsp,nsd,Sp)
	vec(2:nsd+1)=x_inter(:,i)-x_den(:,j)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount)
	end do
	norm_fluid_temp(:,j)=norm_fluid_temp(:,j)+norm_inter(:,i)*arc_inter(i)*temp*Sp/(hsp**nsd)*vol_nn(j)

!	sur_fluid_temp(:,j)=sur_fluid_temp(:,j)+arc_inter(i)*sur_tension*curv_inter(i)* &
!			norm_inter(:,i)*temp*Sp/(hsp**nsd)*vol_nn(j)*den_f/den_p
!	
     end do
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(norm_fluid_temp(1,1),norm_fluid(1,1),nsd*nn_den,mpi_double_precision, &
                mpi_sum,0,mpi_comm_world,ierror)
  call mpi_bcast(norm_fluid(1,1),nsd*nn_den,mpi_double_precision,0,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)



end subroutine norm_distribution
