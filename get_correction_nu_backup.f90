!=====================================
! get correction term
!=====================================

subroutine get_correction_nu(x_inter,x_center,hg,corr_Ip,I_fluid_center)

  use interface_variables
  use fluid_variables, only:nn,ne,nsd
  use allocate_variables, only:center_domain,nn_center_domain
  use mpi_variables
  include 'mpif.h'
  real(8) x_inter(nsd,maxmatrix)
  real(8) x_center(nsd,ne)
  real(8) hg(ne)
  real(8) corr_Ip(maxmatrix)
  real(8) I_fluid_center(ne)
!  real(8) x(nsd,nn),I_fluid(nn)

  integer i,j,ie,icount,jcount
  real(8) dx(nsd),hs,Sp
  real(8) A(nn_inter,nn_inter),A_temp(nn_inter,nn_inter)
  real(8) BB(nn_inter), BB_temp(nn_inter)  !Ax=B
  integer IPIV(nn_inter), INFO

  real(8) M(nsd+1,nsd+1),B(nsd+1),P(nsd+1) ! used for nonuniform
  real(8) vec(nsd+1),hsg !M*B=P from rkpm
  integer IP(nsd+1)
  real(8) temp

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

  A(:,:)=0.0
  A_temp(:,:)=0.0
  BB(:)=0.0
  BB_temp(:)=0.0

!  do i=1,nn_inter
  do loc_index=1,nn_inter_loc
     i=myid+1+(loc_index-1)*ncpus
     hs=hsp
     M(:,:)=0.0
     B(:)=0.0
     P(:)=0.0
     P(1)=1.0
     vec(1)=1.0
     do j=1,nn_center_domain
	ie=center_domain(j)
	hsg=hg(ie)
        dx(:)=abs(x_inter(:,i)-x_center(:,ie))
	call B_Spline(dx,hs,nsd,Sp)
	vec(2:nsd+1)=x_inter(:,i)-x_center(:,ie)
	do icount=1,nsd+1
	   do jcount=1,nsd+1
	      M(icount,jcount)=M(icount,jcount)+vec(icount)*vec(jcount)*Sp/(hs**nsd)*(hsg**nsd)
	   end do
	end do
     end do
     call DGESV(nsd+1,1,M,nsd+1,IP,P,nsd+1,INFO)
     B(:)=P(:)

     do j=1,nn_inter
	dx(:)=abs(x_inter(:,i)-x_inter(:,j))
	call B_Spline(dx,hs,nsd,Sp)
	vec(2:nsd+1)=x_inter(:,i)-x_inter(:,j)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount)
	end do
	A_temp(i,j)=Sp*temp
     end do

     do j=1,nn_center_domain
	ie=center_domain(j)
	hsg=hg(ie)
	dx(:)=abs(x_inter(:,i)-x_center(:,ie))
	call B_Spline(dx,hs,nsd,Sp)

	vec(2:nsd+1)=x_inter(:,i)-x_center(:,ie)
	temp=0.0
	do icount=1,nsd+1
	   temp=temp+vec(icount)*B(icount)
	end do
	BB_temp(i)=BB_temp(i)+I_fluid_center(ie)*temp*Sp/(hs**nsd)*(hsg**nsd)
     end do
  end do
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(A_temp(1,1),A(1,1),nn_inter*nn_inter,mpi_double_precision, &
                mpi_sum,mpi_comm_world,ierror)
  call mpi_allreduce(BB_temp(1),BB(1),nn_inter,mpi_double_precision, &
                mpi_sum,mpi_comm_world,ierror)

  call mpi_barrier(mpi_comm_world,ierror)

     BB(:)=0.5-BB(:)
     call gmres_correction(A,BB,nn_inter)
  !   call DGESV(nn_inter,1,A,nn_inter,IPIV,BB,nn_inter,INFO)

     corr_Ip(1:nn_inter)=BB(1:nn_inter) !calculate delta Ip
end subroutine get_correction_nu
