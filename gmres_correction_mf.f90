
subroutine gmres_correction_mf(x_inter,B,w,RW,nsd,nn_inter,nn_inter_loc)

  use fluid_variables, only:inner,outer
  use interface_variables, only:maxmatrix
  use mpi_variables
  implicit none
  include 'mpif.h'

  real(8) x_inter(nsd,maxmatrix)
  integer nn_inter,nsd
  real(8) B(nn_inter)
  real(8) w(nn_inter)
  real(8) RW(nsd+1,nn_inter)
  integer loc_index,nn_inter_loc

  real(8) Hm(inner+1,inner) ! henssenberg matrix
  real(8) Vm(nn_inter,inner+1) ! Krylov space matrix

  integer i,j,iouter,icount,jcount,INFO
  real(8) x0(nn_inter)
  real(8) beta(inner+1)
  real(8) eps, r0(nn_inter)
  real(8) rnorm,rnorm0,err
  real(8) dv(nn_inter)
  real(8) Vy(nn_inter)
  real(8) avloc(nn_inter)
  real(8) temp(nn_inter)
  character(1) TRAN
  real(8) workls(2*inner)
  real(8) tmp,tmp1
  real(8) B_temp(nn_inter)
  real(8) space1(nn_inter),space2(nn_inter)

  eps=1.0e-6
  x0(:)=0.0
  iouter=1
  r0(:)=B(:)
  B_temp(:)=B(:)
  TRAN = 'N'

!  call getnorm(r0,r0,nn_inter,rnorm0)
  rnorm0=0.0
  tmp=0.0
  do loc_index=1,nn_inter_loc
     i=myid+1+(loc_index-1)*ncpus
     tmp=tmp+r0(i)*r0(i)
  end do
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(tmp,rnorm0,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
  rnorm=sqrt(rnorm0)
if(myid==0)write(*,*)'begin gmres correction,rnorm=',rnorm
!+++++++++++++start outer loop+++++++++++++++++

  do 111,while((iouter.le.outer).and.(rnorm.ge.1.0e-16))

     Vm(:,:)=0.0
!     do icount=1,nn_inter
     do loc_index=1,nn_inter_loc
	icount=myid+1+(loc_index-1)*ncpus
        Vm(icount,1)=r0(icount)/rnorm
     end do !get V1

     beta(:)=0.0
     beta(1)=rnorm
     Hm(:,:)=0.0

!+++++++++++start inner loop +++++++++++++++++++
     do j=1,inner

        dv(:)=0.0
	temp(:)=0.0

!	do icount=1,nn_inter
	do loc_index=1,nn_inter_loc
	   icount=myid+1+(loc_index-1)*ncpus
	   temp(icount)=1/w(icount)*Vm(icount,j)
	end do

	call mpi_barrier(mpi_comm_world,ierror)
	call mpi_allreduce(temp,dv,nn_inter,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

	avloc(:)=0.0
	call blockgmres_correction(x_inter,RW,nn_inter,dv,avloc,nn_inter_loc,nsd)
!	do icount=1,nn_inter
!	do loc_index=1,nn_inter_loc
!	   icount=myid+1+(loc_index-1)*ncpus
!	   do jcount=1,nn_count(icount)
!	      avloc(icount)=avloc(icount)+A(icount,jcount)*dv(A_index(icount,jcount))
!	   end do
!	end do
	space1(:)=0.0
	space2(:)=0.0

	do i=1,j

!	   do icount=1,nn_inter
	   do loc_index=1,nn_inter_loc
	      icount=myid+1+(loc_index-1)*ncpus
!	      Hm(i,j)=Hm(i,j)+avloc(icount)*Vm(icount,i)
!	      tmp=tmp+avloc(icount)*Vm(icount,i)
	      space1(i)=space1(i)+avloc(icount)*Vm(icount,i)
	   end do
!	   call mpi_barrier(mpi_comm_world,ierror)
!	   call mpi_allreduce(tmp,tmp1,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
!	   Hm(i,j)=tmp1
	end do  ! construct Avj and hi,j

	call mpi_barrier(mpi_comm_world,ierror)
	call mpi_allreduce(space1(1),space2(1),j,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
	Hm(1:j,j)=space2(1:j)


!	do icount=1,nn_inter
	do loc_index=1,nn_inter_loc
	   icount=myid+1+(loc_index-1)*ncpus
	   do i=1,j
	      Vm(icount,j+1)=Vm(icount,j+1)-Hm(i,j)*Vm(icount,i)
	   end do
	   Vm(icount,j+1)=Vm(icount,j+1)+avloc(icount)
	end do   !construct v(j+1)
	   
  
!	do icount=1,nn_inter
!	   temp(icount)=Vm(icount,j+1)
!	end do
!	call getnorm(temp,temp,nn_inter,rnorm0) 
	tmp=0.0
	do loc_index=1,nn_inter_loc
	   icount=myid+1+(loc_index-1)*ncpus
	   tmp=tmp+Vm(icount,j+1)*Vm(icount,j+1)
	end do
        rnorm0=0.0
	call mpi_barrier(mpi_comm_world,ierror)
	call mpi_allreduce(tmp,rnorm0,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)


	Hm(j+1,j)=sqrt(rnorm0)

!	do icount=1,nn_inter
	do loc_index=1,nn_inter_loc
	   icount=myid+1+(loc_index-1)*ncpus
	   Vm(icount,j+1)=Vm(icount,j+1)/Hm(j+1,j)
	end do


    end do
! end innner loop

!    call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
    call givens(Hm,inner,beta)
    Vy(:)=0.0
!    do icount=1,nn_inter
    do loc_index=1,nn_inter_loc
	icount=myid+1+(loc_index-1)*ncpus
	do i=1,inner
	   Vy(icount)=Vy(icount)+Vm(icount,i)*beta(i)
	end do
	x0(icount)=x0(icount)+Vy(icount)
    end do ! calculate Xm

!    do icount=1,nn_inter
    temp(:)=0.0
    dv(:)=0.0
    do loc_index=1,nn_inter_loc
	icount=myid+1+(loc_index-1)*ncpus
	temp(icount)=1.0/w(icount)*x0(icount)
    end do

    call mpi_barrier(mpi_comm_world,ierror)
    call mpi_allreduce(temp,dv,nn_inter,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

    avloc(:)=0.0
    call blockgmres_correction(x_inter,RW,nn_inter,dv,avloc,nn_inter_loc,nsd)
!    do icount=1,nn_inter
!    do loc_index=1,nn_inter_loc
!	icount=myid+1+(loc_index-1)*ncpus
!	do jcount=1,nn_count(icount)
!	   avloc(icount)=avloc(icount)+A(icount,jcount)*dv(A_index(icount,jcount))
!	end do
!    end do !calculate Axm

    tmp=0.0
!    do icount=1,nn_inter
    do loc_index=1,nn_inter_loc
	icount=myid+1+(loc_index-1)*ncpus
	r0(icount)=B(icount)-avloc(icount)
	tmp=tmp+r0(icount)*r0(icount)
    end do !update r0=f-Ax0

    rnorm0=0.0
    call mpi_barrier(mpi_comm_world,ierror)
    call mpi_allreduce(tmp,rnorm0,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

!    call getnorm(r0,r0,nn_inter,rnorm0)
    err=sqrt(rnorm0)
    rnorm=sqrt(rnorm0)
    iouter=iouter+1
if(myid==0) then
	write(*,*)'err for correction=',err
end if

111 continue ! end outer loop

!    B(1:nn_inter)=1.0/w(1:nn_inter)*x0(1:nn_inter)
   temp(:)=0.0
   B(:)=0.0
   do loc_index=1,nn_inter_loc
	icount=myid+1+(loc_index-1)*ncpus
	temp(icount)=1.0/w(icount)*x0(icount)
   end do
   call mpi_barrier(mpi_comm_world,ierror)
   call mpi_allreduce(temp(1),B(1),nn_inter,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

end subroutine gmres_correction_mf

    





















