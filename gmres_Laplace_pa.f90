

subroutine gmres_Laplace_pa(x,d,w,bg,dg,ien,id,&
		ne_local,ien_local,nn_local,node_local,send_address,ad_length,&
		global_com,nn_global_com,local_com,nn_local_com)
	use fluid_variables, only: nsd,nn,ne,nen
	use mpi_variables
	implicit none
	include 'mpif.h'
	integer inner, outer
	parameter (inner = 100) ! Laplace equation inner 50 should be sufficient
	parameter (outer = 10)  ! Laplace equation outer 5 should be sufficient
	integer ien(nen,ne)
	real* 8 x(nsd,nn)
	real* 8 d(nn)
	real* 8 bg(nn), dg(nn), w(nn)
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
	real* 8 Vm(nn_local, inner+1) ! Krylov space matrix
!--------------------------------
	integer id(nn)
!--------------------------------
	integer i,j,iouter,icount,INFO
	real* 8  e1(inner+1)
	real* 8 x0(nn)
	real* 8 beta(inner+1)
	real* 8 eps
	real* 8 r0(nn)
	real* 8 rnorm, rnorm0,err
        real* 8 dv(nn)
	real* 8 Vy(nn)
        real* 8 avloc(nn)
	real* 8 temp
	character(1) TRAN
	real* 8 workls(2*inner)
!-------------------------------
! MPI varibalbes
  integer ne_local ! # of element on each processor
  integer ien_local(ne_local) ! subregion-->wholeregion element index
  integer ie_local ! loop parameter
  integer node_local(nn_local)
  integer nn_local
  integer node
  integer nn_global_com
  integer global_com(nn_global_com)  ! global node index for communication
  integer nn_local_com
  integer local_com(nn_local_com)  ! local index in the communication region on each processor
  integer ad_length
  integer send_address(ad_length,2)
  real(8) dg_sent(nn)
  real(8) space1(inner)
  real(8) space2(inner)

	eps = 1.0e-6
	e1(:) = 0.0
	e1(1) = 1
	x0(:) = 0.0
	iouter = 1
	r0(:) = bg(:)
	TRAN = 'N'
	avloc(:) = 0.0
!	w(:) = 1
        call getnorm_pa(r0,1,nn,node_local,nn_local,rnorm0)
        rnorm = sqrt(rnorm0)
if (myid == 0) write(*,*) 'rnorm0', rnorm0
!do icount=1,nn
!	write(*,*) 'node id',icount, id(icount)
!end do


!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. eps))

	Vm(:,:) = 0.0
	do icount = 1, nn_local
	node=node_local(icount)
	   Vm(icount,1) = r0(node)/rnorm
	end do ! get V1

	   beta(:) = rnorm*e1(:) ! get beta*e1
	   Hm(:,:) = 0.0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
	   do j=1,inner
	  

		dv(:) = 0.0d0

		 do icount=1, nn_local
		 node=node_local(icount)
		    dv(node) = 1.0/w(node)*Vm(icount,j)
		 end do!!!!!!!!!!calcule eps*inv(P)*V1
		

                call communicate_res_ad_sub(dv,1,nn,send_address,ad_length)





		avloc(:) = 0.0d0
		call blockgmres_Laplace(x,dv,avloc,ien,ien_local,ne_local)
!		avloc(:) = (-avloc(:)+bg(:))/eps ! get Av,bg=-r(u)
                call communicate_res_ad_sub(avloc,1,nn,send_address,ad_length)
		call setid_pa(avloc,1,nn,id,node_local,nn_local)
		continue
	 

	space1(:)=0.0d0
	space2(:)=0.0d0



		do i=1,j
		 do icount = 1, nn_local
		   node=node_local(icount)
		   space1(i) = space1(i)+avloc(node)*Vm(icount,i)
		 end do
		end do

		call mpi_barrier(mpi_comm_world,ierror)
		call mpi_allreduce(space1(1),space2(1),j,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
		Hm(1:j,j)=space2(1:j)
	         
	      do icount = 1, nn_local
	      	node=node_local(icount)
		 do i=1,j
		   Vm(icount,j+1) = Vm(icount,j+1)-Hm(i,j)*Vm(icount,i)
		 end do
		 Vm(icount,j+1)=Vm(icount,j+1)+avloc(node)
	      end do  ! construct v(j+1)

temp=0.0d0
	      do icount = 1, nn_local
		 temp= temp + Vm(icount,j+1)*Vm(icount,j+1)
	      end do

                call mpi_barrier(mpi_comm_world,ierror)
		call mpi_allreduce(temp,rnorm0,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)


	      Hm(j+1,j) = sqrt(rnorm0)

	      do icount = 1, nn_local
 		 Vm(icount,j+1)=Vm(icount,j+1)/Hm(j+1,j)
	      end do
	  end do  ! end inner loop
		
!==========================================================================
!	call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
	call mpi_barrier(mpi_comm_world,ierror)


	call givens(Hm,inner,beta)
!!!!!!!!!!!!beta(1:inner) is ym, the solution!!!!!!!!!!!!!!!!!!!!!!!!!
	Vy(:) = 0
	do icount=1,nn_local
	node=node_local(icount)
	   do i=1,inner
	      Vy(node)=Vy(node)+Vm(icount,i)*beta(i)
	   end do
	   x0(node)=x0(node)+Vy(node)
	end do ! calculate Xm
!write(*,*)'x0=',x0(:)

!        do icount=1,nn_local
!	node=global_com(local_com(icount))
!	dv(node)=0.0d0
!	end do

	dv(:)=0.0d0

	do icount = 1, nn_local
	   node = node_local(icount)
	   dv(node) = 1.0/w(node)*x0(node)
	end do


        call communicate_res_ad_sub(dv,1,nn,send_address,ad_length)

	
	avloc(:) = 0.0d0
	call blockgmres_Laplace(x,dv,avloc,ien,ien_local,ne_local)
!	avloc(:) = (-avloc(:)+bg(:))/eps
        call communicate_res_ad_sub(avloc,1,nn,send_address,ad_length)
        call setid_pa(avloc,1,nn,id,node_local,nn_local)


!!!!!!!!!!calculate AXm

	do icount=1,nn_local
	node=node_local(icount)
	   r0(node) = bg(node)-avloc(node)
	end do !update r0=f-AX0
	call getnorm_pa(r0,1,nn,node_local,nn_local,rnorm0)
	err = sqrt(rnorm0)

	rnorm = sqrt(rnorm0)
	iouter = iouter + 1        

if (myid == 0)	write(*,*) '***** Indicator Laplace Eq err******=',err
111     continue  ! end outer loop
!write(*,*)'x0_2=',x0(:)

	dg_sent(:)=0.0d0
	do icount=1,nn_local
		node=node_local(icount)
		dg_sent(node)=1.0/w(node)*x0(node)
	end do
	dg(:)=0.0d0
	call mpi_barrier(mpi_comm_world,ierror)
	call mpi_allreduce(dg_sent(1),dg(1),nn,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
	call mpi_barrier(mpi_comm_world,ierror)
	return
	end
