

subroutine gmres_new(x,w,bg,dg,ien,id,jac, &
			ne_local,ien_local,node_local,nn_local, &
			global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length)
	use fluid_variables, only: nsd,nn,ne,nen,ndf,inner,outer,nquad
        use mpi_variables
	implicit none
     	 include 'mpif.h'
! use GMRES to solve mesh update equation  with diagonal preconditioner
	real* 8 x(nsd,nn),id(nsd,nn)
	real* 8 ien(nen,ne)
	real* 8 bg(nsd*nn), dg(nsd*nn), w(nsd*nn)
	real* 8 jac(nquad,ne)
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
!=========================================================
! reduce dimension to save memory
        real* 8 Vm(nsd*nn_local, inner+1) ! Krylov space matrix
!        real* 8 Vm(nsd*nn, inner+1) ! Krylov space matrix
!=========================================================


	integer i,j,iouter,icount,INFO
	integer e1(inner+1)
	real* 8 x0(nsd*nn)
	real* 8 beta(inner+1)
	real* 8 eps
	real* 8 r0(nsd*nn)
	real* 8 rnorm, rnorm0,err
        real* 8 dv(nsd*nn)
	real* 8 Vy(nsd*nn)
        real* 8 vloc(nsd,nn), avloc(nsd*nn)
!==============================
! reduce dimension save memory
!       real* 8 temp(ndf*nn)
        real* 8 temp
!==============================

	character(1) TRAN
	real* 8 workls(2*inner)
	real* 8 av_tmp(nsd,nn)

	integer time_arrary_0(8)
        integer time_arrary_1(8)
	real(8) start_time
	real(8) end_time
	real(8) err_givens
!============================
! MPI varibalbes
  integer ne_local ! # of element on each processor
  integer ien_local(ne_local) ! subregion-->wholeregion element index
  integer ie_local ! loop parameter
  integer node_local(nn_local)
  integer nn_local
  integer jcount
  integer kcount
  integer node
  integer nn_global_com
  integer global_com(nn_global_com)  ! global node index for communication
  integer nn_local_com
  integer local_com(nn_local_com)  ! local index in the communication region on each processor
  integer ad_length
  integer send_address(ad_length,2)
  integer flag
  real(8) dg_sent(nsd*nn)
  real(8) space1(inner)
  real(8) space2(inner)
!---------------------------------------------
	eps = 1.0e-6
	e1(:) = 0
	e1(1) = 1
	x0(:) = 0
	iouter = 1
	r0(:) = bg(:)
	TRAN = 'N'
	av_tmp(:,:) = 0
	avloc(:) = 0
	dv(:) = 0
!	w(:) = 1.0
        call getnorm_pa(r0,nsd,nn,node_local,nn_local,rnorm0)
        rnorm = sqrt(rnorm0)

!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. eps))

!====================================
	Vm(:,:) = 0.0d0
! Already is a local variable defined on each proc
!===============================

        do icount = 1, nn_local
           node=node_local(icount)
           do jcount=1,nsd
                Vm((icount-1)*nsd+jcount,1) = r0((node-1)*nsd+jcount)/rnorm
           end do
        end do ! get V1

	   beta(:) = rnorm*e1(:) ! get beta*e1
	   Hm(:,:) = 0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
!        call date_and_time(values=time_arrary_0)  
	   do j=1,inner
	  

                 do icount=1, nn_local
                    node=node_local(icount)
                    do jcount=1,nsd
                        dv((node-1)*nsd+jcount) = 1/w((node-1)*nsd+jcount)*Vm((icount-1)*nsd+jcount,j)
                    end do
                 end do!!!!!!!!!!calcule eps*inv(P)*V1
!===============================================
!		vloc(:,:) = 0.0d0
! Clear the matrix locally to aviod too many loops 
                do icount=1,nn_local
                   node=node_local(icount)
                   vloc(1:nsd,node)=0.0d0
                end do
                ! clear all the processor internal nodes
                do icount=1,nn_local_com
                   node=global_com(local_com(icount))
                   vloc(1:nsd,node)=0.0d0
                end do
                ! clear the processor boundary nodes
!=========================================


                call equal_pa(dv,vloc,nsd,nn,node_local,nn_local)
!                call communicate_res(global_com,nn_global_com,local_com,nn_local_com,vloc,nsd,nn)
                call communicate_res_ad(vloc,nsd,nn,send_address,ad_length)
!		vloc(:,:) = vloc(:,:)+d(:,:)  ! calculate u+eps*inv(P)*V1
!===============================================
!		av_tmp(:,:) = 0.0d0
! Clear matrix local first interal then boundary processor nodes
                do icount=1,nn_local
                   node=node_local(icount)
                   av_tmp(1:nsd,node)=0.0d0
                end do
                do icount=1,nn_local_com
                   node=global_com(local_com(icount))
                   av_tmp(1:nsd,node)=0.0d0
                end do
!================================================


        call date_and_time(values=time_arrary_0)

		call blockgmresm(x,vloc,av_tmp,ien,jac,ne_local,ien_local)

        call date_and_time(values=time_arrary_1)
        start_time=time_arrary_0(5)*3600+time_arrary_0(6)*60+time_arrary_0(7)+time_arrary_0(8)*0.001
        end_time=time_arrary_1(5)*3600+time_arrary_1(6)*60+time_arrary_1(7)+time_arrary_1(8)*0.001
!        if (myid == 0) write(*,*) 'Block GMRESN time', end_time - start_time




!                call communicate_res(global_com,nn_global_com,local_com,nn_local_com,av_tmp,nsd,nn)
                call communicate_res_ad(av_tmp,nsd,nn,send_address,ad_length)

!		avloc(:)=0.0d0
                call equal_pa(av_tmp,avloc,nsd,nn,node_local,nn_local)
!		avloc(:) = (-avloc(:)+bg(:))/eps ! get Av,bg=-r(u)
!		call setid(avloc,id,nsd)
                call setid_pa(avloc,nsd,nn,id,node_local,nn_local)

        call date_and_time(values=time_arrary_0)
		space1(:)=0.0d0
		space2(:)=0.0d0
	      do i=1,j

                 do icount=1,nn_local
                        node=node_local(icount)
                        do jcount=1,nsd
                                kcount=(node-1)*nsd+jcount
                                space1(i)=space1(i)+avloc(kcount)*Vm((icount-1)*nsd+jcount,i)
                        end do
                 end do
!                call vector_dot_pa(avloc,Vm(:,i),nsd,nn,nn_local,node_local,space1(i))
	      end do
		call mpi_barrier(mpi_comm_world,ierror)
		call mpi_allreduce(space1(1),space2(1),j,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
		Hm(1:j,j)=space2(1:j)
!		if (myid == 0) write(*,*) Hm(:,:)
        call date_and_time(values=time_arrary_1)
        start_time=time_arrary_0(5)*3600+time_arrary_0(6)*60+time_arrary_0(7)+time_arrary_0(8)*0.001
        end_time=time_arrary_1(5)*3600+time_arrary_1(6)*60+time_arrary_1(7)+time_arrary_1(8)*0.001
!        if (myid == 0) write(*,*) ' In GMRESN Comunicate time', end_time - start_time

              do icount = 1, nn_local
                node=node_local(icount)
                do jcount=1,nsd
                        kcount=(icount-1)*nsd+jcount
                
                        do i=1,j
                                Vm(kcount,j+1) = Vm(kcount,j+1)-Hm(i,j)*Vm(kcount,i)
                        end do
                        Vm(kcount,j+1)=Vm(kcount,j+1)+avloc((node-1)*nsd+jcount)
                end do
              end do  ! construct v(j+1)
temp=0.0d0
              do icount = 1, nn_local
                do jcount=1,nsd
                        kcount=(icount-1)*nsd+jcount
                        temp=temp+Vm(kcount,j+1)*Vm(kcount,j+1)
                end do
              end do
                call mpi_barrier(mpi_comm_world,ierror)
                call mpi_allreduce(temp,rnorm0,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
!              call getnorm_pa(temp,nsd,nn,node_local,nn_local,rnorm0)

              Hm(j+1,j) = sqrt(rnorm0)

              do icount = 1, nn_local
!                node=node_local(icount)
                do jcount=1,nsd
                        kcount=(icount-1)*nsd+jcount
                        Vm(kcount,j+1)=Vm(kcount,j+1)/Hm(j+1,j)
                end do
              end do
	  end do  ! end inner loop

!	call date_and_time(values=time_arrary_1)
!	start_time=time_arrary_0(5)*3600+time_arrary_0(6)*60+time_arrary_0(7)+time_arrary_0(8)*0.001
!        end_time=time_arrary_1(5)*3600+time_arrary_1(6)*60+time_arrary_1(7)+time_arrary_1(8)*0.001
!	if (myid == 0) write(*,*) 'Inner loop time', end_time - start_time


!	call date_and_time(values=time_arrary_0)
!===================================================================================
! Use Givens rotation to solve the LS problem
           call mpi_barrier(mpi_comm_world,ierror)
	call givens(Hm,inner,beta)
!====================================================================================           
! Use lapack to solve the LS problem
!	call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
	
!        call date_and_time(values=time_arrary_1)
!        start_time=time_arrary_0(5)*3600+time_arrary_0(6)*60+time_arrary_0(7)+time_arrary_0(8)*0.001
!        end_time=time_arrary_1(5)*3600+time_arrary_1(6)*60+time_arrary_1(7)+time_arrary_1(8)*0.001
!	write(*,*) 'LS time by lapack', end_time-start_time
!!!!!!!!!!!!beta(1:inner) is ym, the solution!!!!!!!!!!!!!!!!!!!!!!!!!
	Vy(:) = 0
        do icount=1, nn_local
           node=node_local(icount)
           do jcount=1,nsd
                kcount=(node-1)*nsd+jcount
                do i=1,inner
                         Vy(kcount)=Vy(kcount)+Vm((icount-1)*nsd+jcount,i)*beta(i)
                end do
                x0(kcount)=x0(kcount)+Vy(kcount)
           end do
        end do ! calculate Xm

        do icount = 1, nn_local
           node=node_local(icount)
           do jcount=1,nsd
                kcount=(node-1)*nsd+jcount
                dv(kcount) = 1/w(kcount)*x0(kcount)
           end do
        end do

!============================
! Clear matirx locally 
                do icount=1,nn_local
                   node=node_local(icount)
                   vloc(1:nsd,node)=0.0d0
                end do
                ! clear all the processor internal nodes
                do icount=1,nn_local_com
                   node=global_com(local_com(icount))
                   vloc(1:nsd,node)=0.0d0
                end do
                ! clear the processor boundary nodes
!	vloc(:,:) = 0
!==========================
        call equal_pa(dv,vloc,nsd,nn,node_local,nn_local)
!        call communicate_res(global_com,nn_global_com,local_com,nn_local_com,vloc,nsd,nn)
        call communicate_res_ad(vloc,nsd,nn,send_address,ad_length)
!	vloc(:,:) = vloc(:,:)+d(:,:)
!===============================
! Clear matrix local first interal then boundary processor nodes
                do icount=1,nn_local
                   node=node_local(icount)
                   av_tmp(1:nsd,node)=0.0d0
                end do
                do icount=1,nn_local_com
                   node=global_com(local_com(icount))
                   av_tmp(1:nsd,node)=0.0d0
                end do
!	av_tmp(:,:) = 0.0d0
!=============================
        call blockgmresm(x,vloc,av_tmp,ien,jac,ne_local,ien_local)
!        call communicate_res(global_com,nn_global_com,local_com,nn_local_com,av_tmp,nsd,nn)
         call communicate_res_ad(av_tmp,nsd,nn,send_address,ad_length)

!avloc(:)=0.0d0
       call equal_pa(av_tmp,avloc,nsd,nn,node_local,nn_local)
!        call setid(avloc,id,nsd)
         call setid_pa(avloc,nsd,nn,id,node_local,nn_local)


!	avloc(:) = (-avloc(:)+bg(:))/eps
!!!!!!!!!!calculate AXm

        do icount=1,nn_local
           node=node_local(icount)
           do jcount=1,nsd
              kcount=(node-1)*nsd+jcount
              r0(kcount) = bg(kcount)-avloc(kcount)
           end do
        end do !update r0=f-AX0

        call getnorm_pa(r0,nsd,nn,node_local,nn_local,rnorm0)
        err = sqrt(rnorm0)


	rnorm = sqrt(rnorm0)
	iouter = iouter + 1        

	if (myid == 0) write(*,*) 'Mesh GMRES ERROR =',err
111     continue  ! end outer loop

        dg_sent(:)=0.0d0
        do icount=1, nn_local
           node=node_local(icount)
           do jcount=1,nsd
              kcount=(node-1)*nsd+jcount
              dg_sent(kcount) = 1/w(kcount)*x0(kcount) ! unscaled x0
           end do
        end do
        dg(:)=0.0d0
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_allreduce(dg_sent(1),dg(1),nsd*nn,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
!        call mpi_bcast(dg(1),nsd*nn,mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_barrier(mpi_comm_world,ierror)

	return
	end
