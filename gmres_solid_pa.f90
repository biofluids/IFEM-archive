

subroutine gmres_solid_pa(x,w,bg,dg,ien,id,nsd,nn,ne,nen,inner,outer,nquad,wq,sq,x_pre1,&
                        solid_prevel,solid_preacc,solid_stress,ne_sbc,ien_sbc,mtype,&
			ne_local,ien_local,node_local,nn_local,&
			global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length,ndf)
        use mpi_variables
	implicit none
      include 'mpif.h'
!==============================
! reduce dimension save memory
!	real* 8 temp(ndf*nn)
! use GMRES to solve mesh update equation  with diagonal preconditioner
        real* 8 x(nsd,nn),id(nsd,nn)
	integer ien(ne,nen)
	real* 8 bg(nsd*nn), dg(nsd*nn), w(nsd*nn)
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
	real* 8 Vm(nsd*nn_local, inner+1) ! Krylov space matrix
	real(8) x_pre1(nsd,nn)
	real(8) solid_prevel(nsd,nn)
	real(8) solid_preacc(nsd,nn)
	!       real(8) x_pre2(nsd,nn)
	real(8) solid_stress(nsd*2,nn) ! fluid stress acting on the solid boundary including pressure and viscous stress
	integer ne_sbc                  ! solid elements on the boundary
	integer ien_sbc(ne_sbc,nen+2)
	integer mtype(ne)
	!-----------------------------------------------
	integer nsd
	integer ndf
	integer nn
	integer ne
	integer nen
	integer inner
	integer outer
	integer nquad
	real(8) wq(8)
	real(8) sq(0:3,8,8)
	!---------------------------------------------
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
	real* 8 temp
!==============================
	character(1) TRAN
	real* 8 workls(2*inner)
	real* 8 av_tmp(ndf,nn)

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
  real(8) dg_sent(ndf*nn)
  real(8) space1(inner)
  real(8) space2(inner)
        integer time_arrary_0(8)
        integer time_arrary_1(8)
        real(8) start_time
        real(8) end_time

!---------------------------------------------



	eps = 1.0d-6
	e1(:) = 0.0
	e1(1) = 1.0
	x0(:) = 0.0
	iouter = 1
	r0(:) = bg(:)
	TRAN = 'N'
	av_tmp(:,:) = 0.0
	avloc(:) = 0.0
	vloc(:,:)=0.0
!	w(:) = 1
        call getnorm_pa(r0,ndf,nn,node_local,nn_local,rnorm0)
        rnorm = sqrt(rnorm0)


!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. eps))

!==============================
	Vm(:,:) = 0.0d0
! Already is a local variable defined on each proc
!===============================
	do icount = 1, nn_local
	   node=node_local(icount)
	   do jcount=1,ndf
	   	Vm((icount-1)*ndf+jcount,1) = r0((node-1)*ndf+jcount)/rnorm
	   end do
	end do ! get V1

	   beta(:) = rnorm*e1(:) ! get beta*e1
	   Hm(:,:) = 0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
	   do j=1,inner
	  

		 do icount=1, nn_local
		    node=node_local(icount)
		    do jcount=1,ndf
		    	dv((node-1)*ndf+jcount) = 1.0/w((node-1)*ndf+jcount)*Vm((icount-1)*ndf+jcount,j)
		    end do
		 end do!!!!!!!!!!calcule eps*inv(P)*V1
		
!===============================================
! Clear the matrix locally to aviod too many loops 
		do icount=1,nn_local
		   node=node_local(icount)
		   vloc(1:ndf,node)=0.0d0
		end do
		! clear all the processor internal nodes
		do icount=1,nn_local_com
		   node=global_com(local_com(icount))
		   vloc(1:ndf,node)=0.0d0
		end do
		! clear the processor boundary nodes

!===============================================
		call equal_pa(dv,vloc,ndf,nn,node_local,nn_local)

!============================
	        call communicate_res_ad_subsolid(vloc,ndf,nn,send_address,ad_length)

!----------------------------------------------------------------------------------------------
!=======================================
! Clear matrix local first interal then boundary processor nodes
		do icount=1,nn_local
                   node=node_local(icount)
                   av_tmp(1:ndf,node)=0.0d0
                end do
                do icount=1,nn_local_com
                   node=global_com(local_com(icount))
                   av_tmp(1:ndf,node)=0.0d0
                end do
!======================================

                call blockgmres_solid_pa(x,vloc,av_tmp,ien,nsd,nen,ne,nn,nquad,wq,sq,x_pre1,&
		                solid_prevel,solid_preacc,ien_sbc,ne_sbc,solid_stress,mtype,&
				ne_local,ien_local)




	        call communicate_res_ad_subsolid(av_tmp,ndf,nn,send_address,ad_length)

		call equal_pa(av_tmp,avloc,ndf,nn,node_local,nn_local)

		!call setid(avloc,id,ndf)
                call setid_pa(avloc,ndf,nn,id,node_local,nn_local)
                space1(:)=0.0d0
                space2(:)=0.0d0

		do i=1,j
		 do icount=1,nn_local
			node=node_local(icount)
			do jcount=1,ndf
				kcount=(node-1)*ndf+jcount
				space1(i)=space1(i)+avloc(kcount)*Vm((icount-1)*ndf+jcount,i)
			end do
		 end do
		! Do the vector product of v_i * v_i set it to be h_i,j
	      end do  ! construct AVj and hi,j

                call mpi_barrier(mpi_comm_world,ierror)
                call mpi_allreduce(space1(1),space2(1),j,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

                Hm(1:j,j)=space2(1:j)

	      do icount = 1, nn_local
		node=node_local(icount)
		do jcount=1,ndf
			kcount=(icount-1)*ndf+jcount
		
		 	do i=1,j
		   		Vm(kcount,j+1) = Vm(kcount,j+1)-Hm(i,j)*Vm(kcount,i)
		 	end do
		 	Vm(kcount,j+1)=Vm(kcount,j+1)+avloc((node-1)*ndf+jcount)
	        end do  
	      end do  ! construct v(j+1)
!
temp=0.0d0
	      do icount = 1, nn_local
!		node=node_local(icount)
		do jcount=1,ndf
			kcount=(icount-1)*ndf+jcount
		 	temp=temp+Vm(kcount,j+1)*Vm(kcount,j+1)
		end do
	      end do
		call mpi_barrier(mpi_comm_world,ierror)
		call mpi_allreduce(temp,rnorm0,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

	      Hm(j+1,j) = sqrt(rnorm0)

	      do icount = 1, nn_local
		do jcount=1,ndf
			kcount=(icount-1)*ndf+jcount
 			Vm(kcount,j+1)=Vm(kcount,j+1)/Hm(j+1,j)
		end do
	      end do


	  end do  ! end inner loop

!if (myid == 1) write(*,*) 'Hm', Hm(:,:)
!====================================================================================		
! Use lapack to solve the LS problem
!	call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
!===================================================================================
! Use Givens rotation to solve the LS problem
           call mpi_barrier(mpi_comm_world,ierror)

call givens(Hm,inner,beta)

!!!!!!!!!!!!beta(1:inner) is ym, the solution!!!!!!!!!!!!!!!!!!!!!!!!!
	Vy(:) = 0
	do icount=1, nn_local
	   node=node_local(icount)
	   do jcount=1,ndf
		kcount=(node-1)*ndf+jcount
	   	do i=1,inner
	       		 Vy(kcount)=Vy(kcount)+Vm((icount-1)*ndf+jcount,i)*beta(i)
	        end do
	        x0(kcount)=x0(kcount)+Vy(kcount)
	   end do
	end do ! calculate Xm

	do icount = 1, nn_local
	   node=node_local(icount)
	   do jcount=1,ndf
		kcount=(node-1)*ndf+jcount
	   	dv(kcount) = 1.0/w(kcount)*x0(kcount)
	   end do
	end do

!write(*,*) dv(:)
!============================
! Clear matirx locally 
                do icount=1,nn_local
                   node=node_local(icount)
                   vloc(1:ndf,node)=0.0d0
                end do
                ! clear all the processor internal nodes
                do icount=1,nn_local_com
                   node=global_com(local_com(icount))
                   vloc(1:ndf,node)=0.0d0
                end do
                ! clear the processor boundary nodes
!==============================


	call equal_pa(dv,vloc,ndf,nn,node_local,nn_local)

        call communicate_res_ad_subsolid(vloc,ndf,nn,send_address,ad_length)


! Clear matrix local first interal then boundary processor nodes
                do icount=1,nn_local
                   node=node_local(icount)
                   av_tmp(1:ndf,node)=0.0d0
                end do
                do icount=1,nn_local_com
                   node=global_com(local_com(icount))
                   av_tmp(1:ndf,node)=0.0d0
                end do
!	av_tmp(:,:) = 0.0d0
!================================

        call blockgmres_solid_pa(x,vloc,av_tmp,ien,nsd,nen,ne,nn,nquad,wq,sq,x_pre1,&
	                solid_prevel,solid_preacc,ien_sbc,ne_sbc,solid_stress,mtype,&
			ne_local,ien_local)


        call communicate_res_ad_subsolid(av_tmp,ndf,nn,send_address,ad_length)
       call equal_pa(av_tmp,avloc,ndf,nn,node_local,nn_local)
       call setid_pa(avloc,ndf,nn,id,node_local,nn_local)

!!!!!!!!!!calculate AXm
	do icount=1,nn_local
	   node=node_local(icount)
           do jcount=1,ndf
              kcount=(node-1)*ndf+jcount
	      r0(kcount) = bg(kcount)-avloc(kcount)
	   end do
	end do !update r0=f-AX0

	call getnorm_pa(r0,ndf,nn,node_local,nn_local,rnorm0)
	err = sqrt(rnorm0)

	rnorm = sqrt(rnorm0)
	iouter = iouter + 1        

	if (myid == 0)	write(*,*) 'err=',err


111     continue  ! end outer loop


	dg_sent(:)=0.0d0
	do icount=1, nn_local
           node=node_local(icount)
           do jcount=1,ndf
              kcount=(node-1)*ndf+jcount
	      dg_sent(kcount) = 1/w(kcount)*x0(kcount) ! unscaled x0
	   end do
	end do
	dg(:)=0.0d0
	call mpi_barrier(mpi_comm_world,ierror)
	call mpi_allreduce(dg_sent(1),dg(1),ndf*nn,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        call mpi_barrier(mpi_comm_world,ierror)


	return
	end
