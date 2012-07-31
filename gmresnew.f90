

subroutine gmres(x,d,dold,w,bg,dg,hg,ien,fext,id, &
		ne_local,ien_local,node_local,nn_local, &
		global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length,&
		sur_fluid,I_fluid,&
		norm_node,spbcnode,spbcele,rngface,pbnode)
	use fluid_variables, only: nsd,nn,ne,nen,ndf,inner,outer,neface,nn_spbc,ne_spbc,nn_pb
 	use solid_variables, only: nn_solid
        use mpi_variables
	implicit none
      include 'mpif.h'
	real* 8 x(nsd,nn),id(ndf,nn)
	real* 8 d(ndf,nn), dold(ndf,nn),hg(ne),fext(ndf,nn)!,ien(nen,ne)
	integer ien(nen,ne)
	real* 8 bg(ndf*nn), dg(ndf*nn), w(ndf*nn)
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
	real* 8 sur_fluid(nsd,nn),I_fluid(nn)
        real* 8 norm_node(nsd,nn)
        integer spbcnode(nn_spbc)
        integer spbcele(ne_spbc)
        integer rngface(neface,ne)
	integer pbnode(2,nn_pb)
	real(8) res_pb(nsd,nn_pb),res_pb_temp(nsd,nn_pb)
!=========================================================
! reduce dimension to save memory
	real* 8 Vm(ndf*nn_local, inner+1) ! Krylov space matrix
!        real* 8 Vm(ndf*nn, inner+1) ! Krylov space matrix
!=========================================================
	integer i,j,iouter,icount,INFO
	integer e1(inner+1)
	real* 8 x0(ndf*nn)
	real* 8 beta(inner+1)
	real* 8 eps
	real* 8 r0(ndf*nn)
	real* 8 rnorm, rnorm0,err
        real* 8 dv(ndf*nn)
	real* 8 Vy(ndf*nn)
        real* 8 vloc(ndf,nn), avloc(ndf*nn)

!==============================
! reduce dimension save memory
!	real* 8 temp(ndf*nn)
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



	eps = 1.0e-6
	e1(:) = 0
	e1(1) = 1
	x0(:) = 0
	iouter = 1
	r0(:) = bg(:)
	TRAN = 'N'
	av_tmp(:,:) = 0
	avloc(:) = 0
!	w(:) = 1.0d-4
        call getnorm_pa(r0,ndf,nn,node_local,nn_local,rnorm0)
        rnorm = sqrt(rnorm0)
        vloc(:,:)=0.0
if(myid==0)write(*,*)'rnorm=',rnorm
!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. 1.0e-6))

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
		    	dv((node-1)*ndf+jcount) = eps/w((node-1)*ndf+jcount)*Vm((icount-1)*ndf+jcount,j)
		    end do
		 end do!!!!!!!!!!calcule eps*inv(P)*V1
		
!===============================================
!		vloc(:,:) = 0.0d0
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
!call getnorm(vloc,vloc,ndf*nn,rnorm)
!if (myid ==0) write(*,*) 'after clear vloc should be zero', rnorm

!===============================================
		call equal_pa(dv,vloc,ndf,nn,node_local,nn_local)

!============================
!		call res_rotation_reverse(vloc,norm_node,nn_spbc,spbcnode)
		do icount=1,nn_local
			node=node_local(icount)
			vloc(1:ndf,node)=vloc(1:ndf,node)+d(1:ndf,node)
		end do


! Let vloc=vloc+d first then communicate, and then it should same # of loop (avoiding loop at the whole domain)
!=============================
!		call communicate_res(global_com,nn_global_com,local_com,nn_local_com,vloc,ndf,nn)
!**************used of periodical************!

if(nn_pb.gt.0) then
res_pb_temp(:,:)=0.0
res_pb(:,:)=0.0
vloc(1:nsd,pbnode(2,1:nn_pb))=0.0
res_pb_temp(1:nsd,1:nn_pb)=vloc(1:nsd,pbnode(1,1:nn_pb))
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(res_pb_temp(1,1),res_pb(1,1),nsd*nn_pb,mpi_double_precision, &
		mpi_sum,mpi_comm_world,ierror)
call mpi_barrier(mpi_comm_world,ierror)
end if

	        call communicate_res_ad(vloc,ndf,nn,send_address,ad_length)

if(nn_pb.gt.0) then
vloc(1:nsd,pbnode(2,1:nn_pb))=res_pb(1:nsd,1:nn_pb)
end if
!----------------------------------------------------------------------------------------------
!vloc(:,:)=vloc(:,:)+d(:,:)

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
!		av_tmp(:,:) = 0.0d0
!======================================





		call blockgmresnew(x,vloc,dold,av_tmp,hg,ien,fext,ne_local,ien_local,node_local,nn_local,&
				sur_fluid,I_fluid)

		call res_slipbc(av_tmp,vloc,x,rngface,ien,spbcnode,spbcele,ne_local,ien_local,I_fluid)
!                call communicate_res(global_com,nn_global_com,local_com,nn_local_com,av_tmp,ndf,nn)


!*****************periodical bc*************************!
if(nn_pb.gt.0) then
res_pb_temp(:,:)=0.0
res_pb(:,:)=0.0
res_pb_temp(1:nsd,1:nn_pb)=av_tmp(1:nsd,pbnode(2,1:nn_pb))
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(res_pb_temp(1,1),res_pb(1,1),nsd*nn_pb,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call mpi_barrier(mpi_comm_world,ierror)
end if

	        call communicate_res_ad(av_tmp,ndf,nn,send_address,ad_length)


if(nn_pb.gt.0) then
av_tmp(1:nsd,pbnode(1,1:nn_pb))=av_tmp(1:nsd,pbnode(1,1:nn_pb))+res_pb(1:nsd,1:nn_pb)
av_tmp(1:nsd,pbnode(2,1:nn_pb))=0.0
end if


!                call res_rotation(av_tmp,norm_node,nn_spbc,spbcnode)
!===================
! avloc(:)=0.0d0
!==================
		call equal_pa(av_tmp,avloc,ndf,nn,node_local,nn_local)


		do icount=1, nn_local
		node=node_local(icount)
			do jcount=1,ndf
				kcount=(node-1)*ndf+jcount
				avloc(kcount) = (-avloc(kcount)+bg(kcount))/eps ! get Av,bg=-r(u)

			end do
		end do

		!call setid(avloc,id,ndf)
                call setid_pa(avloc,ndf,nn,id,node_local,nn_local)
                space1(:)=0.0d0
                space2(:)=0.0d0
end_time=mpi_wtime()
	      do i=1,j
		 do icount=1,nn_local
			node=node_local(icount)
			do jcount=1,ndf
				kcount=(node-1)*ndf+jcount
				space1(i)=space1(i)+avloc(kcount)*Vm((icount-1)*ndf+jcount,i)
			end do
		 end do
!		call vector_dot_pa(avloc,Vm(:,i),ndf,nn,nn_local,node_local,space1(i))
		! Do the vector product of v_i * v_i set it to be h_i,j
	      end do  ! construct AVj and hi,j

                call mpi_barrier(mpi_comm_world,ierror)
                call mpi_allreduce(space1(1),space2(1),j,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
!                call mpi_bcast(space2(1),j,mpi_double_precision,0,mpi_comm_world,ierror)
end_time=mpi_wtime()-end_time
!if ((myid == 0) .and. (j==inner)) write(*,*) 'Time for get one hm vector', end_time

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
! get norm of Vm(:,j+1)
!	      call getnorm_pa(temp,ndf,nn,node_local,nn_local,rnorm0)

	      Hm(j+1,j) = sqrt(rnorm0)
	      do icount = 1, nn_local
!		node=node_local(icount)
		do jcount=1,ndf
			kcount=(icount-1)*ndf+jcount
 			Vm(kcount,j+1)=Vm(kcount,j+1)/Hm(j+1,j)
		end do
	      end do


	  end do  ! end inner loop

!if (myid == 0) write(*,*) 'Hm', Hm(1,1),Hm(2,1),Hm(2,2)
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
	   	dv(kcount) = eps/w(kcount)*x0(kcount)
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

!	vloc(:,:) = 0
!==============================


	call equal_pa(dv,vloc,ndf,nn,node_local,nn_local)

!	call res_rotation_reverse(vloc,norm_node,nn_spbc,spbcnode)
	do icount=1,nn_local
		node=node_local(icount)
		vloc(1:ndf,node)=vloc(1:ndf,node)+d(1:ndf,node)
	end do

!        call communicate_res(global_com,nn_global_com,local_com,nn_local_com,vloc,ndf,nn)
!********************periodical bc*****************!
if(nn_pb.gt.0) then
res_pb_temp(:,:)=0.0
res_pb(:,:)=0.0
!vloc(1:nsd,pbnode(2,1:nn_pb))=0.0
res_pb_temp(1:nsd,1:nn_pb)=vloc(1:nsd,pbnode(1,1:nn_pb))
vloc(1:nsd,pbnode(2,1:nn_pb))=0.0
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(res_pb_temp(1,1),res_pb(1,1),nsd*nn_pb,mpi_double_precision, &
                mpi_sum,mpi_comm_world,ierror)
call mpi_barrier(mpi_comm_world,ierror)
end if
        call communicate_res_ad(vloc,ndf,nn,send_address,ad_length)
if(nn_pb.gt.0) then
vloc(1:nsd,pbnode(2,1:nn_pb))=res_pb(1:nsd,1:nn_pb)
end if

!==============================
!	vloc(:,:) = vloc(:,:)+d(:,:)
!==============================
!===============================
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


	call blockgmresnew(x,vloc,dold,av_tmp,hg,ien,fext,ne_local,ien_local,node_local,nn_local,&
				sur_fluid,I_fluid)
!        call communicate_res(global_com,nn_global_com,local_com,nn_local_com,av_tmp,ndf,nn)
	call res_slipbc(av_tmp,vloc,x,rngface,ien,spbcnode,spbcele,ne_local,ien_local,I_fluid)
if(nn_pb.gt.0)then
res_pb_temp(:,:)=0.0
res_pb(:,:)=0.0
res_pb_temp(1:nsd,1:nn_pb)=av_tmp(1:nsd,pbnode(2,1:nn_pb))
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(res_pb_temp(1,1),res_pb(1,1),nsd*nn_pb,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call mpi_barrier(mpi_comm_world,ierror)
end if
        call communicate_res_ad(av_tmp,ndf,nn,send_address,ad_length)
if(nn_pb.gt.0) then
av_tmp(1:nsd,pbnode(1,1:nn_pb))=av_tmp(1:nsd,pbnode(1,1:nn_pb))+res_pb(1:nsd,1:nn_pb)
av_tmp(1:nsd,pbnode(2,1:nn_pb))=0.0
end if
!        call res_rotation(av_tmp,norm_node,nn_spbc,spbcnode)
!==================
!avloc(:)=0.0d0
!==================

       call equal_pa(av_tmp,avloc,ndf,nn,node_local,nn_local)
       
!       call setid(avloc,id,ndf)
       call setid_pa(avloc,ndf,nn,id,node_local,nn_local)

                do icount=1, nn_local
                node=node_local(icount)
                        do jcount=1,ndf
                                kcount=(node-1)*ndf+jcount
                                avloc(kcount) = (-avloc(kcount)+bg(kcount))/eps
                        end do
                end do


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
!	call mpi_bcast(dg(1),ndf*nn,mpi_double_precision,0,mpi_comm_world,ierror)
!	call res_rotation_reverse(dg,norm_node,nn_spbc,spbcnode)
	return
	end
