

subroutine gmres(x,d,dold,w,bg,dg,hg,ien,fext,id, &
		ne_local,ien_local,mdata,n_mdata,node_local,nn_local, &
		global_com,nn_global_com,local_com,nn_local_com)
	use fluid_variables, only: nsd,nn,ne,nen,ndf,inner,outer
 	use solid_variables, only: nn_solid
        use mpi_variables
	implicit none
      include 'mpif.h'
	real* 8 x(nsd,nn),id(ndf,nn)
	real* 8 d(ndf,nn), dold(ndf,nn),hg(ne),fext(ndf,nn),ien(nen,ne)
	real* 8 bg(ndf*nn), dg(ndf*nn), w(ndf*nn)
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
	real* 8 Vm(ndf*nn, inner+1) ! Krylov space matrix

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
	real* 8 temp(ndf*nn)
	character(1) TRAN
	real* 8 workls(2*inner)
	real* 8 av_tmp(ndf,nn)

!---------------------------------------
  integer mdata(nn_solid)
  integer n_mdata
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
!	w(:) = 1
        call getnorm_pa(r0,ndf,nn,node_local,nn_local,rnorm0)
        rnorm = sqrt(rnorm0)

!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. eps))

	Vm(:,:) = 0.0d0
	do icount = 1, nn_local
	   node=node_local(icount)
	   do jcount=1,ndf
	   	Vm((node-1)*ndf+jcount,1) = r0((node-1)*ndf+jcount)/rnorm
	   end do
	end do ! get V1

	   beta(:) = rnorm*e1(:) ! get beta*e1
	   Hm(:,:) = 0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
	   do j=1,inner
	  

		 do icount=1, nn_local
		    node=node_local(icount)
		    do jcount=1,ndf
		    	dv((node-1)*ndf+jcount) = eps/w((node-1)*ndf+jcount)*Vm((node-1)*ndf+jcount,j)
		    end do
		 end do!!!!!!!!!!calcule eps*inv(P)*V1
		vloc(:,:) = 0.0d0

		call equal_pa(dv,vloc,ndf,nn,node_local,nn_local)



!           call mpi_barrier(mpi_comm_world,ierror)
!av_tmp(:,:)=0.0d0
!           call mpi_reduce(vloc(1,1),av_tmp(1,1),ndf*nn,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
!if (myid == 0) then
! open(unit=8436, file='res3.out', status='unknown')
!write(8436,*) Vm(1:ndf*nn,j)
!stop
!end if



		call communicate_res(global_com,nn_global_com,local_com,nn_local_com,vloc,ndf,nn)
!----------------------------------------------------------------------------------------------
! ???????????????????????????????????????????????????????????


!		do icount=1, nn_local ! For nodes on its own proc only and not shared by any other procs
!		   node=node_local(icount)
!		   do kcount=1, nn_local_com
!			if (global_com(local_com(kcount)) == node) then
!			   flag = 0
!			end if
!		   end do
!
!		   if (flag .ne. 0) then
!			   do jcount=1,ndf
!				vloc(jcount,node) = vloc(jcount,node)+d(jcount,node)  ! calculate u+eps*inv(P)*V1
!		   	   end do
!		    end if
!		end do
!
!		do icount=1, nn_local_com ! For nodes on shared by other procs
!			node=global_com(local_com(icount))
!			do jcount=1,ndf
!				vloc(jcount,node)=vloc(jcount,node)+d(jcount,node)
!			end do
!		end do

vloc(:,:)=vloc(:,:)+d(:,:)
!if (myid == 0) then
! open(unit=8435, file='res2.out', status='unknown')
!write(8435,*) vloc(1:ndf,1:nn)
!end if


!----------------------------------------------------------------------------------------------
! ???????????????????????????????????????????????????????????
		av_tmp(:,:) = 0.0d0


		call blockgmresnew(x,vloc,dold,av_tmp,hg,ien,fext,ne_local,ien_local,mdata,n_mdata,node_local,nn_local)




                call communicate_res(global_com,nn_global_com,local_com,nn_local_com,av_tmp,ndf,nn)
!           call mpi_barrier(mpi_comm_world,ierror)
!           call mpi_reduce(av_tmp(1,1),avloc(1),ndf*nn,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
!if (myid == 0) then
! open(unit=8434, file='res1.out', status='unknown')
!write(8434,*) av_tmp(1:ndf,1:nn)
!end if
avloc(:)=0.0d0
		call equal_pa(av_tmp,avloc,ndf,nn,node_local,nn_local)

!if (myid == 0) then
! open(unit=8433, file='res.out', status='unknown')
!write(8433,*) avloc(1:ndf*nn)
!stop
!end if



		do icount=1, nn_local
		node=node_local(icount)
			do jcount=1,ndf
				kcount=(node-1)*ndf+jcount
				avloc(kcount) = (-avloc(kcount)+bg(kcount))/eps ! get Av,bg=-r(u)
			end do
		end do
!if (myid == 0) then
! open(unit=8433, file='res.out', status='unknown')
!write(8433,*) avloc(1:ndf*nn)
!end if

		call setid(avloc,id,ndf)
                space1(:)=0.0d0
                space2(:)=0.0d0
end_time=mpi_wtime()
	      do i=1,j
		call vector_dot_pa(avloc,Vm(:,i),ndf,nn,nn_local,node_local,space1(i))
		! Do the vector product of v_i * v_i set it to be h_i,j
	      end do  ! construct AVj and hi,j

                call mpi_barrier(mpi_comm_world,ierror)
                call mpi_reduce(space1(1),space2(1),j,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
                call mpi_bcast(space2(1),j,mpi_double_precision,0,mpi_comm_world,ierror)
end_time=mpi_wtime()-end_time
if (myid == 0) write(*,*) 'Time for get one hm vector', end_time

                Hm(1:j,j)=space2(1:j)

	      do icount = 1, nn_local
		node=node_local(icount)
		do jcount=1,ndf
			kcount=(node-1)*ndf+jcount
		
		 	do i=1,j
		   		Vm(kcount,j+1) = Vm(kcount,j+1)-Hm(i,j)*Vm(kcount,i)
		 	end do
		 	Vm(kcount,j+1)=Vm(kcount,j+1)+avloc(kcount)
	        end do  
	      end do  ! construct v(j+1)

	      do icount = 1, nn_local
		node=node_local(icount)
		do jcount=1,ndf
			kcount=(node-1)*ndf+jcount
		 	temp(kcount)=Vm(kcount,j+1)
		end do
	      end do
	      call getnorm_pa(temp,ndf,nn,node_local,nn_local,rnorm0)

	      Hm(j+1,j) = sqrt(rnorm0)

	      do icount = 1, nn_local
		node=node_local(icount)
		do jcount=1,ndf
			kcount=(node-1)*ndf+jcount
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
	       		 Vy(kcount)=Vy(kcount)+Vm(kcount,i)*beta(i)
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

	vloc(:,:) = 0
	call equal_pa(dv,vloc,ndf,nn,node_local,nn_local)
        call communicate_res(global_com,nn_global_com,local_com,nn_local_com,vloc,ndf,nn)
	vloc(:,:) = vloc(:,:)+d(:,:)

	av_tmp(:,:) = 0.0d0
	call blockgmresnew(x,vloc,dold,av_tmp,hg,ien,fext,ne_local,ien_local,mdata,n_mdata,node_local,nn_local)
        call communicate_res(global_com,nn_global_com,local_com,nn_local_com,av_tmp,ndf,nn)
avloc(:)=0.0d0
       call equal_pa(av_tmp,avloc,ndf,nn,node_local,nn_local)
       
       call setid(avloc,id,ndf)

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
	call mpi_reduce(dg_sent(1),dg(1),ndf*nn,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)
	call mpi_bcast(dg(1),ndf*nn,mpi_double_precision,0,mpi_comm_world,ierror)
        call mpi_barrier(mpi_comm_world,ierror)


	return
	end
