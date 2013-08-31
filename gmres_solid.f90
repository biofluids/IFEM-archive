

subroutine gmres_solid(x,w,bg,dg,ien,id,nsd,nn,ne,nen,inner,outer,nquad,wq,sq,x_pre1,&
			solid_prevel,solid_preacc,solid_stress,ne_sbc,ien_sbc,mtype)
	use mpi_variables, only: myid
	implicit none
! use GMRES to solve mesh update equation  with diagonal preconditioner
	real* 8 x(nsd,nn)
        integer id(nsd,nn)
	real* 8 ien(ne,nen)
	real* 8 bg(nsd*nn), dg(nsd*nn), w(nsd*nn)
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
	real* 8 Vm(nsd*nn, inner+1) ! Krylov space matrix
	real(8) x_pre1(nsd,nn)
	real(8) solid_prevel(nsd,nn)
	real(8) solid_preacc(nsd,nn)
!	real(8) x_pre2(nsd,nn)
        real(8) solid_stress(nsd*2,nn) ! fluid stress acting on the solid boundary including pressure and viscous stress
        integer ne_sbc                  ! solid elements on the boundary
        integer ien_sbc(ne_sbc,nen+2)
	integer mtype(ne)
!-----------------------------------------------
	integer nsd
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
	real* 8 e1(inner+1)
	real* 8 x0(nsd*nn)
	real* 8 beta(inner+1)
	real* 8 eps
	real* 8 r0(nsd*nn)
	real* 8 rnorm, rnorm0,err
        real* 8 dv(nsd*nn)
	real* 8 Vy(nsd*nn)
        real* 8 vloc(nsd,nn), avloc(nsd*nn)
	real* 8 temp(nsd*nn)
	character(1) TRAN
	real* 8 workls(2*inner)
	real* 8 av_tmp(nsd,nn)

	integer time_arrary_0(8)
        integer time_arrary_1(8)
	real(8) start_time
	real(8) end_time
	real(8) err_givens
	eps = 1.0e-6
	e1(:) = 0
	e1(1) = 1
	x0(:) = 0
	iouter = 1
	r0(:) = bg(:)
	TRAN = 'N'
	av_tmp(:,:) = 0.0
	avloc(:) = 0.0
	dv(:) =0.0
!	w(:) = 1.0
        call getnorm(r0,r0,nsd*nn,rnorm0)
        rnorm = sqrt(rnorm0)


!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. 1.0e-5))

	Vm(:,:) = 0
	do icount = 1, nsd*nn
	   Vm(icount,1) = r0(icount)/rnorm
	end do ! get V1

	   beta(:) = rnorm*e1(:) ! get beta*e1
	   Hm(:,:) = 0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
        call date_and_time(values=time_arrary_0)  
	   do j=1,inner
	  

		 do icount=1, nsd*nn
		    dv(icount) = 1/w(icount)*Vm(icount,j)
		 end do!!!!!!!!!!calcule eps*inv(P)*V1
		vloc(:,:) = 0.0d0
		call equal(dv,vloc,nsd*nn)
		av_tmp(:,:) = 0.0d0
		call blockgmres_solid(x,vloc,av_tmp,ien,nsd,nen,ne,nn,nquad,wq,sq,x_pre1,&
		solid_prevel,solid_preacc,ien_sbc,ne_sbc,solid_stress,mtype)
		call equal(av_tmp,avloc,nsd*nn)
		call setsolid_id(avloc,id,nsd)
	      do i=1,j
		do icount = 1, nsd*nn
		   Hm(i,j)=Hm(i,j)+avloc(icount)*Vm(icount,i)
		 
		end do
	      end do  ! construct AVj and hi,j

	         
	      do icount = 1, nsd*nn
	! call equal(temp,dg,nsd*nn)
		 do i=1,j
		   Vm(icount,j+1) = Vm(icount,j+1)-Hm(i,j)*Vm(icount,i)
		 end do
		 Vm(icount,j+1)=Vm(icount,j+1)+avloc(icount)
	      end do  ! construct v(j+1)

	      do icount = 1, nsd*nn
		 temp(icount)=Vm(icount,j+1)
	      end do
	      call getnorm(temp, temp, nsd*nn,rnorm0)

	      Hm(j+1,j) = sqrt(rnorm0)

	      do icount = 1, nsd*nn
 		 Vm(icount,j+1)=Vm(icount,j+1)/Hm(j+1,j)
	      end do
	  end do  ! end inner loop
	call date_and_time(values=time_arrary_1)
	start_time=time_arrary_0(5)*3600+time_arrary_0(6)*60+time_arrary_0(7)+time_arrary_0(8)*0.001
        end_time=time_arrary_1(5)*3600+time_arrary_1(6)*60+time_arrary_1(7)+time_arrary_1(8)*0.001
!        write(*,*) 'Inner loop time', end_time-start_time


	call date_and_time(values=time_arrary_0)
!	write(*,*) 'hm-hgivens', Hm(:,:)
	call givens(Hm,inner,beta)
!	call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
	
!	err_givens=maxval(abs(beta(1:inner)-x_givens(1:inner)))
!	write(*,*) 'givens error', err_givens
        call date_and_time(values=time_arrary_1)
        start_time=time_arrary_0(5)*3600+time_arrary_0(6)*60+time_arrary_0(7)+time_arrary_0(8)*0.001
        end_time=time_arrary_1(5)*3600+time_arrary_1(6)*60+time_arrary_1(7)+time_arrary_1(8)*0.001
!	write(*,*) 'LS time by lapack', end_time-start_time
!!!!!!!!!!!!beta(1:inner) is ym, the solution!!!!!!!!!!!!!!!!!!!!!!!!!
	Vy(:) = 0
	do icount=1,nsd*nn
	   do i=1,inner
	      Vy(icount)=Vy(icount)+Vm(icount,i)*beta(i)
	   end do
	   x0(icount)=x0(icount)+Vy(icount)
	end do ! calculate Xm
	do icount = 1, nsd*nn
	   dv(icount) = 1/w(icount)*x0(icount)
	end do

	vloc(:,:) = 0
	call equal(dv,vloc,nsd*nn)
!	vloc(:,:) = vloc(:,:)+d(:,:)
	av_tmp(:,:) = 0.0d0
        call blockgmres_solid(x,vloc,av_tmp,ien,nsd,nen,ne,nn,nquad,wq,sq,x_pre1,&
		solid_prevel,solid_preacc,ien_sbc,ne_sbc,solid_stress,mtype)
    call equal(av_tmp,avloc,nsd*nn)
        call setsolid_id(avloc,id,nsd)

!	avloc(:) = (-avloc(:)+bg(:))/eps
!!!!!!!!!!calculate AXm
	do icount=1,nsd*nn
	   r0(icount) = bg(icount)-avloc(icount)
	end do !update r0=f-AX0

	call getnorm(r0,r0,nsd*nn,rnorm0)
	err = sqrt(rnorm0)

	rnorm = sqrt(rnorm0)
	iouter = iouter + 1        
if (myid == 0) then
	write(*,*) 'In solid dispment solver gmres err=',err
end if
111     continue  ! end outer loop



	dg(:) = 1/w(:)*x0(:) ! unscaled x0
!write(*,*) 'x0', x0(:)
	return
	end
