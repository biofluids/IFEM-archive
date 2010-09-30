

subroutine gmres_normal(x,w,bg,dg,hg,ien,id_inter,I_fluid)
	use fluid_variables, only: nsd,nn,ne,nen,inner,outer
	implicit none

	real* 8 x(nsd,nn)
        integer id_inter(nn)
	real* 8 hg(ne),ien(nen,ne)
	real* 8 bg(nsd*nn), dg(nsd*nn), w(nsd*nn)
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
	real* 8 Vm(nsd*nn, inner+1) ! Krylov space matrix
	real* 8 I_fluid(nn)

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
	real* 8 temp(nsd*nn)
	character(1) TRAN
	real* 8 workls(2*inner)
	real* 8 av_tmp(nsd,nn)

	eps = 1.0e-6
	e1(:) = 0
	e1(1) = 1
	x0(:) = 0
	iouter = 1
	r0(:) = bg(:)
	TRAN = 'N'
	av_tmp(:,:) = 0
	avloc(:) = 0
	w(:) = 1
        call getnorm(r0,r0,nsd*nn,rnorm0)
        rnorm = sqrt(rnorm0)

!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. eps))

	Vm(:,:) = 0
	do icount = 1, nsd*nn
	   Vm(icount,1) = r0(icount)/rnorm
	end do ! get V1

	   beta(:) = rnorm*e1(:) ! get beta*e1
	   Hm(:,:) = 0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
	   do j=1,inner
	  

		 do icount=1, nsd*nn
		    dv(icount) = 1/w(icount)*Vm(icount,j)
		 end do!!!!!!!!!!calculeinv(P)*V1
		av_tmp(:,:) = 0.0d0
		call blockgmres_norm(x,dv,av_tmp,ien,hg)
		call equal(av_tmp,avloc,nsd*nn)
		call set_id_inter(avloc,id_inter)

	      do i=1,j
		do icount = 1, nsd*nn
		   Hm(i,j)=Hm(i,j)+avloc(icount)*Vm(icount,i)
		 
		end do
	      end do  ! construct AVj and hi,j

	         
	      do icount = 1, nsd*nn
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
		
	call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
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


	av_tmp(:,:) = 0.0d0
	call blockgmres_norm(x,dv,av_tmp,ien,hg)
	call equal(av_tmp,avloc,nsd*nn)
	call set_id_inter(avloc,id_inter)

!	avloc(:) = (-avloc(:)+bg(:))/eps
!!!!!!!!!!calculate AXm
	do icount=1,nsd*nn
	   r0(icount) = bg(icount)-avloc(icount)
	end do !update r0=f-AX0

	call getnorm(r0,r0,nsd*nn,rnorm0)
	err = sqrt(rnorm0)

	rnorm = sqrt(rnorm0)
	iouter = iouter + 1        

	write(*,*) 'err=',err
111     continue  ! end outer loop


	dg(:) = 1/w(:)*x0(:) ! unscaled x0

	return
	end
