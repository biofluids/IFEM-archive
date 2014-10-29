

subroutine gmres_Laplace(x,d,w,bg,dg,ien,id,lp_el,n_el)
	use fluid_variables, only: nsd,nn,ne,nen
	implicit none
	integer inner, outer
	parameter (inner = 50) ! Laplace equation inner 50 should be sufficient
	parameter (outer = 5)  ! Laplace equation outer 5 should be sufficient
	integer ien(nen,ne)
	real* 8 x(nsd,nn)
	real* 8 d(nn)
	real* 8 bg(nn), dg(nn), w(nn)
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
	real* 8 Vm(nn, inner+1) ! Krylov space matrix
!--------------------------------
	integer lp_el(n_el)
	integer n_el
	integer id(nn)
!--------------------------------
	integer i,j,iouter,icount,INFO
	integer e1(inner+1)
	real* 8 x0(nn)
	real* 8 beta(inner+1)
	real* 8 eps
	real* 8 r0(nn)
	real* 8 rnorm, rnorm0,err
        real* 8 dv(nn)
	real* 8 Vy(nn)
        real* 8 avloc(nn)
	real* 8 temp(nn)
	character(1) TRAN
	real* 8 workls(2*inner)

	eps = 1.0e-6
	e1(:) = 0.0
	e1(1) = 1
	x0(:) = 0.0
	iouter = 1
	r0(:) = bg(:)
	TRAN = 'N'
	avloc(:) = 0.0
!	w(:) = 1
        call getnorm(r0,r0,nn,rnorm0)
        rnorm = sqrt(rnorm0)

!write(*,*) 'rnorm0', rnorm0

!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. eps))

	Vm(:,:) = 0.0
	do icount = 1, nn
	   Vm(icount,1) = r0(icount)/rnorm
	end do ! get V1

	   beta(:) = rnorm*e1(:) ! get beta*e1
	   Hm(:,:) = 0.0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
	   do j=1,inner
	  

		 do icount=1, nn
		    dv(icount) = 1/w(icount)*Vm(icount,j)
		 end do!!!!!!!!!!calcule eps*inv(P)*V1
		
		avloc(:) = 0.0d0
		call blockgmres_Laplace(x,dv,avloc,ien,lp_el,n_el)
!		avloc(:) = (-avloc(:)+bg(:))/eps ! get Av,bg=-r(u)
		call setid(avloc,id,1)
!do icount=1,nn
!                write(*,*) 'node av', icount,avloc(icount)
!end do
		continue
	      do i=1,j
		do icount = 1, nn
		   Hm(i,j)=Hm(i,j)+avloc(icount)*Vm(icount,i)
		 
		end do
	      end do  ! construct AVj and hi,j

	         
	      do icount = 1, nn
		 do i=1,j
		   Vm(icount,j+1) = Vm(icount,j+1)-Hm(i,j)*Vm(icount,i)
		 end do
		 Vm(icount,j+1)=Vm(icount,j+1)+avloc(icount)
	      end do  ! construct v(j+1)

	      do icount = 1, nn
		 temp(icount)=Vm(icount,j+1)
	      end do
	      call getnorm(temp, temp, nn,rnorm0)

	      Hm(j+1,j) = sqrt(rnorm0)

	      do icount = 1, nn
 		 Vm(icount,j+1)=Vm(icount,j+1)/Hm(j+1,j)
	      end do
	  end do  ! end inner loop
		
!==========================================================================
!	call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
	call givens(Hm,inner,beta)
!!!!!!!!!!!!beta(1:inner) is ym, the solution!!!!!!!!!!!!!!!!!!!!!!!!!
	Vy(:) = 0
	do icount=1,nn
	   do i=1,inner
	      Vy(icount)=Vy(icount)+Vm(icount,i)*beta(i)
	   end do
	   x0(icount)=x0(icount)+Vy(icount)
	end do ! calculate Xm
!write(*,*)'x0=',x0(:)
	do icount = 1, nn
	   dv(icount) = 1/w(icount)*x0(icount)
	end do

!write(*,*) dv(:)

	
	avloc(:) = 0.0d0
	call blockgmres_Laplace(x,dv,avloc,ien,lp_el,n_el)
!	avloc(:) = (-avloc(:)+bg(:))/eps
	call setid(avloc,id,1)

!!!!!!!!!!calculate AXm
	do icount=1,nn
	   r0(icount) = bg(icount)-avloc(icount)
	end do !update r0=f-AX0
	call getnorm(r0,r0,nn,rnorm0)
	err = sqrt(rnorm0)

	rnorm = sqrt(rnorm0)
	iouter = iouter + 1        

!	write(*,*) '***** Indicator Laplace Eq err******=',err
111     continue  ! end outer loop
!write(*,*)'x0_2=',x0(:)

	dg(:) = 1/w(:)*x0(:) ! unscaled x0
!write(*,*)'dg=',dg(:)
!write(*,*)'w=',w(:)
!stop
	return
	end
