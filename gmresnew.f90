

subroutine gmres(x,d,dold,w,bg,dg,hg,ien,fext,id)
	use fluid_variables, only: nsd,nn,ne,nen,ndf,inner,outer
	implicit none

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
        call getnorm(r0,r0,ndf*nn,rnorm0)
        rnorm = sqrt(rnorm0)


!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. eps))

	Vm(:,:) = 0
	do icount = 1, ndf*nn
	   Vm(icount,1) = r0(icount)/rnorm
	end do ! get V1

	   beta(:) = rnorm*e1(:) ! get beta*e1
	   Hm(:,:) = 0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
	   do j=1,inner
	  

		 do icount=1, ndf*nn
		    dv(icount) = eps/w(icount)*Vm(icount,j)
		 end do!!!!!!!!!!calcule eps*inv(P)*V1
		vloc(:,:) = 0.0d0
		call equal(dv,vloc,ndf*nn)
	 	
!open(unit=8436, file='res3.out', status='unknown')
!write(8436,*) Vm(1:ndf*nn,j)


		vloc(:,:) = vloc(:,:)+d(:,:)  ! calculate u+eps*inv(P)*V1
!open(unit=8435, file='res2.out', status='unknown')
!write(8435,*) vloc(1:ndf,1:nn)
!stop
		av_tmp(:,:) = 0.0d0
		call blockgmresnew(x,vloc,dold,av_tmp,hg,ien,fext)

!open(unit=8434, file='res1.out', status='unknown')
!write(8434,*) av_tmp(1:ndf,1:nn)
		call equal(av_tmp,avloc,ndf*nn)

!open(unit=8433, file='res.out', status='unknown')
!write(8433,*) avloc(1:ndf*nn)
!stop

		avloc(:) = (-avloc(:)+bg(:))/eps ! get Av,bg=-r(u)

!open(unit=8433, file='res.out', status='unknown')
!write(8433,*) avloc(1:ndf*nn)


		call setid(avloc,id,ndf)

	      do i=1,j
		do icount = 1, ndf*nn
		   Hm(i,j)=Hm(i,j)+avloc(icount)*Vm(icount,i)
		 
		end do
	      end do  ! construct AVj and hi,j
	         
	      do icount = 1, ndf*nn
!	call equal(temp,dg,ndf*nn)
		 do i=1,j
		   Vm(icount,j+1) = Vm(icount,j+1)-Hm(i,j)*Vm(icount,i)
		 end do
		 Vm(icount,j+1)=Vm(icount,j+1)+avloc(icount)
	      end do  ! construct v(j+1)

	      do icount = 1, ndf*nn
		 temp(icount)=Vm(icount,j+1)
	      end do
	      call getnorm(temp, temp, ndf*nn,rnorm0)

	      Hm(j+1,j) = sqrt(rnorm0)

	      do icount = 1, ndf*nn
 		 Vm(icount,j+1)=Vm(icount,j+1)/Hm(j+1,j)
	      end do
	  end do  ! end inner loop

!write(*,*) 'Hm', Hm(:,:)
!====================================================================================		
! Use lapack to solve the LS problem
!	call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
!===================================================================================
! Use Givens rotation to solve the LS problem
        call givens(Hm,inner,beta)


!!!!!!!!!!!!beta(1:inner) is ym, the solution!!!!!!!!!!!!!!!!!!!!!!!!!
	Vy(:) = 0
	do icount=1,ndf*nn
	   do i=1,inner
	      Vy(icount)=Vy(icount)+Vm(icount,i)*beta(i)
	   end do
	   x0(icount)=x0(icount)+Vy(icount)
	end do ! calculate Xm
	do icount = 1, ndf*nn
	   dv(icount) = eps/w(icount)*x0(icount)
	end do

!write(*,*) dv(:)

	vloc(:,:) = 0
	call equal(dv,vloc,ndf*nn)
	vloc(:,:) = vloc(:,:)+d(:,:)

	av_tmp(:,:) = 0.0d0
	call blockgmresnew(x,vloc,dold,av_tmp,hg,ien,fext)
    call equal(av_tmp,avloc,ndf*nn)
        call setid(avloc,id,ndf)

	avloc(:) = (-avloc(:)+bg(:))/eps
!!!!!!!!!!calculate AXm
	do icount=1,ndf*nn
	   r0(icount) = bg(icount)-avloc(icount)
	end do !update r0=f-AX0

	call getnorm(r0,r0,ndf*nn,rnorm0)
	err = sqrt(rnorm0)

	rnorm = sqrt(rnorm0)
	iouter = iouter + 1        

	write(*,*) 'err=',err
111     continue  ! end outer loop


	dg(:) = 1/w(:)*x0(:) ! unscaled x0

!open(unit=8433, file='res.out', status='unknown')
!write(8433,*) dg(1:ndf*nn)
!stop

	return
	end
