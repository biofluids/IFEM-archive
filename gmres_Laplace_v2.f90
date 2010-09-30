

!subroutine gmres_Laplace(x,d,w,bg,dg,ien,inter_ele,ne_inter,rng,ien_center)
subroutine gmres_Laplace_v2(x,d,w,bg,dg,ne_outer,outer_ele,ne_inner,inner_ele,ien_center)
	use fluid_variables, only: nsd,ne,inner,outer
	use centermesh_variables
	use interface_variables
	implicit none
	
	real* 8 x(nsd,nn_center)
	real* 8 d(nn_center)
	real* 8 bg(nn_center), dg(nn_center), w(nn_center)
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
	real* 8 Vm(nn_center, inner+1) ! Krylov space matrix
	integer ien_center(nen_center,ne_center)
!	integer ne_inter
!	integer inter_ele(ne)
	integer ne_inner,ne_outer
	integer inner_ele(ne),outer_ele(ne)
!	real* 8 flag

	integer i,j,iouter,icount,INFO
	integer e1(inner+1)
	real* 8 x0(nn_center)
	real* 8 beta(inner+1)
	real* 8 eps
	real* 8 r0(nn_center)
	real* 8 rnorm, rnorm0,err
        real* 8 dv(nn_center)
	real* 8 Vy(nn_center)
        real* 8 avloc(nn_center)
	real* 8 temp(nn_center)
	character(1) TRAN
	real* 8 workls(2*inner)

	eps = 1.0e-6
	e1(:) = 0
	e1(1) = 1
	x0(:) = 0
	iouter = 1
	r0(:) = bg(:)
	TRAN = 'N'
	avloc(:) = 0
!	w(:) = 1
!	inner=20
!	outer=5
        call getnorm(r0,r0,nn_center,rnorm0)
        rnorm = sqrt(rnorm0)


!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
	do 111, while((iouter .le. outer) .and. (rnorm .ge. 1.0e-9))

	Vm(:,:) = 0
	do icount = 1, nn_center
	   Vm(icount,1) = r0(icount)/rnorm
	end do ! get V1

	   beta(:) = rnorm*e1(:) ! get beta*e1
	   Hm(:,:) = 0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
	   do j=1,inner
	  

		 do icount=1, nn_center
		    dv(icount) = 1/w(icount)*Vm(icount,j)
		 end do!!!!!!!!!!calcule eps*inv(P)*V1
		
		avloc(:) = 0.0d0
		call blockgmres_Laplace(x,dv,avloc,ien_center)
!		avloc(:) = (-avloc(:)+bg(:))/eps ! get Av,bg=-r(u)
!		flag = 0.0
!		call form_inter_ele(inter_ele,ne_inter,avloc,ien,flag)
!		call form_inter_bc(avloc,rng,ien,flag) !set bc
		do icount=1,ne_outer
		  avloc(outer_ele(icount))=0.0
		end do
		do icount=1,ne_inner
		  avloc(inner_ele(icount))=0.0
		end do


	      do i=1,j
		do icount = 1, nn_center
		   Hm(i,j)=Hm(i,j)+avloc(icount)*Vm(icount,i)
		 
		end do
	      end do  ! construct AVj and hi,j

	         
	      do icount = 1, nn_center
		 do i=1,j
		   Vm(icount,j+1) = Vm(icount,j+1)-Hm(i,j)*Vm(icount,i)
		 end do
		 Vm(icount,j+1)=Vm(icount,j+1)+avloc(icount)
	      end do  ! construct v(j+1)

	      do icount = 1, nn_center
		 temp(icount)=Vm(icount,j+1)
	      end do
	      call getnorm(temp, temp, nn_center,rnorm0)

	      Hm(j+1,j) = sqrt(rnorm0)

	      do icount = 1, nn_center
 		 Vm(icount,j+1)=Vm(icount,j+1)/Hm(j+1,j)
	      end do
	  end do  ! end inner loop
		
	call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
!!!!!!!!!!!!beta(1:inner) is ym, the solution!!!!!!!!!!!!!!!!!!!!!!!!!
	Vy(:) = 0
	do icount=1,nn_center
	   do i=1,inner
	      Vy(icount)=Vy(icount)+Vm(icount,i)*beta(i)
	   end do
	   x0(icount)=x0(icount)+Vy(icount)
	end do ! calculate Xm
!write(*,*)'x0=',x0(:)
	do icount = 1, nn_center
	   dv(icount) = 1/w(icount)*x0(icount)
	end do

!write(*,*) dv(:)

	
	avloc(:) = 0.0d0
	call blockgmres_Laplace(x,dv,avloc,ien_center)
!	avloc(:) = (-avloc(:)+bg(:))/eps

!	flag = 0.0
!	call form_inter_ele(inter_ele,ne_inter,avloc,ien,flag)
!	call form_inter_bc(avloc,rng,ien,flag)
        do icount=1,ne_outer
           avloc(outer_ele(icount))=0.0
        end do
        do icount=1,ne_inner
           avloc(inner_ele(icount))=0.0
        end do

!!!!!!!!!!calculate AXm
	do icount=1,nn_center
	   r0(icount) = bg(icount)-avloc(icount)
	end do !update r0=f-AX0
	call getnorm(r0,r0,nn_center,rnorm0)
	err = sqrt(rnorm0)

	rnorm = sqrt(rnorm0)
	iouter = iouter + 1        

	write(*,*) 'err_2=',err
111     continue  ! end outer loop
!write(*,*)'x0_2=',x0(:)

	dg(:) = 1/w(:)*x0(:) ! unscaled x0
!write(*,*)'dg=',dg(:)
!write(*,*)'w=',w(:)
!stop
	return
	end
