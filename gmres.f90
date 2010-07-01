!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c GMRES.f
!c Lucy Zhang
!c Northwestern University
!c Iterative Solver
!c output increment of d, dg
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine gmres(x,d,dold,id,w,bg,dg,hg,ien, &
	    z,v,zg,avg,sm,vloc,avloc,h,y,cc,ss,fext)
      use fluid_variables, only: nsd,nn,ne,nen,ndf,inner,outer,iscaling
	implicit none

	real* 8 x(nsd,nn)
	real* 8 d(ndf,nn), dold(ndf,nn),hg(ne),dd(ndf,nn)
	real* 8 bg(ndf*nn), dg(ndf*nn), w(ndf*nn)
	integer id(ndf,nn),ien(nen,ne)
	real* 8 h(inner+1,inner)
	real* 8 y(inner+1)
	real* 8 cc(inner), ss(inner)
	real* 8 z(ndf*nn,inner)
	real* 8 v(ndf*nn,inner+1)
	real* 8 zg(ndf*nn), avg(ndf*nn), sm(ndf*nn)
	real* 8 vloc(ndf,nn),avloc(ndf,nn)
	logical assemble
	real* 8 eps, rnorm, rnorm0, order
	real* 8 gam,hsave,ysave,tmpo
	integer i,j,k,ij,jj,i1,j1,k1,l,iqc,igmres
	real*8 fext(nsd,nn)

	eps = 1.0e-12
	assemble = .true.

	do iqc=1,ndf*nn
	   sm(iqc) = 1.0
	enddo

	if(iscaling.eq.0) then
	   do iqc=1,ndf*nn
	      w(iqc) = 1.0
	   enddo
	else
	   do iqc=1,ndf*nn
	      if (w(iqc).lt.0.0) sm(iqc) = -1.0
	      w(iqc) = 1.0 / sqrt(abs(w(iqc)))
	   enddo
	endif

!c       clear arrays
	v(1:ndf*nn,1:inner+1) = 0.0d0

	do iqc=1,ndf*nn
	   v(iqc,1) = bg(iqc)
	  ! v(iqc,1) = w(iqc) * v(iqc,1)
	enddo
	call getnorm(v(1,1),v(1,1),ndf*nn,rnorm0)
	rnorm0 = sqrt(rnorm0)

!c       outer GMRES loop (igmres)
	igmres = 0 
 10	igmres = igmres + 1

!c       convergence check
	call getnorm(v(1,1),v(1,1),ndf*nn,rnorm)
	rnorm = sqrt(rnorm)
	if (rnorm.le.eps.or.igmres.gt.outer) goto 700

!c       first krylov vector
	do iqc=1,ndf*nn
	   v(iqc,1) = v(iqc,1) / rnorm
	enddo

!c       construct krylov vectors 2:inner and hessenberg matrix
	do j=1,inner

!c       insert contribution of preconditioner
	   do iqc=1,ndf*nn
	      zg(iqc) = sm(iqc) * v(iqc,j)
	      z(iqc,j) = zg(iqc)
	   enddo

!c       compute A * v_j
	   do iqc=1,ndf*nn
	      zg(iqc) = w(iqc) * zg(iqc)
	   enddo

	   call equal(zg,vloc,ndf*nn)

       avloc(1:ndf,1:nn) = 0.0d0
       dd(:,:)=d(:,:)+eps*vloc(:,:)
!	   call blockgmres(x,d,dold,vloc,avloc,hg,ien,fext)
           call blockgmresnew(x,dd,dold,avloc,hg,ien,fext)
            
	   call equal(avloc,avg,ndf*nn)
           avg(:)=(-avg(:)+bg(:))/eps
	   call setid(avg,id,ndf)

	   do iqc=1,ndf*nn
	   !   avg(iqc) = w(iqc) * avg(iqc)
	      v(iqc,j+1) = avg(iqc)
	   enddo
	   do i=1,j
	      call getnorm(v(1,j+1),v(1,i),ndf*nn,tmpo)
	      h(i,j) = tmpo
	      do iqc=1,ndf*nn
		 v(iqc,j+1) = v(iqc,j+1) - tmpo * v(iqc,i) 
	      enddo
	   enddo

	   call getnorm(v(1,j+1),v(1,j+1),ndf*nn,tmpo)
	   tmpo = sqrt(tmpo)
	   h(j+1,j) = tmpo
	   do iqc=1,ndf*nn
	      v(iqc,j+1) = v(iqc,j+1) / tmpo
	   enddo

	enddo			! end inner loop

!c       initialize reduced residual
	y(1) = rnorm
	do i=2,inner+1
	   y(i) = 0.0
	enddo

!c       rotations
	do j=1,inner

	   j1 = j + 1

!c       previously computed rotations on colundf j
	   do i=2,j
	      i1 = i - 1
	      hsave = h(i1,j)
	      h(i1,j) = + cc(i1) * hsave + ss(i1) * h(i,j)
	      h(i ,j) = - ss(i1) * hsave + cc(i1) * h(i,j)
	   end do

!c       new rotation on colundf j
	   gam = sqrt(h(j,j)**2 + h(j1,j)**2)
	   cc(j) = h(j,j) / gam
	   ss(j) = h(j1,j) / gam
	   h(j,j) = cc(j) * h(j,j) + ss(j) * h(j1,j)
!c       note: under-diagonal term h(j+1,j) becomes 0

	   y(j1) = - ss(j) * y(j)
	   y(j ) = + cc(j) * y(j)
!c       note: till now y(j+1) = 0
	   rnorm = abs(y(j1))
	end do


!c       back substitution
	j = inner		!should reach here straight from rotation loop
	y(j) = y(j)/h(j,j)
	do jj=2,j
	   k = j - jj + 1
	   k1 = k + 1
	   ysave = y(k)
	   do l=k1,j
	      ysave = ysave - h(k,l) * y(l)
	   end do
	   y(k) = ysave / h(k,k)
	end do

!c       compute dg iterate dg = dg + Z_m * y
	j = inner		!(PVM only fix)
	do jj=1,j
	   tmpo = y(jj)
	   do iqc=1,ndf*nn
	      dg(iqc) = dg(iqc) + tmpo * z(iqc,jj)
	   enddo
	end do

!c       if not done recover residual vector
	if (igmres.le.outer) then
	   do jj=1,j
	      ij = j - jj + 2
	      y(ij-1) = -ss(ij-1)*y(ij)
	      y(ij) = cc(ij-1)*y(ij)
	   end do
	   do jj=1,j+1
	      tmpo = y(jj)
	      if (jj.eq.1) tmpo = tmpo - 1.0
	      do iqc=1,ndf*nn
		 v(iqc,1) = v(iqc,1) + tmpo*v(iqc,jj)
	      enddo
	   enddo
	endif

	goto 10

 700	continue

!c       go back to unscaled system
	do iqc=1,ndf*nn
	   dg(iqc) = w(iqc) * dg(iqc)
	enddo

	order = 0.43429448 * log(rnorm0/rnorm)
	write(6,101) order,rnorm0,rnorm
	write(7,101) order,rnorm0,rnorm

 101	format('Flow     : convergence order = ', f5.2, 2e14.7)
 102	format('Flow     : zero residual')

	return
	end
