c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine gmresf(x,u,id,w,bg,dg,hg,ien,hn,hm,
     &                 z,v,zg,avg,sm,vloc,avloc,h,y,cc,ss)

      implicit none
	include "global.h"

      real* 8 x(nsd,nn_loc),u(nsd,nn_loc)
	real* 8 hg(nec)
	real* 8 bg(nnc), dg(nnc), w(nnc)
      integer id(nnc),ien(nen,nec)
	real* 8 hn(nnc),hm(nn_loc)

	real* 8 h(inner+1,inner)
	real* 8 y(inner+1)
	real* 8 cc(inner), ss(inner)
	
	real* 8 z(nnc,inner)
     	real* 8 v(nnc,inner+1)
	real* 8 zg(nnc), avg(nnc), sm(nnc)
	real* 8 vloc(nn_loc),avloc(nn_loc)
      logical assemble
	real* 8 eps, rnorm, rnorm0, order
	real* 8 gam,hsave,ysave,tmpo
	integer i,j,k,ij,jj,i1,j1,k1,l,iqc,igmres


      eps = 1.0e-12
      assemble = .true.

        do iqc=1,nnc
        sm(iqc) = 1.0
	  enddo

        if(iscaling.eq.0) then
	  do iqc=1,nnc
        w(iqc) = 1.0
	  enddo
        else
        do iqc=1,nnc
        if (w(iqc).lt.0.0) sm(iqc) = -1.0
        w(iqc) = 1.0 / sqrt(abs(w(iqc)))
        enddo
        endif

c	clear arrays
        call fclear (v,nnc*(inner+1))
	  call fclear (vloc,nn_loc)

c     compute residual as r = W**(-1/2) * (b - A * d)

      do iqc=1,nnc 
	v(iqc,1) = bg(iqc)
	v(iqc,1) = w(iqc) * v(iqc,1)
	enddo
	call getnorm(v(1,1),v(1,1),nnc,rnorm0)
	rnorm0 = sqrt(rnorm0)

	if(rnorm0.le.eps) then
	if(myid.eq.0) write(6,102)
	if(myid.eq.0) write(7,102)
	return
	endif


c	outer GMRES loop (igmres)
	igmres = 0 
 10	igmres = igmres + 1

c     convergence check
      call getnorm(v(1,1),v(1,1),nnc,rnorm)
	rnorm = sqrt(rnorm)
	if (rnorm.le.eps.or.igmres.gt.outer) goto 700

c	first krylov vector
      do iqc=1,nnc
	v(iqc,1) = v(iqc,1) / rnorm
	enddo

c	construct krylov vectors 2:inner and hessenberg matrix
	do j=1,inner

c	insert contribution of preconditioner
            do iqc=1,nnc
		zg(iqc) = sm(iqc) * v(iqc,j)
		z(iqc,j) = zg(iqc)
		enddo

c     compute A * v_j
            do iqc=1,nnc
	      zg(iqc) = w(iqc) * zg(iqc)
		enddo
            call gather (vloc, zg,  1, hn, hm)
		call fclear (avloc,nn_loc)
            call blockgmresf(x,u,vloc,avloc,hg,ien)
	      call scatter(avloc, avg,  1, assemble, hn, hm)
            call setid(avg,id, 1)

            do iqc=1,nnc
		avg(iqc) = w(iqc) * avg(iqc)
		v(iqc,j+1) = avg(iqc)
		enddo
		do i=1,j
                call getnorm(v(1,j+1),v(1,i),nnc,tmpo)
                h(i,j) = tmpo
		    do iqc=1,nnc
		      v(iqc,j+1) = v(iqc,j+1) - tmpo * v(iqc,i) 
		    enddo
		enddo

                call getnorm(v(1,j+1),v(1,j+1),nnc,tmpo)
		tmpo = sqrt(tmpo)
		h(j+1,j) = tmpo
		do iqc=1,nnc
		v(iqc,j+1) = v(iqc,j+1) / tmpo
		enddo

	enddo				! end inner loop

c	compute y(1:inner+1) from local hessenberg linear system
c	H_m * y = beta * e_1

c		initialize reduced residual
		y(1) = rnorm
		do i=2,inner+1
		y(i) = 0.0
		enddo

c		rotations
		do j=1,inner

			j1 = j + 1

c			previously computed rotations on column j
			do i=2,j
				i1 = i - 1
				hsave = h(i1,j)
				h(i1,j) = + cc(i1) * hsave + ss(i1) * h(i,j)
				h(i ,j) = - ss(i1) * hsave + cc(i1) * h(i,j)
			end do

c			new rotation on column j
			gam = sqrt(h(j,j)**2 + h(j1,j)**2)
			cc(j) = h(j,j) / gam
			ss(j) = h(j1,j) / gam
			h(j,j) = cc(j) * h(j,j) + ss(j) * h(j1,j)
c			note: under-diagonal term h(j+1,j) becomes 0

			y(j1) = - ss(j) * y(j)
			y(j ) = + cc(j) * y(j)
c			note: till now y(j+1) = 0
			rnorm = abs(y(j1))
		end do


c		back substitution
		j = inner !should reach here straight from rotation loop
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

c	compute dg iterate dg = dg + Z_m * y
	j = inner !(PVM only fix)
	do jj=1,j
		tmpo = y(jj)
		do iqc=1,nnc
		dg(iqc) = dg(iqc) + tmpo * z(iqc,jj)
		enddo
	end do

c	if not done recover residual vector
	if (igmres.le.outer) then
			do jj=1,j
				ij = j - jj + 2
				y(ij-1) = -ss(ij-1)*y(ij)
				y(ij) = cc(ij-1)*y(ij)
			end do
		do jj=1,j+1
			tmpo = y(jj)
			if (jj.eq.1) tmpo = tmpo - 1.0
			do iqc=1,nnc
			v(iqc,1) = v(iqc,1) + tmpo*v(iqc,jj)
			enddo
		end do
	end if

	goto 10

700   continue

c	go back to unscaled system
      do iqc=1,nnc
	dg(iqc) = w(iqc) * dg(iqc)
	enddo

	order = 0.43429448 * log(rnorm0/rnorm)
	if(myid.eq.0) write(6,101) order
	if(myid.eq.0) write(7,101) order

 101  format('Interface: convergence order = ', f5.2)
 102  format('Interface: zero residual')

	return
	end
