	subroutine initialize

	include "global.h"

      call facemap

        alpha = 1.0
           tt = 0.0
           dt = 1.0
      t_start = 0.0
      ref_lgt = 1.0
      ref_vel = 1.0
      ref_den = 1.0
      vis_liq = 1.0
      den_liq = 1.0
      vis_gas = 1.0
      den_gas = 1.0

      do isd = 1,nsd
      gravity(isd) = 0.0
      interface(isd) = -999.0
      enddo
           
      turb_kappa = 0.41

      hydro = 0 
      do i=0,21
      surf(i) = 0
      enddo

       ndf = 4
       nsd = 3
        nn = 0
        ne = 0
        nq = 0
       nqf = 0
       nec = 0
      nrng = 0

        iquad = 2
          nts = 1
          nit = 1
      ntsbout = 1
        idisk = 0

             inner = 1
             outer = 1
          iscaling = 1

         restart  = .false.
          stokes  = .false.
          steady  = .false.
          hg_vol  = .false.
          static  = .false.
	     taudt  = .false.
	   conserve = .false.

        do j=1,21
        do i=1,ndfpad
        bc(i,j) = 0
        bv(i,j) = -999.0
        enddo
          bcf(j) = 0
          bvf(j) = -999.0
        enddo

          do i=1,ndfpad
          ic(i) = -999.0
          enddo
          icf = -999

        do idel=0,21
        delta(idel) = 0.0
        end do

	do i=15
	   totaltime(i)=0.0
	   starttime(i)=0.0
	enddo

        mass = 1.0

        return
        end

      subroutine facemap
	include "global.h"

	data map /
c	triangle         |                 |                 |                 |
     & 1, 2, 3, 0, 0, 0, 2, 3, 1, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c	quadrilateral    |                 |                 |                 |
     & 1, 2, 3, 4, 0, 0, 2, 3, 4, 1, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c	tetrahedron      |                 |                 |                 |
     & 3, 1, 2, 3, 0, 0, 2, 2, 3, 1, 0, 0,
     & 1, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c	hexahedron       |                 |                 |                 |
     & 1, 1, 2, 3, 4, 5, 4, 2, 3, 4, 1, 6,
     & 3, 6, 7, 8, 5, 7, 2, 5, 6, 7, 8, 8,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c	space-time triangle                |                 |                 |
     & 1, 2, 3, 0, 0, 0, 2, 3, 1, 0, 0, 0,
     & 5, 6, 4, 0, 0, 0, 4, 5, 6, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c	space-time quadrilateral           |                 |                 |
     & 1, 2, 3, 4, 0, 0, 2, 3, 4, 1, 0, 0,
     & 6, 7, 8, 5, 0, 0, 5, 6, 7, 8, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c	space-time tetrahedron             |                 |                 |
     & 3, 1, 2, 3, 0, 0, 2, 2, 3, 1, 0, 0,
     & 1, 4, 4, 4, 0, 0, 7, 5, 6, 7, 0, 0,
     & 6, 6, 7, 5, 0, 0, 5, 8, 8, 8, 0, 0,
     & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
c	space-time hexahedron              |                 |                 |
     & 1, 1, 2, 3, 4, 5, 4, 2, 3, 4, 1, 6,
     & 3, 6, 7, 8, 5, 7, 2, 5, 6, 7, 8, 8,
     & 9, 9,10,11,12,13,12,10,11,12, 9,14,
     &11,14,15,16,13,15,10,13,14,15,16,16/

	return
	end
