!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	quad3d8n.fcm                                                         c
!	---------------------------------------------------------------------c
!     quadrature definitions in three dimensional quad reference domain    c
!	---------------------------------------------------------------------c
!	920422 - written                                                     c
!	930311 - removed initialization of xquad                             c
!	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine quad3d8n(iquad, nquad, xquad, wquad, maxnsd, maxnquad)

	real* 8 xquad(maxnsd,maxnquad)
	real* 8 wquad(maxnquad)

      logical     twod
      common /twodimension/ twod

	if ((iquad.eq.1).or.(iquad.eq.111111)) then
!		1 x 1 x 1 Gaussian quadrature - code 1.1,1.1,1.1 - old code 1
		nquad = 1

		xquad(1,1) = + 0.000000000000000
		xquad(2,1) = + 0.000000000000000
		xquad(3,1) = + 0.000000000000000

		wquad(1) = 8.0
		return

	else if ((iquad.eq.2).or.(iquad.eq.121212)) then
!		2 x 2 x 2 Gaussian quadrature - code 1.2,1.2,1.2 - old code 2
		nquad = 8
		xquad(1,1) = - 0.577350269189626
		xquad(2,1) = - 0.577350269189626
		xquad(3,1) = - 0.577350269189626
		xquad(1,2) = + 0.577350269189626
		xquad(2,2) = - 0.577350269189626
		xquad(3,2) = - 0.577350269189626
		xquad(1,3) = - 0.577350269189626
		xquad(2,3) = + 0.577350269189626
		xquad(3,3) = - 0.577350269189626
		xquad(1,4) = + 0.577350269189626
		xquad(2,4) = + 0.577350269189626
		xquad(3,4) = - 0.577350269189626
		xquad(1,5) = - 0.577350269189626
		xquad(2,5) = - 0.577350269189626
		xquad(3,5) = + 0.577350269189626
		xquad(1,6) = + 0.577350269189626
		xquad(2,6) = - 0.577350269189626
		xquad(3,6) = + 0.577350269189626
		xquad(1,7) = - 0.577350269189626
		xquad(2,7) = + 0.577350269189626
		xquad(3,7) = + 0.577350269189626
		xquad(1,8) = + 0.577350269189626
		xquad(2,8) = + 0.577350269189626
		xquad(3,8) = + 0.577350269189626

		wquad(1) = 1.0
		wquad(2) = 1.0
		wquad(3) = 1.0
		wquad(4) = 1.0
		wquad(5) = 1.0
		wquad(6) = 1.0
		wquad(7) = 1.0
		wquad(8) = 1.0
		return
	else
		call error("quad3d8n: unknown quadrature code ", iquad, .true.)
	end if

	end
