c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	quad3d4n.fcm                                                         c
c	---------------------------------------------------------------------c
c     quadrature definitions in three dimensional tetrahedral ref domain   c
c	---------------------------------------------------------------------c
c	930329 - converted from quad2d3n.fcm                                 c
c	         rules taken from Hughes Table 3.I.2                         c
c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine quad3d4n(iquad, nquad, xquad, wquad, maxnsd, maxnquad)

	real* 8 xquad(maxnsd,maxnquad)
	real* 8 wquad(maxnquad)

	if ((iquad.eq.1).or.(iquad.eq.313131)) then

c	1 point quadrature - code 3.1,3.1,3.1 - old code 1 - precision 2
		nquad = 1

		xquad(1,1) = + 0.250000000000000
		xquad(2,1) = + 0.250000000000000
		xquad(3,1) = + 0.250000000000000

		wquad(1) = 0.166666666666667
		return

	else if ((iquad.eq.2).or.(iquad.eq.323232)) then
c		4 point quadrature - code 3.2,3.2,3.2 - old code 2 - precision 3
		nquad = 4

		xquad(1,1) = + 0.585410200000000
		xquad(2,1) = + 0.138196600000000
		xquad(3,1) = + 0.138196600000000
		xquad(1,2) = + 0.138196600000000
		xquad(2,2) = + 0.585410200000000
		xquad(3,2) = + 0.138196600000000
		xquad(1,3) = + 0.138196600000000
		xquad(2,3) = + 0.138196600000000
		xquad(3,3) = + 0.585410200000000
		xquad(1,4) = + 0.138196600000000
		xquad(2,4) = + 0.138196600000000
		xquad(3,4) = + 0.138196600000000

		wquad(1) = 0.041666666666667
		wquad(2) = 0.041666666666667
		wquad(3) = 0.041666666666667
		wquad(4) = 0.041666666666667
		return

	else if ((iquad.eq.3).or.(iquad.eq.333333)) then
c		5 point quadrature - code 3.3,3.3,3.3 - old code 3 - precision 4
		nquad = 5

		xquad(1,1) = + 0.250000000000000
		xquad(2,1) = + 0.250000000000000
		xquad(3,1) = + 0.250000000000000
		xquad(1,2) = + 0.333333333333333
		xquad(2,2) = + 0.166666666666667
		xquad(3,2) = + 0.166666666666667
		xquad(1,3) = + 0.166666666666667
		xquad(2,3) = + 0.333333333333333
		xquad(3,3) = + 0.166666666666667
		xquad(1,4) = + 0.166666666666667
		xquad(2,4) = + 0.166666666666667
		xquad(3,4) = + 0.333333333333333
		xquad(1,5) = + 0.166666666666667
		xquad(2,5) = + 0.166666666666667
		xquad(3,5) = + 0.166666666666667

		wquad(1) = - 0.133333333333333
		wquad(2) = + 0.075000000000000
		wquad(3) = + 0.075000000000000
		wquad(4) = + 0.075000000000000
		wquad(5) = + 0.075000000000000
		return

	else if ((iquad.eq.4).or.(iquad.eq.222222)) then
c	4 point Lobatto quadrature - element corners - code 2.2,2.2,2.2 - old 4
		nquad = 4

		xquad(1,1) = + 0.0
		xquad(2,1) = + 0.0
		xquad(3,1) = + 0.0
		xquad(1,2) = + 1.0
		xquad(2,2) = + 0.0
		xquad(3,2) = + 0.0
		xquad(1,3) = + 0.0
		xquad(2,3) = + 1.0
		xquad(3,3) = + 0.0
		xquad(1,4) = + 0.0
		xquad(2,4) = + 0.0
		xquad(3,4) = + 1.0
	
		wquad(1) = 0.041666666666667
		wquad(2) = 0.041666666666667
		wquad(3) = 0.041666666666667
		wquad(4) = 0.041666666666667
		return

	else
		call error("quad3d4n: unknown quadrature code ", iquad, .true.)
	end if

	end
