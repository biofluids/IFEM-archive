subroutine quad2d3n(igau,ngau,gau,wei,maxnsd,maxnquad)

	real*8  gau(maxnsd,maxnquad), wei(maxnquad)

 
	DO I = 1,maxnquad
		WEI(I) = 0.0
		DO J = 1,maxnsd
			GAU(J,I) = 0.0
		ENDDO
	ENDDO
!  
	IF(IGAU.EQ.1) THEN
		!     1-POINT GAUSS
		NGAU = 1
		GAU(1,1)=0.333333333333333
		GAU(2,1)=0.333333333333333
		WEI(  1)=0.5


		RETURN
 
	elseif (IGAU.EQ.2) then
		!     3-POINT GAUSS
		NGAU = 3
		GAU(1,1)= 0.166666666666667
        GAU(2,1)= 0.166666666666667
        GAU(1,2)= 0.66666666666667
        GAU(2,2)= 0.166666666666667
        GAU(1,3)= 0.166666666666667
        GAU(2,3)= 0.666666666666667
        WEI(1)  = 0.166666666666667
        WEI(2)  = 0.166666666666667
        WEI(3) = 0.166666666666667

        RETURN
	endif
	
END
 
