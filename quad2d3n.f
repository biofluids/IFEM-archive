C
	subroutine quad2d3n(igau,ngau,gau,wei,maxnsd,maxnquad)
C
      real*8  gau(maxnsd,maxnquad), wei(maxnquad)
C
      GO TO (1,2),IGAU
C
C     1-POINT GAUSS
C
 1    CONTINUE
      NGAU = 1
      GAU(1,1)=0.333333333333333
      GAU(2,1)=0.333333333333333
      WEI(  1)=0.5
 
      RETURN
 
C     2-POINT GAUSS
 
 2    CONTINUE
      NGAU = 3
      GAU(1,1)= 0.166666666666667
      GAU(2,1)= 0.166666666666667
      GAU(1,2)= 0.166666666666667
      GAU(2,2)= 0.666666666666667
      GAU(1,3)= 0.666666666666667
      GAU(2,3)= 0.166666666666667

      WEI(1)  = 0.166666666666667
      WEI(2)  = 0.166666666666667
      WEI(3)  = 0.166666666666667
 
      RETURN
      END
 
