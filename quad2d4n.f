C
	subroutine quad2d4n(igau,ngau,gau,wei,maxnsd,maxnquad)
C
      real*8 gau(maxnsd,maxnquad), wei(maxnquad)
C
	DO I = 1,maxnquad
	WEI(I) = 0.0
	DO J = 1,maxnsd
	GAU(J,I) = 0.0
	ENDDO
	ENDDO
C
      IF(IGAU.EQ.1) THEN
	NGAU = 1  
	WEI(1) = 4.0
	RETURN  
	ENDIF
 
      IF(IGAU.EQ.2) THEN 
      NGAU = 4
      GC      =+.577350269189626
      GAU(1,1)=-GC
      GAU(2,1)=-GC
      GAU(1,2)=+GC
      GAU(2,2)=-GC
      GAU(1,3)=-GC
      GAU(2,3)=+GC
      GAU(1,4)=+GC
      GAU(2,4)=+GC
      WEI(1)= 1.0
      WEI(2)= 1.0
      WEI(3)= 1.0
      WEI(4)= 1.0
      RETURN
      ENDIF

      END
 
