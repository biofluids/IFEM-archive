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
		GAU(3,1) = -1.0
		GAU(2,2) = -1.0
		GAU(1,3) = 1.0
		GAU(2,4) = 1.0
		GAU(1,5) = -1.0
		GAU(3,6) = 1.0

		WEI(1:6) = 4.0
		RETURN  
	  ENDIF
	  
      IF(IGAU.EQ.2) THEN 
		NGAU = 4
		GC      =+.577350269189626

        GAU(1,1)=-GC
        GAU(2,1)=-GC
        GAU(3,1)=-1.0
        GAU(1,2)=+GC
        GAU(2,2)=-GC
        GAU(3,2)=-1.0
        GAU(1,3)=-GC
        GAU(2,3)=+GC
        GAU(3,3)=-1.0
        GAU(1,4)=+GC
        GAU(2,4)=+GC
        GAU(3,4)=-1.0
		WEI(1)= 1.0
		WEI(2)= 1.0
		WEI(3)= 1.0
		WEI(4)= 1.0

c        GAU(1,5)=-GC
c        GAU(2,5)=-1.0
c        GAU(3,5)=-GC
c        GAU(1,6)=+GC
c        GAU(2,6)=-1.0
c        GAU(3,6)=-GC
c        GAU(1,7)=-GC
c        GAU(2,7)=-1.0
c        GAU(3,7)=+GC
c        GAU(1,8)=+GC
c        GAU(2,8)=-1.0
c        GAU(3,8)=+GC
c		WEI(9)= 1.0
c		WEI(10)= 1.0
c		WEI(11)= 1.0
c		WEI(12)= 1.0
c
c        GAU(1,9)=1.0
c        GAU(2,9)=-GC
c        GAU(3,9)=-GC
c        GAU(1,10)=1.0
c        GAU(2,10)=+GC
c        GAU(3,10)=-GC
c        GAU(1,11)=1.0
c        GAU(2,11)=-GC
c        GAU(3,11)=+GC
c        GAU(1,12)=1.0
c        GAU(2,12)=+GC
c        GAU(3,12)=+GC
c		WEI(9)= 1.0
c		WEI(10)= 1.0
c		WEI(11)= 1.0
c		WEI(12)= 1.0
c
c        GAU(1,13)=-GC
c        GAU(2,13)=+1.0
c        GAU(3,13)=-GC
c        GAU(1,14)=+GC
c        GAU(2,14)=+1.0
c        GAU(3,14)=-GC
c        GAU(1,15)=-GC
c        GAU(2,15)=+1.0
c        GAU(3,15)=+GC
c        GAU(1,16)=+GC
c        GAU(2,16)=+1.0
c        GAU(3,16)=+GC
c		WEI(13)= 1.0
c		WEI(14)= 1.0
c		WEI(15)= 1.0
c		WEI(16)= 1.0
c
c        GAU(1,17)=-1.0
c        GAU(2,17)=-GC
c        GAU(3,17)=-GC
c        GAU(1,18)=-1.0
c        GAU(2,18)=+GC
c        GAU(3,18)=-GC
c        GAU(1,19)=-1.0
c        GAU(2,19)=-GC
c        GAU(3,19)=+GC
c        GAU(1,20)=-1.0
c        GAU(2,20)=+GC
c        GAU(3,20)=+GC
c		WEI(17)= 1.0
c		WEI(18)= 1.0
c		WEI(19)= 1.0
c		WEI(20)= 1.0
c
c        GAU(1,22)=-GC
c        GAU(2,21)=-GC
c        GAU(3,21)=1.0
c        GAU(1,22)=+GC
c        GAU(2,22)=-GC
c        GAU(3,22)=1.0
c        GAU(1,23)=-GC
c        GAU(2,23)=+GC
c        GAU(3,23)=1.0
c        GAU(1,24)=+GC
c        GAU(2,24)=+GC
c        GAU(3,24)=1.0
c		WEI(21)= 1.0
c		WEI(22)= 1.0
c		WEI(23)= 1.0
c		WEI(24)= 1.0

		RETURN
      ENDIF

      END
 
