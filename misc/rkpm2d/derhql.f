      SUBROUTINE  DERHQL (NODEDG, LOCATE, NEDGE, LEDGES, NSPACE,
     &                    RST, DERIV)
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C        SHAPE FUNCTION DERIVATIVES FOR GENERAL SERENDIPITY 
C               LINE, QUAD, OR OR HEXAHEDRON WITH AN 
C              ARBITRARY NUMBER OF NODES ON EACH EDGE
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
CDP   IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER  ( MAXDEG = 20 )
      DIMENSION  RST(3),     BLKCRD(3,8), POLI2(3),  CDRFN(3), 
     &           FARSID(3),  DERIV(3),    CRDEDG(3,MAXDEG+1), 
     &           NODEDG(12), NEATC(3,8),  NODEOP(2,12), 
     &           NODATC(3),  LOCAL(12)
      DATA BLKCRD
     &/ -1.,-1.,-1.,   1.,-1.,-1.,   1.,1.,-1.,   -1.,1.,-1.,
     &  -1.,-1., 1.,   1.,-1., 1.,   1.,1., 1.,   -1.,1., 1./
      DATA NEATC /  1,4,9,   1,2,10,  3,2,11,  3,4,12,
     &              5,8,9,   5,6,10,  7,6,11,  7,8,12 /
      DATA NODEOP / 7,8,  8,5,  5,6,   6,7,
     &              3,4,  4,1,  1,2,   2,3,
     &              3,7,  4,8,  1,5,   2,6 /
      DATA LOCAL / -1, -2, 1, 2, -1, -2, 1, 2, -3, -3, -3, -3 /
C     BLKCRD  = BLOCK CORNER LOCAL COORDINATES
C     CRDEDG  = LOCAL COORDINATES OF SIDE NODES JOINING CORNER
C     FARSID  = FAR SIDE LOCAL COORDINATE
C     LEDGES  = NUMBER OF ELEMENT EDGES, 1, 4, OR 12
C     LOCAL   = LOCAL COORDINATE PARALLEL TO EACH EDGE
C     LOCATE  = POSITION NUMBER ON EDGE, 0 IF CORNER
C     MAXDEG  = MAXIMUM PLOYNOMIAL DEGREE ON ANY SIDE
C     NEATC   = THE 1, 2, OR 3 EDGES AT A CORNER
C     NEDGE   = EDGE NUMBER OR CORNER NUMBER OF THE NODE COMPUTED
C     NODATC  = NUMBER OF SIDE NODES JOINING A CORNER
C     NODEDG  = NUMBER ON NODES ON 1,4, OR 12 EDGES
C     NODEOP  = 2 DIAGONALLY OPPOSITE NODES FOR EACH EDGE
C     NSPACE  = NUMBER OF SPATIAL DIMENSIONS
C     RST     = LOCAL COORDINATES FOR EVALUATION
C     VALUE   = SHAPE FUNCTION VALUE (RETURNED)
C
C     VALUE = A(R,S,T)*( P1(R) + P2(S) + P3(T) + CONSTANT )
C     DERIV = DA(R,S,T)*( P1(R) + P2(S) + P3(T) + CONSTANT )
C           + A(R,S,T)*( DP1(R) + DP2(S) + DP3(T) )
C
C REF: G. ZAVARISE, ET AL, "AN ALGORITHM FOR GENERATION OF SHAPE
C      FUNCTIONS IN SERENDIPITY ELEMENTS", ENG COMP,8,19-31,1991
C
C   T:  S      C8 *---E7----* C7     T:  S         8---15----7   
C    : /         /.        /:         : /         /.        /:
C    :/         / .       / :         :/         / 22      / :
C    *---R     / E12     /  E11       *---R     /  .      /  20 
C            E8   .     E6  :                 16   21    /   :
C            /    .    /    :                 /    .    /    :
C           /   C4*.../.E3..* C3             /     4.13/.12..3   
C          /     .   /     /                /     .   /     /
C      C5 *--E5-----* C6  /                5---------6     /
C         :    .    :    /                 :    .    :    11
C         :  E4     :   E2                 :  14    19   /   
C        E9  .    E10  /                  17  .      :  10
C         : .       : /                    : .      18 /
C         :.        :/                     :.        :/
C      C1 *---E1----* C2                   1----9----2
C CORNER NODE & EDGE NUMBERS.   22 NODES: CORNERS, THEN BY EDGES. 
C                               CCW IF |T|=1, ELSE IN POSITIVE T.
C                       === 3-D FORM ===
C
C      C4 *---E3----* C3                   4----8----3     
C         :    .    :         :S           :         :   
C         :         :         :            :         :
C        E4        E2         :            9         :  
C         :         :         *---R        :         :
C         :         :                      :         :
C      C1 *---E1----* C2                   1--5-6-7--2 
C CORNER NODE & EDGE NUMBERS.  9 NODES: CORNERS, THEN BY EDGE ORDER.
C                        === 2-D FORM ===
C
C      C1 *---E1----* C2                   1--2-3-4--5 
C CORNER NODE & EDGE NUMBERS.       9 NODES NUMBERED BY EDGE ORDER.
C                        === 1-D FORM ===
      POLI1 = 1.
      IF ( LOCATE .EQ. 0 ) THEN
C
C  SHAPE FUNCTION FOR CORNER NODES
C
        DO 100 ICORD = 1,NSPACE
          POLI1 = POLI1*(RST(ICORD) + BLKCRD(ICORD,NEDGE))
     &            /(2*BLKCRD(ICORD,NEDGE))
  100   CONTINUE
        CPNUL = 1.
        POLI2(1) = 0.
        POLI2(2) = 0.
        POLI2(3) = 0.
        DO 200 ICORD = 1,NSPACE
          NSIDE = NEATC(ICORD,NEDGE)
          NODATC(ICORD) = NODEDG(NSIDE) - 2
          IF ( NODATC(ICORD) .GT. 0 ) THEN
            IF ( NODATC(ICORD) .GT. MAXDEG ) STOP 'MAXDEG, DERSHAFN'
            CPNUL = CPNUL - 1.
            POLI2(ICORD) = 1.
            FARSID(ICORD) = 2./(NODEDG(NSIDE) - 1)
            DO 300 INODE = 1,NODATC(ICORD)
              CRDEDG(ICORD,INODE) =  -1. + FARSID(ICORD)*INODE
              POLI2(ICORD) = POLI2(ICORD)*(RST(ICORD) 
     &                     - CRDEDG(ICORD,INODE))/(BLKCRD(ICORD,NEDGE)
     &                     - CRDEDG(ICORD,INODE))
  300       CONTINUE
          ENDIF
  200   CONTINUE
C       VALUE = POLI1*(POLI2(1) + POLI2(2) + POLI2(3) + CPNUL)
      ELSE
C
C  SHAPE FUNCTION FOR EDGE NODES
C
        NOPV1 = NODEOP(1,NEDGE)
        NOPV2 = NODEOP(2,NEDGE)
        ISRFN = ABS(LOCAL(NEDGE))
        FARSID(1) = 2./(NODEDG(NEDGE) - 1)
        CDRFN(1) =  -BLKCRD(1,NOPV1)
        CDRFN(2) =  -BLKCRD(2,NOPV1)
        CDRFN(3) =  -BLKCRD(3,NOPV1)
        CDRFN(ISRFN) = (1. - FARSID(1)*LOCATE)*LOCAL(NEDGE)/ISRFN
        DO 400 ICORD = 1,NSPACE
          POLI1 = POLI1*(RST(ICORD) - BLKCRD(ICORD,NOPV1))
     &           /(CDRFN(ICORD) - BLKCRD(ICORD,NOPV1))
  400   CONTINUE
        PLAN2 = (RST(ISRFN) - BLKCRD(ISRFN,NOPV2))
     &         /(CDRFN(ISRFN) - BLKCRD(ISRFN,NOPV2))
        POLI3 = 1.
        NODATC(1) = NODEDG(NEDGE) - 2
        IF ( NODATC(1) .GT. 0 ) THEN
          IF ( NODATC(1) .GT. MAXDEG ) STOP 'MAXDEG, DERSHAFN'
          DO 500 INODE = 1,NODATC(1)
            CRDEDG(1,INODE) =  -1. + FARSID(1)*INODE
            IF ( ABS(CRDEDG(1,INODE) - CDRFN(ISRFN)) .GT. 0.0001)
     &      THEN
              POLI3 = POLI3*(RST(ISRFN) - CRDEDG(1,INODE))
     &               /(CDRFN(ISRFN) - CRDEDG(1,INODE))
            ENDIF
  500     CONTINUE
        ENDIF
C       VALUE = POLI1*PLAN2*POLI3
      ENDIF
C
C  DERIVATIVES OF SHAPE FUNCTIONS
C
      DO 600 ICOR1 = 1,NSPACE
        IF ( LOCATE .EQ. 0 ) THEN
C
C  DERIVATIVES FOR CORNER NODES
C
          DPOL1 = POLI2(1) + POLI2(2) + POLI2(3) + CPNUL
          DO 700 ICOR2 = 1,NSPACE
            IF ( ICOR2 .NE. ICOR1 ) THEN
              DPOL1 = DPOL1*(RST(ICOR2) + BLKCRD(ICOR2,NEDGE))
     &               /(2*BLKCRD(ICOR2,NEDGE))
            ELSE
              DPOL1 = DPOL1/(2*BLKCRD(ICOR2,NEDGE))
            ENDIF
  700     CONTINUE
          DPOL2 = 0.
          DO 800 INOD1 = 1,NODATC(ICOR1)
            DETP2 = 1.
            DO 900 INOD2 = 1,NODATC(ICOR1)
              IF ( INOD2 .NE. INOD1 ) THEN
                DETP2 = DETP2*(RST(ICOR1) - CRDEDG(ICOR1,INOD2))
     &                 /(BLKCRD(ICOR1,NEDGE) - CRDEDG(ICOR1,INOD2))
              ELSE
                DETP2 = DETP2/(BLKCRD(ICOR1,NEDGE)
     &                - CRDEDG(ICOR1,INOD2))
              ENDIF
  900       CONTINUE
            DPOL2 = DPOL2 + DETP2
  800     CONTINUE
          DPOL2 = DPOL2*POLI1
          DERIV(ICOR1) = DPOL1 + DPOL2
        ELSE
C
C  DERIVATIVES FOR EDGE NODES
C
          DPOL1 = POLI3*PLAN2
          DO 1000 ICOR2 = 1,NSPACE
            IF ( ICOR2 .NE. ICOR1 ) THEN
              DPOL1 = DPOL1*(RST(ICOR2) - BLKCRD(ICOR2,NOPV1))
     &               /(CDRFN(ICOR2) - BLKCRD(ICOR2,NOPV1))
            ELSE
              DPOL1 = DPOL1/(CDRFN(ICOR2) - BLKCRD(ICOR2,NOPV1))
            ENDIF
 1000     CONTINUE
          DPLA2 = 0.
          DPOL3 = 0.
          IF ( ICOR1 .EQ. ISRFN ) THEN
            DPLA2 = POLI1*POLI3/(CDRFN(ISRFN) - BLKCRD(ISRFN,NOPV2))
            DO 1100 INOD1 = 1,NODATC(1)
              IF ( ABS(CRDEDG(1,INOD1) - CDRFN(ISRFN)) .GT. 0.0001)
     &        THEN
                DETP3 = 1.
                DO 1200 INOD2 = 1,NODATC(1)
                  IF ( ABS(CRDEDG(1,INOD2) - CDRFN(ISRFN)) .GT.
     &                 0.0001 ) THEN
                    IF ( INOD2 .NE. INOD1 ) THEN
                      DETP3 = DETP3*(RST(ISRFN) - CRDEDG(1,INOD2))
     &                        /(CDRFN(ISRFN) - CRDEDG(1,INOD2))
                    ELSE
                      DETP3 = DETP3/(CDRFN(ISRFN) - CRDEDG(1,INOD2))
                    ENDIF
                  ENDIF
 1200           CONTINUE
                DPOL3 = DPOL3 + DETP3
              ENDIF
 1100       CONTINUE
            DPOL3 = DPOL3*POLI1*PLAN2
          ENDIF
          DERIV(ICOR1) = DPOL1 + DPLA2 + DPOL3
	ENDIF
  600   CONTINUE
      RETURN
      END
