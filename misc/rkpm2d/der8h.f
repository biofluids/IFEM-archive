      SUBROUTINE  DER8H (R,S,T,DH)
C     * * * * * * * * * * * * * * * * * * * * * * * * * *
C     LOCAL DERIVATIVES FOR EIGHT NODE HEXAHEDRON
C     * * * * * * * * * * * * * * * * * * * * * * * * * *
CDP   IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  DH(3,8)
C     R,S,T = LOCAL COORDINATES OF THE POINT
C     DH(1,K)=DH/DR, DH(2,K)=DH/DS, DH(3,K)=DH/DT
C     H = ELEMENT SHAPE FUNCTIONS, SEE SHP8H
      RP = 1. + R
      RM = 1. - R
      SP = 1. + S
      SM = 1. - S
      TP = 1. + T
      TM = 1. - T
      DH(1,1) =  0.125*SP*TP
      DH(2,1) =  0.125*RP*TP
      DH(3,1) =  0.125*RP*SP
      DH(1,2) =  0.125*SM*TP
      DH(2,2) = -0.125*RP*TP
      DH(3,2) =  0.125*RP*SM
      DH(1,3) =  0.125*SM*TM
      DH(2,3) = -0.125*RP*TM
      DH(3,3) = -0.125*RP*SM
      DH(1,4) =  0.125*SP*TM
      DH(2,4) =  0.125*RP*TM
      DH(3,4) = -0.125*RP*SP
      DH(1,5) = -0.125*SP*TP
      DH(2,5) =  0.125*RM*TP
      DH(3,5) =  0.125*RM*SP
      DH(1,6) = -0.125*SM*TP
      DH(2,6) = -0.125*RM*TP
      DH(3,6) =  0.125*RM*SM
      DH(1,7) = -0.125*SM*TM
      DH(2,7) = -0.125*RM*TM
      DH(3,7) = -0.125*RM*SM
      DH(1,8) = -0.125*SP*TM
      DH(2,8) =  0.125*RM*TM
      DH(3,8) = -0.125*RM*SP
      RETURN
      END
