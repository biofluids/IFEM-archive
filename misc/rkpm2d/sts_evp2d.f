c
      subroutine sts_evp2d(disp,vel,
     &            sts_KH,sts_PK1,
     &            eff_ste,eff_sts,dt,
     &            theta,mgk,ik,
     &            xkc1,xkc2)
c
c-----------------------------------------------------------
c
c    Documentation:
c
c    This subroutine is designated for integration of
c    elasto-viscoplastic continuum.
c    
c    Instead of adopting radial return, the so-called
c    tangent modulus method proposed by Peirce et. al.
c    is adopted.
c
c    References:
c    
c    1. Peirce, D., Shih, C. F. and Needleman, A. (1984)
c       Computers and Structures, Vol. 18, pp. 875-887
c
c    2. Needleman, A. (1989) 
c       Journal of Applied Mechanics, Vol. 56, pp. 1-9
c
c    3. Nacar, A., Needleman, A. and Ortiz, M. (1989)
c       Computer Methods in Applied Mechanics and Engineering,
c       Vol. 73, 235-258
c
c    4. Belytschko, T. and Tabbara, M. (1993)
c       International Journal for Numerical Methods in Engineering,
c       Vol. 36, 4245-4265
c    
c    
c    Tasks:
c
c --- Update the stress status of each Gauss point 
c     for large deformation rate-dependent
c     elasto-viscoplastic problem (tangent modulus method)
c
c
c    Note:
c
c    (1)  There are two ways to integrate the consititutive laws:
c         A. via Kirchhoff stress update; (Nacar et al. [1989])
c         B. via Cauchy stress update;    (Belytschko et al.[1993])
c
c    (2)  The tangent modulus method proposed by Peirce et al.
c         is a forward gradient method.
c         It may not have the capacity to accommodate any
c         predictor/corrector algorithm without modification.
c     
c
c    Variables:
c
c --- Input parameters:
c
c    Lgp :
c    Lmgp:
c    dt  : the time increament;
c    veli: the velocity array at   T = dt * n  
c          for velocity other than T = dt * n, 
c          there is an inconsistency with the tangent
c          modulus algorithm used.
c
c    sts_PK1: of LAST time step at this Gauss point 
c             ( The First Piola-Kirchhoff );
c    (xkc1,xkc2): the coordinate of the Gauss point;
c
c    eff_ste:  effective strain of LAST time step at this Gauss point;
c    eff_sts:  effective stress of LAST time step at this Gauss point;
c    mgk    :  total number of Gauss points;
c    ik     :  the # of current gauss point;
c    maxNumnp: the upper limit of maximum nodal points;
c    maxGP:    the upper limit of maximum Gauss points;
c
c
c --- Local parameters:
c
c    sts_KH: of LAST time step at this Gauss point (Kirchhoff stress);
c    dmat: an array that contants material constants;
c
c
c --- output parameters:
c
c    sts_PK1: of CURRENT time step at this Gauss point 
c             ( The First Piola-Kirchhoff )
c
c    eff_sts: of CURRENT time step
c    eff_ste: of CURRENT time step
c
c--------------------------------------------------------------
c
      implicit none
      include 'parameter.h'
c                                                 c
c............. I&O parameters ....................c
c                                                 c
      integer ik,mgk
c
      real*8 disp(2,maxNumnp),vel(2,maxNumnp)
      real*8 sts_PK1(3,3),sts_KH(3,3),
     &       xkc1,xkc2,theta
c
      real*8 dt,eff_sts, eff_ste
c      
c............. Local arrays ...............
c
      real*8 Dtot(4),Wt12,               !
     &       Dsts_KH(4),                 ! Devitoric part of stress;
     &       dvD(2,2),dmat(4,4),
     &       smlp(4),PP(4), AtA(4,4),EP_tan(4),
     &       sig_jrate(4),dLtan(4,4),sts_dot(4) 
c
c......... Local variables ..................
c
      real*8 F_matx(3,3),Finv_matx(3,3),detF
      real*8 trace_KH
      real*8 Sb,Eb                ! (Sb := eff_sts_n; Eb := eff_ste_n)
      real*8 Ed,Es,G0,C1111,C1122,C2211,
     &       C2222,C3333,C1212,C3311,C3322
      real*8 Eb1,gtem,             ! (gtem := 1 + (Eb/Eb1)^2
     &       Eb0,sig0,bigN,gsoft,  ! (gsoft:= sig0 (1+ Eb/Eb0)^N/gtem )
     &       rc,                   ! (rc = xkc1^2 + xkc2^2
     &       Ebtot_c,Ebtot,dG_Eb,  ! (dG_Eb = d gsoft / d Eb )
     &       xxii,                 ! (xxii = theta*dt*hh*dEbtot_Eb) 
     &       dEbtot_Sb,hh,         ! (hh = p:L:p - dEbtot_Eb/dEbtot_Sb)
     &       reg,                  ! (reg := Sb/gsoft)
     &       smm,sminv             ! (sminv := 1/smm)
c
c..............................................
c     
c    rc   : the distance of the Gauss point from the origin;
c    gsoft: the rate-dependent relation of eff_strain with
c           the softening effect;
c
c.........Local constant .......................
c
      real*8 eps1,
     &       Eb0tot,Sb0,ecc,
     &       pCp,ctcon,
     &       r0,xc1,xc2,zimp,cr1,cr2
c      
      integer ip,ipt,i,j,mLoop,n_three
c
      logical elastic
c
c........ common blocks ..................
c
      integer lnods(maxNode,maxElem),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        Lgp(maxGP),Lmgp(mnsch,maxGP)
c
      real*8 el_E,el_V,el_G,el_K
c
      real*8 shpk(mnsch,maxGP),shpkdx(mnsch,maxGP),
     &       shpkdy(mnsch,maxGP),shpn(mnsch,maxNumnp)
c
      common /elastic/el_E, el_V, el_G, el_K
      common /shear/r0,xc1,xc2,zimp,cr1,cr2
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /shapeK/shpk,shpkdx,shpkdy,shpn
c
c------------------------------------------------------
c------------------------------------------------------
c
c......Assign the material constants, and
c      the intrinsic constants
c
c-------------------------------------------------------
c
      n_three = 3
      eps1    = 1.0d-16    ! ( the criteria for Jacobian)
c      
c----------------------------------------------------------
c * Eb1 = 1.0d4*Eb0 (hardening), = 1.0d2 * Eb0 (softening) *
c------------------------------------------------------------
c
      Eb0     = 2.18d-3
      Eb1     = 100.0 * Eb0
      Sb0     = 460.0d6
      bigN    = 0.10d0	   !
      Eb0tot  = 2.0d-3	   ! eb0 dot			
      smm     = 0.010d0	   ! m
      sminv   = 100.0d0    ! 1./m
c
c......Implement the imperfection into the viscoplasticity
c      stress accumulate function 
c      
c
      rc   = (xc1 - xkc1)**2. + (xc2 - xkc2)**2.
      sig0 = Sb0*(1.0d0 - zimp*dexp(-rc/(r0*r0)))
c       
c........Find the number of node surrouding
c        the Gauss point, ik, i.e. mLoop
c
      mLoop = Lgp(ik)
c         
c........... Calculate the deformation gradient
c            and calculate the rate of deformation 
c            at T = dt * n ( the last time step ? )
c
         do i = 1, 2
	    do j = 1,2
	       dvD(i,j)  = 0.0d0
	    enddo
         enddo
c
	 do i = 1,4
	    do j = 1,4
	       AtA(i,j)   = 0.0d0
	       dLtan(i,j) = 0.0d0
	    enddo
	    Dtot(i)       = 0.0d0
	    sig_jrate(i)  = 0.0d0
	    EP_tan(i)     = 0.0d0
	    PP(i)         = 0.0d0
	    sts_dot(i)    = 0.0d0
	    smlp(i)       = 0.0d0
	    Dsts_KH(i)    = 0.0d0
	 enddo
c
         do i=1,3
            do j=1,3       
               F_matx(i,j)    = 0.0d0
               Finv_matx(i,j) = 0.0d0
            enddo
         enddo
c 
         do ip = 1, mLoop
            ipt = Lmgp(ip,ik)
c
c.......(1)................
c
            F_matx(1,1) = F_matx(1,1)
     &                  + shpkdx(ip,ik) * disp(1,ipt)
            F_matx(1,2) = F_matx(1,2) 
     &                  + shpkdy(ip,ik) * disp(1,ipt)
            F_matx(2,1) = F_matx(2,1) 
     &                  + shpkdx(ip,ik) * disp(2,ipt)
            F_matx(2,2) = F_matx(2,2) 
     &                  + shpkdy(ip,ik) * disp(2,ipt)
c
c.......(2)..............
c
            dvD(1,1) = dvD(1,1) + shpkdx(ip,ik)*vel(1,ipt)
            dvD(1,2) = dvD(1,2) + shpkdy(ip,ik)*vel(1,ipt)
            dvD(2,1) = dvD(2,1) + shpkdx(ip,ik)*vel(2,ipt)
            dvD(2,2) = dvD(2,2) + shpkdy(ip,ik)*vel(2,ipt)
c
         enddo
c
c.......Add the Keroneker delta part:
c
         do i=1,3
            F_matx(i,i) = F_matx(i,i) + 1.0d0
         enddo
c         
c.........Find the Jacobian
c
        detF = F_matx(1,1)*F_matx(2,2) - F_matx(1,2)*F_matx(2,1)
c
        if (detF .le. 1.0d-8) then
	   print *, 'Negative Jacobian', detF,' at Gauss pt:', ik
	   stop
        endif
c
c----------------------------------------------------------
c
c................ End of Calculation .............. 
c
c----------------------------------------------------------
c
c.....Manule Inversion of Deformation Gradient
c
      Finv_matx(1,1) =   F_matx(2,2)/detF 
      Finv_matx(1,2) = - F_matx(1,2)/detF 
      Finv_matx(2,1) = - F_matx(2,1)/detF 
      Finv_matx(2,2) =   F_matx(1,1)/detF 
      Finv_matx(3,3) =   1.0d0
c
c
c.....Rate of Deformation ( The symmetry part )
c
c     1   2   3   4
c     11  22  33  12
c
      Dtot(1) = dvD(1,1)*Finv_matx(1,1) + dvD(1,2)*Finv_matx(2,1)
      Dtot(2) = dvD(2,1)*Finv_matx(1,2) + dvD(2,2)*Finv_matx(2,2)
      Dtot(3) = 0.00d0
      Dtot(4) = dvD(1,1)*Finv_matx(1,2) + dvD(1,2)*Finv_matx(2,2)
     &        + dvD(2,1)*Finv_matx(1,1) + dvD(2,2)*Finv_matx(2,1)
c
      Dtot(4) = 0.50d0 * Dtot(4)
c
c.....Spin.........
c
      Wt12  = dvD(1,1)*Finv_matx(1,2) + dvD(1,2)*Finv_matx(2,2)
     &      - dvD(2,1)*Finv_matx(1,1) - dvD(2,2)*Finv_matx(2,1)
      Wt12  = 0.50d0 * Wt12
c
c----------------------------------------------------------
c	 End of calculation rate of deformation for large deformation
c        or small deformation
c----------------------------------------------------------
c	                        
c
c......Calculate the Elastic Muduli tensor
c
c      For plane strain
c
c      E0 = el_E/(1-el_V**2)
c      V0 = el_V/(1-el_V)
c      D0 = E0/(1-V0*V0)
c
c.....................................................
c
c    The Voigt Notation: 
c    
c    1     2     3     4     5     6
c
c    11    22    33   12     23    31
c
c    C1111 = C11 :   dmat(1,1)
c    C1122 = C12 :   dmat(1,2)
c    C2211 = C21 :   dmat(2,1)
c    C2222 = C22 :   dmat(2,2)
c    C3333 = C33 :   dmat(3,3)
c    C1212 = C44 :   dmat(4,4)
c
c    For plane strain, we take
c
c    C3333 = C1133 = C3311 = C2233 = C3322 = 0
c    C1313 = C1112 = C2212 = C3312 = 0
c      
c.......................................................
c
      Ed = el_E *(1.0 - el_V) /((1.0 + el_V)*(1.0d0 - 2.0d0 * el_V))
      Es = el_E * el_V /((1.0d0 + el_V)*(1.0d0 - 2 * el_V))
      G0 = el_E/(1.0d0 + el_V)
c	 
      C1111 = Ed
      C1122 = Es
      C2211 = Es
      C2222 = Ed
      C1212 = G0
      C3333 = 0.0  ! This doesn't matter (plane strain)
      C3311 = Es
      C3322 = Es
c
      dmat(1,1) = C1111
      dmat(1,2) = C1122
      dmat(1,3) = 0.0d0
      dmat(1,4) = 0.0d0
c
      dmat(2,1) = dmat(1,2)
      dmat(2,2) = C2222
      dmat(2,3) = 0.0d0
      dmat(2,4) = 0.0d0
c
      dmat(3,1) = C3311
      dmat(3,2) = C3322
      dmat(3,3) = C3333
      dmat(3,4) = 0.0d0
c
      dmat(4,1) = 0.0d0
      dmat(4,2) = 0.0d0
      dmat(4,3) = 0.0d0
      dmat(4,4) = C1212
c
c----------------------------
c
c     Analysis of Stress
c
c----------------------------
c
       trace_KH = 0.0d0
c
       do i = 1, 3
          trace_KH = trace_KH + sts_KH(i,i)
       enddo
       trace_KH = trace_KH/3.0d0
c
c------------------------------------------------
c    
c    Voigt notation
c
c    1     2     3     4     5     6
c
c    11    22    33   12     23    31
c
c------------------------------------------------
c
       Dsts_KH(1) = sts_KH(1,1) - trace_KH
       Dsts_KH(2) = sts_KH(2,2) - trace_KH
       Dsts_KH(3) = sts_KH(3,3) - trace_KH 
       Dsts_KH(4) = sts_KH(1,2) 
c
c  *NOTE*:
c  There is a very subtle point about initialing
c  Sb. Since Sb = 0.0 at beginning, it will cause
c  singularity at computing "smlp(i)"
c
c
       Sb = eff_sts
       Eb = eff_ste
c
       if (Sb .le. eps1) then
	  xxii  = 0.0d0
	  ctcon = 0.0d0
	  ecc   = 0.0d0
	  Ebtot = 0.0d0
c
	  elastic = .true.
	  go to 138
       endif
c
c------------------------------------------------
c
c.......... Calculate C:p and p:C:p .............
c
c------------------------------------------------
c
       do i = 1, 4
          smlp(i) = 3.0d0 * Dsts_KH(i)/(2.0 * Sb)
       enddo
c
      do i = 1, 4
	 PP(1) = PP(1) + dmat(1,i)*smlp(i)
	 PP(2) = PP(2) + dmat(2,i)*smlp(i)
	 PP(3) = PP(3) + dmat(3,i)*smlp(i)
	 PP(4) = PP(4) + dmat(4,i)*smlp(i)
      enddo
c
c....... pCp  Remark: should we use pCp = 3 G ?
c        in 2d, we better calulate it.
c
      pCp = 0.0
      do i = 1, 4
	 pCp = pCp + smlp(i)* PP(i)
      enddo
c
c------------------------------------------------
c
c...... Calculate the evolution function
c
c    dEbtot_Eb := d Ebtot / d Ebar;   
c    dEbtot_Sb := d Ebtot / d Sbar;
c    dG_Eb     := d gsoft / d Ebar
c
c--------------------------------------------------
c
c
      gtem  = 1.0d0 + (Eb/Eb1)**2
      gsoft = sig0 * (1.0 + Eb/Eb0)**bigN
      gsoft = gsoft/gtem
c
      reg   = Sb/gsoft
c
      Ebtot = Eb0tot * (reg)**(sminv)
c
      dG_Eb = gsoft * (
     &        bigN /(Eb + Eb0)
     &      - 2.0 * Eb/(Eb1**2 + Eb**2)
     &        ) 
c
c
      dEbtot_Sb = sminv*Ebtot/Sb
c
c......hh    = pCp - dEbtot_Eb/dEbtot_Sb
c
        hh    = pCp + reg * dG_Eb
        xxii  = theta * dt * hh * dEbtot_Sb
        ctcon = xxii/((1.0 + xxii)*hh)
c
c------------------------------------------------
c
c       Calculate the Tangent Moduli
c
c-----------------------------------------------
c
      do i = 1,4
	 AtA(1,i) = PP(1) * PP(i)
	 AtA(2,i) = PP(2) * PP(i)
	 AtA(3,i) = PP(3) * PP(i)
	 AtA(4,i) = PP(4) * PP(i)
      enddo
c
      do i = 1,	4
         EP_tan(i)  = Ebtot * PP(i) / (1.0 + xxii ) 
	 dLtan(1,i) = dmat(1,i) - ctcon * AtA(1,i)
	 dLtan(2,i) = dmat(2,i) - ctcon * AtA(2,i)
	 dLtan(3,i) = dmat(3,i) - ctcon * AtA(3,i)
	 dLtan(4,i) = dmat(4,i) - ctcon * AtA(4,i)
      enddo
c
c--------------------------------------------------------
c
c  Calculate the Jaumann rate based on consititutive law
c
c---------------------------------------------------------
c
c.........Calculate the effective plastic strain rate, and
c         and update the effective plastic strain
c
      ecc = 0.0
      do i = 1, 4
	 ecc = ecc + PP(i) * Dtot(i)
      enddo
c
 138  continue
c
      Ebtot_c = Ebtot/(1.0 + xxii) + ctcon * ecc
      eff_ste = Eb + Ebtot_c * dt 
c
c
      do i = 1,4
	 do j = 1, 4
            sig_jrate(i) = sig_jrate(i) + dLtan(i,j) * Dtot(j) 
	 enddo
      enddo
c
      if (.not. elastic) then
         do i = 1,4
	    sig_jrate(i) = sig_jrate(i) - EP_tan(i)
         enddo
      endif
c
c..........Update the Kirchoff stress.............
c
         sts_dot(1) = sig_jrate(1) + 2.0d0 * sts_KH(1,2) * Wt12 
         sts_dot(4) = sig_jrate(4) 
     &              - Wt12*(sts_KH(1,1) - sts_KH(2,2)) 
         sts_dot(2) = sig_jrate(2) - 2.0d0 * sts_KH(1,2) * Wt12 
         sts_dot(3) = sig_jrate(3)
c
         sts_KH(1,1) = sts_KH(1,1) + sts_dot(1)*dt
         sts_KH(1,2) = sts_KH(1,2) + sts_dot(4)*dt
         sts_KH(2,1) = sts_KH(1,2)
         sts_KH(2,2) = sts_KH(2,2) + sts_dot(2)*dt
         sts_KH(3,3) = sts_KH(3,3) + sts_dot(3)*dt
c
c-----------------------------------------------
c........Update the effective stress
c
c-----------------------------------------------
c
c------------------------------------------------
c
c    Voigt notation
c
c    1     2     3     4     5     6
c
c    11    22    33   12     23    31
c
c------------------------------------------------
c
c 
       Sb  = 0.5 * ((sts_KH(1,1) - sts_KH(2,2))**2
     1        +     (sts_KH(2,2) - sts_KH(3,3))**2 
     2        +     (sts_KH(3,3) - sts_KH(1,1))**2 ) 
     3        + 3.0 * sts_KH(1,2)**2
c
       eff_sts = sqrt(Sb)
c
c........Update the PK1
c
      call mul_matx(Finv_matx,sts_KH,sts_PK1,n_three,n_three,n_three)
c
      end
c
c
