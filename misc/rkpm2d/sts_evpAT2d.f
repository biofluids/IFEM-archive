c
      subroutine sts_evpAT2d(disp,vel,
     &           sts_KH,sts_PK1,Temp,
     &           eff_ste,eff_sts,dt,
     &           theta,ik,
     &           xkc1,xkc2)
c
c-----------------------------------------------------------
c
c    Documentation: (without damage)
c
c    This subroutine is designated for integration of
c    thermo-elasto-viscoplastic continuum.
c    Because the strain rate is high, we only consider
c    the adiabatic situation, that is, the heat conduction
c    is neglected.
c
c    
c    The so-called tangent modulus method proposed 
c    by Peirce et. al.  is extended to adiabatic case.
c    A monolithic or simultaneous integration scheme
c    is implemented, which is an extension of 
c    Perirce's tangent modulus method to the adiabatic case.
c    (The algorithm has not been published yet)
c
c
c    References:
c    
c    1. Peirce, D., Shih, C. F. and Needleman, A. (1984)
c       Computers and Structures, Vol. 18, pp. 875-887
c
c    2. Needleman, A. (1989) 
c       Journal of Applied Mechanics, Vol. 56, pp. 1-9
c
c    3. Lemonds, J. & Needleman, A (1986)
c       Mechanics of Materials, Vol. 5, 339-361
c
c    4. Zhou, M. Ravichardran, G. and Rosakis, A.J. (1996)
c       Journal of Mechanics and Physics of Solids, Vol. 44, No,6
c       pp. 1007-1032
c
c    5. Shaofan Li (1999)
c       Research Note:  A Thermo-elasto-viscoplastic Solid
c                       Under Large Deformation
c       NU-MECH 99-03
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
c    (1)  Kirchhoff stress update is adopted. 
c
c
c    (2)  The tangent modulus method proposed by Peirce et al.
c         is a forward gradient method, which is modified in
c         this code to accomodate the temperature field.
c     
c
c    Variables:
c
c --- Input parameters:
c
c    Lgp:
c    Lmgp:
c    dt : the time increament;
c    vel: the velocity array at   T = dt * n  
c         for velocity other than T = dt * n, 
c         there is an inconsistency with the tangent
c         modulus algorithm used.
c
c    sts_PK1: of LAST time step at this Gauss point 
c             ( The First Piola-Kirchhoff );
c    sts_KH: of LAST time step at this Gauss point (Kirchhoff stress);
c    (xkc(3)): the coordinate of the Gauss point;
c
c    eff_ste:  effective strain of LAST time step at this Gauss point;
c    eff_sts:  effective stress of LAST time step at this Gauss point;
c    ik     :  the # of current gauss point;
c    maxNumnp: the upper limit of maximum nodal points;
c    maxGP:    the upper limit of maximum Gauss points;
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
c    sts_KH : of CURRENT time step at this Gauss pt (Cauchy stress)
c             ( for ploting)
c
c    eff_sts: of CURRENT time step
c    eff_ste: of CURRENT time step
c
c
c    Author: Shaofan Li
c    September, 1998
c 
c--------------------------------------------------------------
c
      implicit none
      include 'parameter.h'
c
c.......... I&O parameters .............
c
      integer ik
      integer Lgp(maxGP),Lmgp(mnsch,maxGP),
     &        Lnp(maxNumnp),Lmnp(mnsch,maxNumnp),
     &        lnods(maxNode,maxElem)
c
      real*8 shpk(mnsch,maxGP),
     &       shpkdx(mnsch,maxGP),
     &       shpkdy(mnsch,maxGP),
     &       shpn(mnsch,maxNumnp)
c
      real*8 disp(2,maxNumnp),vel(2,maxNumnp),
     &       sts_PK1(3,3),sts_KH(3,3)
c 
      real*8 Temp,kappa,delta,Temp_rate,
     &       alpha_T,Cp,xkc1,xkc2,theta,
     &       eff_sts,eff_ste
c      
c............. Local arrays ...............
c
      real*8 Wt12,sts_dot(4)
      real*8 Dsts_KH(4)             ! Devitoric part of velocity;
      real*8 dvD(2,2),dmat(4,4)
      real*8 smlp(4),PP(4),AtA(4,4),EP_tan(4),
     &       sig_jrate(4),dLtan(4,4),Dtot(4),
     &       diagU(4),LtA(4,4)
c
c......... Local variables ..................
c
      real*8 F_matx(3,3),Finv_matx(3,3),detF,
     &       trace_KH,
     &       Sb,Eb             ! (Sb := eff_sts_n; Eb := eff_ste_n)
      real*8 C1111,C1122,C2211,
     &       C2222,C3333,C1212,
     &       C3311,C3322,
     &       Ed,Es,G0
      real*8 Eb0,bigN,gsoft       ! (gsoft:= Sb0 (1+ Eb/Eb0)^N/gtem )
      real*8 Ebtot_c,Ebtot,dG_Eb  ! (dG_Eb = d gsoft / d Eb )
      real*8 xxii,                ! (xxii = theta*dt*hh*pE_ptau) 
     &       hh,                  ! (hh is newly derived)
     &       reg,                 ! (reg := Sb/gsoft)
     &       smm
c
c..............................................
c     
c    gsoft: the rate-dependent relation of eff_strain with
c           the thermal softening effect;
c
c.........Local constant .......................
c
      real*8 eps1,eps2,
     &       Eb0tot,Sb0,ecc,chi,xi,
     &       pCp,ctcon,dG_T,pE_ptau,pE_evt,
     &       pE_Tvt,
     &       rho0,cc1,cc2,rambda,Temp_0
c      
      integer jp,jpt,i,j
      integer mLoop,n_three
c
      logical elastic
c
c........ common blocks ..................
c
      real*8 dt,
     &       el_E,el_V,el_G,el_K,
     &       pl_k0,pl_H,pl_EP,
     &       rc,r0,xc1,xc2,zimp,cr1,cr2
c
      common /connect/Lnp,Lmnp,Lgp,Lmgp,lnods
      common /shapeK/shpk,shpkdx,shpkdy,shpn
      common /elastic/ el_E,el_V,el_G,el_K
      common /plastic/pl_k0,pl_H,pl_EP
      common /hyperelast/rho0,cc1,cc2,rambda
      common /temperature/Temp_0,Cp
      common /shear/r0,xc1,xc2,zimp,cr1,cr2
c
c------------------------------------------------------
c
c......Assign the material constants, and
c      the intrinsic constants
c
c-------------------------------------------------------
c
      n_three = 3
      eps1    = 1.0d-16  ! ( the criteria for Jacobian)
      eps2    = 1.0d-16  ! ( the criteria for Sb)
c
c.....Good redundency................
c
      el_E    = 2.0d+11
      el_V    = 0.30
      el_G    = el_E/(2.0*(1.0 + el_V))
      el_K    = el_E/(3.0*(1.0 - 2.0 * el_V))
      Temp_0  = 293.00    ! K Temp_0
      Sb0     = 2.0d+9    ! pl_k0
      rho0    = 7830      ! kg m^{-3}
c
      elastic = .false.
c
      Eb0     = Sb0/el_E ! pl_k0/el_E
      Eb0tot  = 0.001
      bigN    = 0.01
      smm     = 70.0
c
c.....Temperature parameters (heat conductivity k = 34.6)
c
      alpha_T = 11.20d-6  ! K^{-1}
      delta   =  0.800
      theta   =  0.50
      kappa   = 500.00
      Cp      = 448.00    ! J(kg.K)^{-1}
      chi     =   0.95
c
c
c......Implement the imperfection into the viscoplasticity
c      stress accumulate function
c
c
      rc  = (xc1 - xkc1)**2 + (xc2 - xkc2)**2
      Sb0 = Sb0*(1.0d0 - zimp*exp(-rc/(r0*r0)))
c
c
c........... Calculate the deformation gradient
c            and calculate the rate of deformation 
c            at T = dt * n ( the last time step ? )
c
         do i = 1, 3
	    do j = 1,3
	       sts_PK1(i,j)   = 0.0
               F_matx(i,j)    = 0.0
               Finv_matx(i,j) = 0.0
	    enddo
         enddo
c
         do i = 1,2
	    do j = 1,2
	       dvD(i,j) = 0.0d0
	    enddo
	 enddo
c
	 do i = 1,4
	    do j = 1,4
	       AtA(i,j)   = 0.0
	       LtA(i,j)   = 0.0
	       dLtan(i,j) = 0.0
               dmat(i,j)  = 0.0
	    enddo
	    Dtot(i)       = 0.0
	    sig_jrate(i)  = 0.0
	    EP_tan(i)     = 0.0
	    PP(i)         = 0.0
            smlp(i)       = 0.0
            Dsts_KH(i)    = 0.0
	 enddo
c
      mLoop = Lgp(ik)
c
         do jp = 1, mLoop
            jpt = Lmgp(jp,ik)
c
            F_matx(1,1) = F_matx(1,1)
     &                  + shpkdx(jp,ik)*disp(1,jpt)
            F_matx(1,2) = F_matx(1,2) 
     &                  + shpkdy(jp,ik)*disp(1,jpt)
c
            F_matx(2,1) = F_matx(2,1) 
     &                  + shpkdx(jp,ik)*disp(2,jpt)
            F_matx(2,2) = F_matx(2,2) 
     &                  + shpkdy(jp,ik)*disp(2,jpt)
c
c.......(2)..............
c
            dvD(1,1) = dvD(1,1) + shpkdx(jp,ik)*vel(1,jpt)
            dvD(1,2) = dvD(1,2) + shpkdy(jp,ik)*vel(1,jpt)
c
            dvD(2,1) = dvD(2,1) + shpkdx(jp,ik)*vel(2,jpt)
            dvD(2,2) = dvD(2,2) + shpkdy(jp,ik)*vel(2,jpt)
c
         enddo
c
c.......Add the Kernonical delta part:
c
         do i = 1, 3
            F_matx(i,i) = F_matx(i,i) + 1.0
         enddo
c         
c         
c----------------------------------------------------------
c
c................ End of Calculation .............. 
c
c----------------------------------------------------------
c
c.....Inversion of Deformation Gradient
c
      call invers_matx33(F_matx,Finv_matx,detF)
c
      if (detF .le. eps1) then
          print *, 'Negative Jacobian', detF,' at Gauss pt:', ik
	  stop
      endif
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
c	 End of calculation rate of deformation 
c        for large deformation
c        or small deformation
c----------------------------------------------------------
c	                        
c
c......Calculate the Elastic Muduli tensor
c
c
c      lambda = K - (2/3)*G
c
c.....................................................
c
c    The Voigt Notation:  (See Malvern page 285)
c    
c    1     2    3    4     5     6
c
c    11    22   33   23    31   12
c
c.......................................................
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
c
      Ed = el_E *(1.0 - el_V) /
     &     ((1.0 + el_V)*(1.0d0 - 2.0d0 * el_V))
      Es = el_E * el_V /
     &     ((1.0d0 + el_V)*(1.0d0 - 2 * el_V))
      G0 = el_E/(1.0d0 + el_V)
c
      C1111 = Ed
      C1122 = Es
      C2211 = Es
      C2222 = Ed
      C1212 = G0
      C3333 = Ed  
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
c.....Initialize ``dLtan''................
c
      do i = 1, 4
	 do j = 1,4
	    dLtan(i,j) = dmat(i,j) 
	 enddo
      enddo
c
      diagU(1) = 1.00
      diagU(2) = 1.00
      diagU(3) = 1.00
      diagU(4) = 0.00
c
c
      Eb      = eff_ste
      Sb      = eff_sts
      xi      = 1.00/(rho0*Cp)
c
c.....Initial or elastic value........
c
      if (Sb .le. eps2) then
         xxii    = 0.0d0
         ctcon   = 0.0d0
         ecc     = 0.0d0
         Ebtot   = 0.0d0
c
         elastic = .true.
	 go to 138
      endif
c
c
c----------------------------
c
c  Analysis of Stress
c
c---------------------------
c
       trace_KH = 0.0
c
       do i = 1, 3
          trace_KH = trace_KH + sts_KH(i,i)
       enddo
       trace_KH = trace_KH/3.0
c
c------------------------------------------------
c    
c    Voigt notation
c
c    1     2     3     4     5     6
c
c    11    22    33   23     31    12
c
c------------------------------------------------
c
       Dsts_KH(1) = sts_KH(1,1) - trace_KH
       Dsts_KH(2) = sts_KH(2,2) - trace_KH
       Dsts_KH(3) = sts_KH(3,3) - trace_KH 
       Dsts_KH(4) = sts_KH(1,2) 
c
c
c------------------------------------------------
c
c.......... Calculate C:p and p:C:p .............
c
c------------------------------------------------
c
      do i = 1, 4
         smlp(i)  = 3.0 * Dsts_KH(i)/(2.0 * Sb)
      enddo
c
      do i = 1, 4
	 do j = 1, 4
	    PP(i) = PP(i) + dmat(i,j)*smlp(j)
         enddo
      enddo
c
c..Note.... pCp .........
c
      do i = 1,4
         pCp = pCp + smlp(i)*PP(i)  
      enddo
c
c------------------------------------------------
c
c  pE_ptau:= partial Ebtot / partial tau
c  dG_Eb  := d gsoft / d Ebar
c  dG_T   := d gsoft / d T
c
c  pE_evt := (partial Ebtot/partial Ebar)/(partial Ebtot/partial tau)
c  pE_Tvt := (partial Ebtot/partial T )/(partial Ebtot/partial tau)
c
c--------------------------------------------------
c
c......Implement the imperfection into the viscoplasticity
c      stress accumulate function and calculate
c      the evolution function
c
      gsoft = Sb0 *  (1.0 + Eb/Eb0)**bigN
      gsoft = gsoft* (1.0 - delta*(dexp((Temp-Temp_0)/kappa)
     &      - 1.0d0))
c
      reg   = Sb/gsoft
      Ebtot = Eb0tot * (reg)**(smm)
c
      pE_ptau = smm * Ebtot/Sb
c
      dG_Eb   = bigN * gsoft/(Eb + Eb0)
      dG_T    = - (Sb0 * delta/kappa)*((1.0 + Eb/Eb0)**bigN)
     &        * dexp((Temp - Temp_0)/kappa)
c
      pE_evt  = - dG_Eb * (Sb/gsoft)
      pE_Tvt  = - dG_T  * (Sb/gsoft)
c
c
      hh    = pCp 
     &      - pE_evt
     &      - pE_Tvt*chi*xi*Sb
c     
      xxii  = theta * dt * hh * pE_ptau
      ctcon = xxii/((1.0 + xxii)*hh)
c
c------------------------------------------------
c
c       Calculate the Tangent Moduli
c
c-----------------------------------------------
c
      do i = 1,4
	 do j = 1,4
	    AtA(i,j) = PP(i) * PP(j)
	    LtA(i,j) = diagU(i) * PP(j)
	 enddo
      enddo
c
      do i = 1,	4
         EP_tan(i)  = (Ebtot/ (1.0 + xxii ))
     &              *( PP(i) 
     &              + 3.0*el_K * alpha_T *xi * chi * Sb 
     &              * diagU(i)) 
c
	 do j = 1, 4
	    dLtan(i,j) = dLtan(i,j) - ctcon * (AtA(i,j)
     &                 + 3.0*el_K * alpha_T * chi * xi
     &                 * Sb * LtA(i,j)) 
	 enddo
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
      ecc  = 0.0
      do i = 1, 4
	 ecc = ecc + PP(i) * Dtot(i)
      enddo
c
c
 138  continue   ! ................... !
c
c......Strain rate update...................
c
      Ebtot_c = Ebtot/(1.0 + xxii) + ctcon * ecc
      eff_ste = Eb + Ebtot_c * dt
c
c......Stress update.......
c
      do i = 1,4
	 do j = 1, 4
            sig_jrate(i) = sig_jrate(i) 
     &                   + dLtan(i,j)*Dtot(j)
	 enddo
      enddo
c
c
      if (.not. elastic) then
	 do i = 1, 4
            sig_jrate(i) = sig_jrate(i) - EP_tan(i)
	 enddo
      endif    
c
c-------------------------------------------------------
c
c..........Update the Kirchhoff stress.............
c
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
c
c-----------------------------------------------
c........Update the effective stress
c
c-----------------------------------------------
c
 238   continue
c
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
       Sb  = 0.5 * (
     &           (sts_KH(1,1)-sts_KH(2,2))**2
     &     +     (sts_KH(2,2)-sts_KH(3,3))**2 
     &     +     (sts_KH(3,3)-sts_KH(1,1))**2
     &             ) 
     &     + 3.0 *sts_KH(1,2)**2 
c
       Sb  = sqrt(Sb)       
c    
       eff_sts = Sb
c
c......Temperature update..................
c
      Temp_rate = chi * xi * Sb * Ebtot_c 
      Temp      = Temp + Temp_rate * dt
c
c
c........Update the PK1 = F^{-1} * sts_KH
c
      call mul_matx(Finv_matx,sts_KH,sts_PK1,
     &              n_three,n_three,n_three)
c
      return
      end
c
c
