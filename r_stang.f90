subroutine r_stang(solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_accel,solid_pave,solid_stress,solid_strain)
  use run_variables, only: ntsbout,its
  use solid_variables, only: nsd_solid,ne_solid,nn_solid,nen_solid,nsurface,nquad_solid,xq_solid,wq_solid,nquadpad_solid
  use r_common
  implicit none

  integer,dimension(1:ne_solid,1:nen_solid) :: solid_fem_con   !...connectivity for solid FEM mesh
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_vel         !...velocity
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_accel       !...acceleration
  real(8),dimension(nn_solid)   :: solid_pave  !...averaged solid pressure (from mixed formulation -> ???)
  real(8),dimension(1:nsd_solid*2,nn_solid) :: solid_stress  !...solid stress (Voigt notation)
  real(8),dimension(1:nsd_solid*2,nn_solid) :: solid_strain  !...solid strain (Voigt notation)

  !...ox means X
  !... x means x
  real(8) :: x(nsd_solid,nen_solid)    !...element nodes initial position
  real(8) :: y(nsd_solid,nen_solid)    !...element nodes current position
  real(8) :: vel(nsd_solid,nen_solid)  !...element nodes velocity
  real(8) :: acc(nsd_solid,nen_solid)  !...element nodes acceleration
  real(8) :: toc(nsd_solid,nsd_solid)  !...Cauchy Green Deformation Tensor at t=0
  real(8) :: xto(nsd_solid,nsd_solid)  !...F    = dx/dX
  real(8) :: xot(nsd_solid,nsd_solid)  !...F^-1 = dX/dx
  real(8) :: xj(nsd_solid,nsd_solid)   !...
  real(8) :: xji(nsd_solid,nsd_solid)  !...
  real(8) :: rs(nsd_solid)             !...iso-parametric shapefunction value at quadrature point
  real(8) :: toxj(nsd_solid,nsd_solid)
  real(8) :: toxji(nsd_solid,nsd_solid)
  real(8) :: xmj(3),xmi(3),dxmj(3,6),ddxmj(3,6,6),obc(6,6),ocuu(6,6),ocup(6)
  real(8) :: xfrtem(6,6),tem(6),ten(6),ttm(6)
  real(8) :: xkup(3*nen_solid,nup,ne_solid)
  real(8) :: xkpp(nup,nup,ne_solid)
  real(8) :: xfp(nup,ne_solid)
  real(8) :: ge(1:nsd_solid*2,ne_solid,nquadpad_solid)    !...Green strain
  real(8) :: cstr(1:nsd_solid*2,ne_solid,nquadpad_solid)  !...Cauchy stress
  real(8) :: cstr_element(1:nsd_solid*2)   !...Cauchy stress in element
  real(8) :: pre(nup,ne_solid) !...pressure in solid (only used for almost compressible material)
  real(8) :: tot_vol_init,tot_vol_curr
  real(8) :: det    !...Jacobian Determinante of Deformation Gradient
  real(8) :: todet  !...Jacobian Determinante of Deformation Gradient at t=0
  real(8) :: w_init !...Gauss weight pluss Jacobien Det
  real(8) :: w_curr !...Gauss weight pluss Jacobien Det
  real(8) :: wto    !...Potential W
  real(8) :: ocpp
  integer :: nos,ntem         !...counter
  integer :: isd,iq  !...counter
  integer :: ine,in,nu1,mu1,ip,jp  !...counter

  write(*,*) " calculate internal + inertial forces (r_stang)"

  predrf(1:nsd_solid*nn_solid) = 0.0d0
  tot_vol_init = 0.0d0
  tot_vol_curr = 0.0d0

  element: do ine=1,ne_solid

     xfp(1:nump,ine)=0.0d0
     xkup(1:nsd_solid*nen_solid,1:nump,ine)=0.0d0
     xkpp(1:nump,1:nump,ine)=0.0d0
        cstr_element(:)=0.0
     do nos=1,nen_solid
        ntem=solid_fem_con(ine,nos) !...connectivity
        x(1:nsd_solid,nos)   = solid_coor_init(1:nsd_solid,ntem)
        y(1:nsd_solid,nos)   = solid_coor_curr(1:nsd_solid,ntem)
        vel(1:nsd_solid,nos) = solid_vel(1:nsd_solid,ntem)
        acc(1:nsd_solid,nos) = solid_accel(1:nsd_solid,ntem)
     enddo

    !...gauss integration
    !...update 06.03.2003, Axel G.: can handle tetrahedral elements as well
    !...       07.07.2003, Axel G.: integration and shape function the same as fluid
     gauss_int: do iq = 1,nquad_solid

        rs(1:nsd_solid) = xq_solid(1:nsd_solid,iq)

!     isoparametric interpolation
        call r_element(rs)
!     y-(r,s)
        call r_jacob(y,xj,xji,det)
		
!     x-(r,s)
        call r_jacob(x,toxj,toxji,todet)

!     derivative about ox and x
        call r_bdpd_curr(xji)
        call r_bdpd_init(toxji)
!     deformation gradient
        call r_stoxc(xto,xot,xj,xji,toxj,toxji,toc,ine)

!================================================
! Hyperelastic Material --> Option material_type=1
	if (material_type==1) then
!     material j
		call r_smaterj(wto,toc,xmi,xmj,dxmj,ddxmj)
!     discretized pressure
		call r_spress(rs,ine)
!     continuous pressure
        call r_sbpress(dxmj,ddxmj,xmj)
!     material c
        call r_sboc(obc,ocpp,ocuu,ocup,xmj,dxmj,ddxmj)
!     strain
        call r_sstrain(toc,xto,iq,ine,ge)
!     First Piola-Kirchoff stress - P
        call r_spiola(xmj,dxmj,xto) !hyperelastic material
!     correction for viscous fluid stress
        call r_spiola_viscous(xot,vel)  
!     calculate cauchy stress for output
!     ! for force calculation in current configuration (updated Lagrangian) -> remove "if" condition!!!
      if (mod(its,ntsbout) == 0) then
        call r_scauchy(det,todet,xto,cstr_element)
        cstr(1:nsd_solid*2,ine,iq) = cstr_element(1:nsd_solid*2)
      endif
!==========================================================
! Linear elastic Cauchy stress - sigma
	elseif (material_type==2) then
!     strain
        call r_sstrain(toc,xto,iq,ine,ge)
!	  Calculate cauchy stress then transform to 1st PK stress
		call r_spiola_elastic(det,xot,ge,iq,ine,cstr_element)
!     correction for viscous fluid stress
        call r_spiola_viscous(xot,vel)  
!     assemble cauchy stress for output
        if (mod(its,ntsbout) == 0) then
			cstr(1:nsd_solid*2,ine,iq) = cstr_element(1:nsd_solid*2)
		endif
	endif
!===========================================================

!     gauss weight and volume
        w_init = wq_solid(iq) * todet
        w_curr = wq_solid(iq) * det
!     volume
        tot_vol_init = tot_vol_init + w_init
        tot_vol_curr = tot_vol_curr + w_curr
!     internal force and stiffness matrix
	     call r_sstif(ocpp,ocup,xkup,xkpp,xfp,ine,w_init,vel,acc,solid_fem_con)

     enddo gauss_int  
  enddo element

  write(*,'("  total solid volume (init) = ",f12.6)') tot_vol_init
  write(*,'("  total solid volume (curr) = ",f12.6)') tot_vol_curr

 !     pressure condensation, inverse kpp
  if (1 == 0) then 
  do ine=1,ne_solid
     tem(1:nump) = xfp(1:nump,ine)
     xfrtem(1:nump,1:nump) = xkpp(1:nump,1:nump,ine)
     call gaussj(xfrtem,nump,nsd_solid*2,tem,1,1)
     do ip=1,nump
        ttm(ip)=0.0d0
        do nos=1,nen_solid
           do isd=1,nsd_solid
              nu1 = (isd-1)*nn_solid + solid_fem_con(ine,nos)
              mu1 = (isd-1)*nen_solid + nos
              ttm(ip) = ttm(ip) + xkup(mu1,ip,ine)*du(isd,solid_fem_con(ine,nos))
              predrf(nu1) = predrf(nu1) + xkup(mu1,ip,ine)*tem(ip)
              write(*,*) xkup(mu1,ip,ine)*tem(ip)
           enddo
        enddo
     enddo

!...  storage
     do ip=1,nump
        ten(ip) = 0.0d0
        do jp=1,nump
           ten(ip) = ten(ip)+xfrtem(ip,jp)*ttm(jp)
        enddo
        pre(ip,ine) = -tem(ip)-ten(ip)
     enddo
  enddo
  endif

 !...calculate pressure and stress for the structure output
  if (mod(its,ntsbout) == 0) then
   write(*,*) "  calculate stress and strain for output"
   do in=1,nn_solid
     solid_stress(1:nsd_solid*2,in)=0.0d0
     solid_strain(1:nsd_solid*2,in)=0.0d0
     ntem=0
     solid_pave(in)=0
     if (nsd_solid .ne. 0) then  !do not calculate if it is a point
      	do ine = 1,ne_solid
           do nos=1,nen_solid
              if (solid_fem_con(ine,nos) == in) then
                 ntem = ntem + 1
                 solid_pave(in) = solid_pave(in) + pre(1,ine)
			 	 do iq = 1,nquad_solid
				   solid_stress(1:nsd_solid*2,in) = solid_stress(1:nsd_solid*2,in) + wq_solid(iq)*cstr(1:nsd_solid*2,ine,iq) !...constant stress and strain in element
                   solid_strain(1:nsd_solid*2,in) = solid_strain(1:nsd_solid*2,in) + wq_solid(iq)*ge(1:nsd_solid*2,ine,iq)
			  	 enddo
                 goto 541
              endif
           enddo
 541    enddo

        solid_stress(1:nsd_solid*2,in) = solid_stress(1:nsd_solid*2,in)/ntem
        solid_strain(1:nsd_solid*2,in) = solid_strain(1:nsd_solid*2,in)/ntem
        solid_pave(in) = solid_pave(in)/ntem
		endif
   enddo
  endif
  write(*,*) " done                                 (r_stang)"
  return
end subroutine r_stang
