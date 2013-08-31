!c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Based on Huge's book take solid acceleration as the primary variable
!       solve it and then get the displacement for dynamic problems
!	Only linear elastic, I use \alpha - method
!       Visco-linear elastic, I use Newmark method
!c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine block_solid(xloc,dloc, w, p, ien,nsd,nen,ne,nn,nquad,wq,sq,&
				x_pre1,solid_prevel,solid_preacc,ien_sbc,ne_sbc,solid_stress,mtype)
!	use fluid_variables, only: nsd,nen,ne,nn,nquad,wq,sq
	use global_constants
	use run_variables, only: dt
	use r_common, only: group_young, Poisson, density_solid
	use fluid_variables, only: den_liq,gravity
	use solid_variables, only: damp_solid
	implicit none
!-----------------------------------------------
! xloc --- the initial configuration
! dloc --- solid acceleration at n+1 to be solved
! solid_preacc --- solid acceleration at n 
! solid_vel --- solid velocity at n+1
! solid_prevel --- solid velocity at n
! alpha, beta, gamma --- constant parameters definition seen Huge 532

	integer ien(ne,nen)
	real* 8 xloc(nsd,nn), dloc(nsd,nn)
	real* 8 x(nsd,nen), d(nsd,nen), e(nsd,nen)
	real* 8 p(nsd,nn)
	real(8) x_pre1(nsd,nn)
	real(8) solid_prevel(nsd,nn)
	real(8) solid_preacc(nsd,nn)
	integer mtype(ne)
!        real(8) x_pre2(nsd,nn)
	real(8) acc(nsd,nen), acc_pre(nsd,nen)
	real(8) d1(nsd,nen), d2(nsd,nen)
	real(8) vel(nsd,nen), vel_pre(nsd,nen)

!----------------------------------------
	integer nsd
	integer nen
	integer ne
	integer nn
	integer nquad
	real(8) wq(8)
	real(8) sq(0:3,8,8)
!----------------------------------------
	real* 8 eft0,det,eft1
	real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
	real* 8 xr(nsd,nsd), cf(nsd,nsd),sx(nsd,nsd)

	real* 8 drx(nsd),dry(nsd),drz(nsd)
	real* 8 ttt,txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz

	real* 8 mu,la
	real* 8 erx,ery,erz,ers

!c....   stiffness matrix
	real* 8 w(nsd,nn)  
	integer inl, ie, isd, iq, node, jsd
	real(8) ax, ay, az
!-----------------------------------	
	real(8) rho_solid	! solid density should be passed from input file 
				! same as blockgmres_solid.f90
	real(8) alpha
	real(8) beta
	real(8) gama 
!-----------------------------------
	real(8) solid_bcforce(nsd,nn) ! external force acting at the solid boundary
	real(8) solid_stress(nsd*2,nn) ! fluid stress acting on the solid boundary including pressure and viscous stress
	integer ne_sbc			! solid elements on the boundary
	integer ien_sbc(ne_sbc,nen+2)
!------------------------------------
! For viscoelastic model
	real(8) dsdt(nsd,nsd) ! time derivative of strain and then times damping

	real(8) g(nsd)

g(1:nsd)=gravity(1:nsd)

!------------------------------------
p(:,:) = 0.0d0
w(:,:) = 0.0d0

!---------------------------------
! define the numerical parameters
alpha = -0.05 ! -1/3 < alpha < 0 and alpha == 0 is Newmark method
gama = (1.0 - 2.0 * alpha) * 0.5
beta = ( (1.0 - alpha)**2 ) * 0.25
!beta = (1.0 - alpha **2) * 0.25
!----------------------------------


! set solid densitiy
	rho_solid = density_solid + den_liq

!	rho_solid = 0.0
        do ie=1,ne 
! change E, \nv to \lammda and \mu (Lame parameters)
        mu = group_young(mtype(ie))/((1+Poisson)*(1-2*Poisson))
        la = group_young(mtype(ie))/(2*(1+Poisson))


	   do inl=1,nen
	      do isd=1,nsd
		 x(isd,inl) = xloc(isd,ien(ie,inl)) ! local initial configuration
		 acc(isd,inl) = dloc(isd,ien(ie,inl)) ! local acc to be solved at n+1
		 acc_pre(isd,inl) = solid_preacc(isd,ien(ie,inl)) ! local acc at n passed in
		 vel_pre(isd,inl) = solid_prevel(isd,ien(ie,inl)) ! local vel at n
		 d1(isd,inl) = x_pre1(isd,ien(ie,inl)) - x(isd,inl) ! displacement at n

!                vel(isd,inl) = solid_vel(isd,ien(ie,inl))        ! local vel at n+1
!                d2(isd,inl) = x_pre2(isd,ien(ie,inl)) - x(isd,inl) ! displacement at n-1

	      enddo
	   enddo
!-------------------------------------------------
! solid displacement at n+1 to be solved
	d(:,:) = d1(:,:) + dt*vel_pre(:,:) + (dt**2)*0.5*( (1.0-2.0*beta)*acc_pre(:,:) + 2.0*beta*acc(:,:) )
! solid velocity at n+1 to be solved
	vel(:,:) = vel_pre(:,:) + dt*( (1-gama)*acc_pre(:,:) + gama*acc(:,:) )
! do weighted average using alpha
	d(:,:) = (1 + alpha)*d(:,:) - alpha*d1(:,:)
	vel(:,:) = (1 + alpha)*vel(:,:) - alpha*vel_pre(:,:)
! adding damping term to displacement for viscoelastic model just comment it out
!	d(:,:) = d(:,:) + damp_solid * vel(:,:) 

!-------------------------------------------------
	   do iq=1,nquad
	if (nsd == 3) then ! 3-D case
	      if (nen.eq.4) then
		 include "sh3d4n.h"
	      else if (nen.eq.8) then
		 include "sh3d8n.h"
	      end if
	end if
	if (nsd == 2) then ! 2-D case
	      if (nen.eq.3) then !calculate shape function at quad point
                 include "sh2d3n.h"
              elseif (nen.eq.4) then
                 include "sh2d4n.h"
              endif
	end if
	      eft0 = abs(det) * wq(iq)

!=============================================
!c.... no jacobian
!c	      eft0 = wq(iq)   
!c............ Calculate local Stiffness Matrix k (kd=f)
! Use the diagonal to be the pre-conditioner, see nodes
	      do inl=1,nen
	if (nsd == 3) then ! 3-D case
		 txx = sh(1,inl)**2
		 tyy = sh(2,inl)**2
		 tzz = sh(3,inl)**2
		 ttt = txx + tyy + tzz 
		 node = ien(ie,inl)
		 w(xsd,node)=w(xsd,node)+mu*(ttt+txx)*eft0
		 w(ysd,node)=w(ysd,node)+mu*(ttt+tyy)*eft0
		 w(zsd,node)=w(zsd,node)+mu*(ttt+tzz)*eft0 
		 
		 w(xsd,node)=w(xsd,node)+la*txx*eft0
		 w(ysd,node)=w(ysd,node)+la*tyy*eft0
		 w(zsd,node)=w(zsd,node)+la*tzz*eft0 

		 w(xsd,node)=w(xsd,node)*(dt**2)*beta+rho_solid*eft0*sh(0,inl)
                 w(ysd,node)=w(ysd,node)*(dt**2)*beta+rho_solid*eft0*sh(0,inl)
                 w(zsd,node)=w(zsd,node)*(dt**2)*beta+rho_solid*eft0*sh(0,inl)

	end if

	if (nsd == 2) then ! 2-D case
		node = ien(ie,inl)
		w(xsd,node)=w(xsd,node)+(sh(1,inl)**2)*(la+2*mu)+(sh(2,inl)**2)*mu
		w(ysd,node)=w(ysd,node)+(sh(1,inl)**2)*mu+(sh(2,inl)**2)*(la+2*mu) 

                 w(xsd,node)=w(xsd,node)*(dt**2)*beta+rho_solid*eft0*sh(0,inl)
                 w(ysd,node)=w(ysd,node)*(dt**2)*beta+rho_solid*eft0*sh(0,inl)
	end if                
	      enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        CALCULATE kd-f
!c                  in this problem f=0
!c        ONLY CALCULATE  kd in local coordinates
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------------------------
!
!Add gravity for f
!Chu Wang 2013.3
!----------------------------------
!c.........  initialize variables
	      erx = 0.0
	      ery = 0.0
	      erz = 0.0

	      do isd = 1,nsd
		 drx(isd) = 0.0
		 dry(isd) = 0.0
	if (nsd == 3) then ! To be justified for both  2, 3 - D cases
		 drz(isd) = 0.0
	end if 
	      enddo


! Accerlearation at G point
		ax = 0.0
		ay = 0.0
		az = 0.0


	      do inl=1,nen
		 ax = ax + sh(0,inl)*acc(1,inl)*rho_solid
		 ay = ay + sh(0,inl)*acc(2,inl)*rho_solid
	if (nsd == 3) then
		 az = az + sh(0,inl)*acc(3,inl)*rho_solid
	end if
!c............... calculate the first derivative
		 do isd=1,nsd
		    drx(isd)=drx(isd)+sh(1,inl)*d(isd,inl)      
		    dry(isd)=dry(isd)+sh(2,inl)*d(isd,inl)     
        if (nsd == 3) then ! To be justified for both  2, 3 - D cases 
		    drz(isd)=drz(isd)+sh(3,inl)*d(isd,inl)
	end if
		 enddo
	      end do
 

!--------------------------------------
! calculate damping term for viscoelastic material
		dsdt(:,:) = 0.0
		do inl=1,nen
			do isd=1,nsd
				do jsd=1,nsd
				dsdt(isd,jsd) = dsdt(isd,jsd) + &
				( sh(isd,inl)*vel(jsd,inl) + sh(jsd,inl)*vel(isd,inl) ) * 0.5 * damp_solid
				end do
			end do
		end do

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	      do inl=1,nen
		 ph(:,inl) = sh(:,inl)*eft0
!		 ph(1,inl) = sh(1,inl)*eft0
!		 ph(2,inl) = sh(2,inl)*eft0
!		 ph(3,inl) = sh(3,inl)*eft0
	      enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.....Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	if (nsd == 3) then ! 3-D case
	      ttt = la*(drx(xsd)+dry(ysd)+drz(zsd))
	      txx = mu*(drx(xsd)+drx(xsd))
	      tyx = mu*(dry(xsd)+drx(ysd))
	      tzx = mu*(drz(xsd)+drx(zsd))
	      txy = mu*(drx(ysd)+dry(xsd))
	      tyy = mu*(dry(ysd)+dry(ysd))
	      tzy = mu*(drz(ysd)+dry(zsd))
	      txz = mu*(drx(zsd)+drz(xsd))
	      tyz = mu*(dry(zsd)+drz(ysd))
	      tzz = mu*(drz(zsd)+drz(zsd))
	      do inl=1,nen
		 node=ien(ie,inl)
!c.....Elastic Equation (calculate residual: r=kd=p)
		 p(xsd,node) = p(xsd,node) + &
		      ph(xsd,inl) * ttt + &
		      ph(xsd,inl) * txx + &
		      ph(ysd,inl) * tyx + &
		      ph(zsd,inl) * tzx   &
		     +ph(0,inl) * ax-ph(0,inl)*rho_solid*g(xsd)
		 p(ysd,node) = p(ysd,node) + &
		      ph(ysd,inl) * ttt + &
		      ph(xsd,inl) * txy + &
		      ph(ysd,inl) * tyy + &
		      ph(zsd,inl) * tzy   &
		     +ph(0,inl) * ay-ph(0,inl)*rho_solid*g(ysd)
		 p(zsd,node) = p(zsd,node) + &
		      ph(zsd,inl) * ttt + &
		      ph(xsd,inl) * txz + &
		      ph(ysd,inl) * tyz + &
		      ph(zsd,inl) * tzz   &
		     +ph(0,inl) * az-ph(0,inl)*rho_solid*g(zsd)

	      enddo
	end if 

	if (nsd == 2) then ! 2-D case, plain strain model
	   do inl=1,nen
		node=ien(ie,inl)
		p(xsd,node)=p(xsd,node) + &
			ph(xsd,inl)*(la+2*mu)*drx(xsd) + &
			ph(ysd,inl)*mu*dry(xsd) + &
			ph(xsd,inl)*la*dry(ysd) + &
			ph(ysd,inl)*mu*drx(ysd)   &
			+ph(0,inl) * ax-ph(0,inl)*rho_solid*g(xsd)

		p(ysd,node)=p(ysd,node) + &
			ph(ysd,inl)*la*drx(xsd) + &
			ph(xsd,inl)*mu*dry(xsd) + &
			ph(ysd,inl)*(la+2*mu)*dry(ysd) + &
			ph(xsd,inl)*mu*drx(ysd)   &
                        +ph(0,inl) * ay-ph(0,inl)*rho_solid*g(ysd)

	  end do
	end if

!---------------------------------
! Adding damping term for viscoelastic model
	do inl=1,nen
	node = ien(ie,inl)
	do isd = 1, nsd
		do jsd = 1,nsd
		p(isd,node) = p(isd,node) + ph(jsd,inl) * dsdt(isd,jsd)
		end do
	end do
	end do
!----------------------------------


	   enddo ! Gausian point loop
	enddo ! element loop

!======================================
! Apply 2nd type boundary
	if (nsd ==  2) then 
	call apply_2ndbc_solid2d(x_pre1,nsd,nn,ien_sbc,ne_sbc,nen,ien,ne,solid_bcforce,solid_stress)
	p(:,:) = p(:,:) + solid_bcforce(:,:)
	else
	call apply_2ndbc_solid(x_pre1,nsd,nn,ien_sbc,ne_sbc,nen,ien,ne,solid_bcforce,solid_stress)
!~---------------------------------
! Manually apply B.C.
!	solid_bcforce(:,:) =0.0 
!	do inl = 1, 1376, 11
!	solid_bcforce(2,inl) = 1.0  
!	end do
!---------------------------------
	
	p(:,:) = p(:,:) + solid_bcforce(:,:)
	end if


      return
      end





