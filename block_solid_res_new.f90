!c	cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c       S. K. ALIABADI
!c       modified for linear elastic equation, Lucy Zhang 4/22/99
!c       cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine block_solid_res(xloc,dloc, w, p, ien,nsd,nen,ne,nn,nquad,wq,sq,&
				x_pre1,x_pre2,solid_acc,pre,solid_bcvel,solid_bcvel_old,solid_vel)
!	use fluid_variables, only: nsd,nen,ne,nn,nquad,wq,sq
	use global_constants
	use run_variables, only: dt
	use r_common, only: group_young, Poisson, density_solid
	use fluid_variables, only : den_liq, vis_liq,gravity
	use mpi_variables,only: myid
	implicit none

	integer ien(ne,nen)
	real* 8 xloc(nsd,nn), dloc(nsd,nn)
	real* 8 x(nsd,nen), d(nsd,nen), e(nsd,nen)
	real* 8 p(nsd,nn)
	real(8) x_pre1(nsd,nn)
        real(8) x_pre2(nsd,nn)
	real(8) acc(nsd,nen)
	real(8) d1(nsd,nen), d2(nsd,nen)
	real(8) solid_acc(nsd,nn)
	real(8) pre(nn)
	real(8) solid_bcvel(nsd,nn)
        real(8) solid_bcvel_old(nsd,nn)
	real(8) solid_vel(nsd,nn)
        real(8) g(1:nsd)
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
	real(8) dpre(nsd)
	real(8) pre_loc(nen)
	real(8) acc_fake(nsd,nen)
	real(8) a_fake(nsd)
	real(8) us(nsd)
	real(8) uf(nsd)
	real(8) vel_s(nsd,nen)
	real(8) vel_f(nsd,nen)
	real(8) dusdx(nsd,nsd)
	real(8) dufdx(nsd,nsd)
	real(8) taus(nsd,nsd)
        real(8) tauf(nsd,nsd)


!-----------------------------------	
	real(8) rho_solid	! solid density should be passed from input file 
				! same as blockgmres_solid.f90
!-----------------------------------

!	p(:,:) = 0.0d0
	w(:,:) = 0.0d0

! change E, \nv to \lammda and \mu (Lame parameters)
!	mu = group_young(1)/((1+Poisson)*(1-2*Poisson))
!	la = group_young(1)/(2*(1+Poisson))
! set solid densitiy to be zero
! Here is the differenc between block_solid and blockgmres_solid, since now I do not consider 
! the density related terms in the FSI force which is the semi-implicit FSI algorithm
!	rho_solid =  1.0d0
	rho_solid=density_solid+den_liq
if(myid==0)write(*,*)'rho_solid=',rho_solid,'den_liq=',den_liq
g(1:nsd)=gravity(1:nsd)
!------------------------
	mu = 0.0d0
	la = 0.0d0

!------------------------
        do ie=1,ne 
	   do inl=1,nen
	      do isd=1,nsd
		 x(isd,inl) = xloc(isd,ien(ie,inl))
		 d(isd,inl) = dloc(isd,ien(ie,inl))
		 d1(isd,inl) = x_pre1(isd,ien(ie,inl)) 
		 d2(isd,inl) = x_pre2(isd,ien(ie,inl)) 
!		 acc(isd,inl) = solid_acc(isd,ien(ie,inl)) - solid_accel_old(isd,ien(ie,inl))
                 acc(isd,inl) = solid_acc(isd,ien(ie,inl))

                 acc_fake(isd,inl) = (solid_bcvel(isd,ien(ie,inl)) - solid_bcvel_old(isd,ien(ie,inl))) / dt
		vel_s(isd,inl) = solid_vel(isd,ien(ie,inl))
		vel_f(isd,inl) =solid_bcvel(isd,ien(ie,inl))
	      enddo
		pre_loc(inl) = pre(ien(ie,inl))
	   enddo

!	   acc(1:nsd,1:nen) =  (d(1:nsd,1:nen) - 2.0*d1(1:nsd,1:nen) + d2(1:nsd,1:nen)) / (dt**2)

!	if (acc(3,1) .ne. 0) then 
!write(*,*) 'acc3', acc(3,1), 'z coor', x(3,1)
!continue
!end if
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c        CALCULATE kd-f
!c                  in this problem f=0
!c        ONLY CALCULATE  kd in local coordinates
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
	
		dpre(:) =0.0
		a_fake(:) =0.0
		us(:) = 0.0
		uf(:) = 0.0
		dusdx(:,:) = 0.0
		dufdx(:,:) = 0.0
		

	      do inl=1,nen
		 ax = ax + sh(0,inl)*acc(1,inl)*rho_solid
		 ay = ay + sh(0,inl)*acc(2,inl)*rho_solid
	if (nsd == 3) then
		 az = az + sh(0,inl)*acc(3,inl)*rho_solid
	end if


		do isd = 1,nsd
			dpre(isd) = dpre(isd) + sh(isd,inl) * pre_loc(inl)
			a_fake(isd) = a_fake(isd) + sh(0,inl) * acc_fake(isd,inl) * rho_solid
			us(isd) = us(isd) + sh(0,inl) * vel_s(isd,inl)
			uf(isd) = uf(isd) + sh(0,inl) * vel_f(isd,inl)
			do jsd=1,nsd
				dusdx(isd,jsd) = dusdx(isd,jsd) + sh(jsd,inl) * vel_s(isd,inl) ! du_i / dx_j
				dufdx(isd,jsd) = dufdx(isd,jsd) + sh(jsd,inl) * vel_f(isd,inl) ! du_i / dx_j
			end do
		end do

		do isd = 1,nsd
			do jsd = 1,nsd
				taus(isd,jsd) = vis_liq * (dusdx(isd,jsd) + dusdx(isd,jsd))
                                tauf(isd,jsd) = vis_liq * (dufdx(isd,jsd) + dufdx(isd,jsd))
			end do
		end do



!c............... calculate the first derivative
		 do isd=1,nsd
		    drx(isd)=drx(isd)+sh(1,inl)*d(isd,inl)      
		    dry(isd)=dry(isd)+sh(2,inl)*d(isd,inl)     
        if (nsd == 3) then ! To be justified for both  2, 3 - D cases 
		    drz(isd)=drz(isd)+sh(3,inl)*d(isd,inl)
	end if
		 enddo
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
		 p(xsd,node) = p(xsd,node)  &
!		      ph(xsd,inl) * ttt + &
!		      ph(xsd,inl) * txx + &
!		      ph(ysd,inl) * tyx + &
!		      ph(zsd,inl) * tzx   &
		     -ph(0,inl) * ax	  &
!		     +ph(0,inl) * dpre(xsd) &
!---------------------------------------------------------------------------------
                        +ph(0,inl) * a_fake(xsd)  &
                        +ph(0,inl) * (us(xsd)*dusdx(xsd,xsd) + us(ysd)*dusdx(xsd,ysd) + us(zsd)*dusdx(xsd,zsd)) &
                        -ph(0,inl) * (uf(xsd)*dufdx(xsd,xsd) + uf(ysd)*dufdx(xsd,ysd) + uf(zsd)*dufdx(xsd,zsd)) &
                        +ph(xsd,inl)*taus(xsd,xsd) + ph(ysd,inl)*taus(xsd,ysd) + ph(zsd,inl)*taus(xsd,zsd) &
                        -(ph(xsd,inl)*tauf(xsd,xsd) + ph(ysd,inl)*tauf(xsd,ysd) + ph(zsd,inl)*tauf(xsd,zsd))
!-----------------------------------------------------------------------------------



		 p(ysd,node) = p(ysd,node)  &
!		      ph(ysd,inl) * ttt + &
!		      ph(xsd,inl) * txy + &
!		      ph(ysd,inl) * tyy + &
!		      ph(zsd,inl) * tzy   &
		     -ph(0,inl) * ay	  &
!                     +ph(0,inl) * dpre(ysd) &
!---------------------------------------------------------------------------------
                        +ph(0,inl) * a_fake(ysd)  &
                        +ph(0,inl) * (us(xsd)*dusdx(ysd,xsd) + us(ysd)*dusdx(ysd,ysd) + us(zsd)*dusdx(ysd,zsd)) &
                        -ph(0,inl) * (uf(xsd)*dufdx(ysd,xsd) + uf(ysd)*dufdx(ysd,ysd) + uf(zsd)*dufdx(ysd,zsd)) &
                        +ph(xsd,inl)*taus(ysd,xsd) + ph(ysd,inl)*taus(ysd,ysd) + ph(zsd,inl)*taus(ysd,zsd) &
                        -(ph(xsd,inl)*tauf(ysd,xsd) + ph(ysd,inl)*tauf(ysd,ysd) + ph(zsd,inl)*tauf(ysd,zsd))
!-----------------------------------------------------------------------------------



		 p(zsd,node) = p(zsd,node)  &
!		      ph(zsd,inl) * ttt + &
!		      ph(xsd,inl) * txz + &
!		      ph(ysd,inl) * tyz + &
!		      ph(zsd,inl) * tzz   &
		     -ph(0,inl) * az	  &
!                     +ph(0,inl) * dpre(zsd) &
!---------------------------------------------------------------------------------
                        +ph(0,inl) * a_fake(zsd)  &
                        +ph(0,inl) * (us(xsd)*dusdx(zsd,xsd) + us(ysd)*dusdx(zsd,ysd) + us(zsd)*dusdx(zsd,zsd)) &
                        -ph(0,inl) * (uf(xsd)*dufdx(zsd,xsd) + uf(ysd)*dufdx(zsd,ysd) + uf(zsd)*dufdx(zsd,zsd)) &
                        +ph(xsd,inl)*taus(zsd,xsd) + ph(ysd,inl)*taus(zsd,ysd) + ph(zsd,inl)*taus(zsd,zsd) &
                        -(ph(xsd,inl)*tauf(zsd,xsd) + ph(ysd,inl)*tauf(zsd,ysd) + ph(zsd,inl)*tauf(zsd,zsd))
!-----------------------------------------------------------------------------------





	      enddo
	end if 

	if (nsd == 2) then ! 2-D case, plain strain model
	   do inl=1,nen
		node=ien(ie,inl)
		p(xsd,node)=p(xsd,node)  &
!			ph(xsd,inl)*(la+2*mu)*drx(xsd) + &
!			ph(ysd,inl)*mu*dry(xsd) + &
!			ph(xsd,inl)*la*dry(ysd) + &
!			ph(ysd,inl)*mu*drx(ysd)   &
			+ph(0,inl) * ax		  &
!                        +ph(0,inl) * dpre(xsd)    &
!---------------------------------------------------------------------------------
			-ph(0,inl) * a_fake(xsd)  &
!			+ph(0,inl) * (us(xsd)*dusdx(xsd,xsd) + us(ysd)*dusdx(xsd,ysd)) &
			-ph(0,inl) * (uf(xsd)*dufdx(xsd,xsd) + uf(ysd)*dufdx(xsd,ysd))*rho_solid &
!			-ph(0,inl) * g(1)*rho_solid &
			-ph(xsd,inl)*taus(xsd,xsd) + ph(ysd,inl)*taus(xsd,ysd) &
			+(ph(xsd,inl)*tauf(xsd,xsd) + ph(ysd,inl)*tauf(xsd,ysd))
!-----------------------------------------------------------------------------------


		p(ysd,node)=p(ysd,node)  &
!			ph(ysd,inl)*la*drx(xsd) + &
!			ph(xsd,inl)*mu*dry(xsd) + &
!			ph(ysd,inl)*(la+2*mu)*dry(ysd) + &
!			ph(xsd,inl)*mu*drx(ysd)   &
                        +ph(0,inl) * ay           &
!                        +ph(0,inl) * dpre(ysd)    &
!----------------------------------------------------------------------------------
                        -ph(0,inl) * a_fake(ysd)  &
!                        +ph(0,inl) * (us(xsd)*dusdx(ysd,xsd) + us(ysd)*dusdx(ysd,ysd)) &
                        -ph(0,inl) * (uf(xsd)*dufdx(ysd,xsd) + uf(ysd)*dufdx(ysd,ysd))*rho_solid &
!			-ph(0,inl) * g(2) * rho_solid &
                        -ph(xsd,inl)*taus(ysd,xsd) + ph(ysd,inl)*taus(ysd,ysd) &
                        +(ph(xsd,inl)*tauf(ysd,xsd) + ph(ysd,inl)*tauf(ysd,ysd))
!----------------------------------------------------------------------------------


	  end do
	end if

	   enddo ! Gausian point loop
	enddo ! element loop
      return
      end





