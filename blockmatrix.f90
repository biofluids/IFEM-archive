!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  L. Zhang
!  Northwestern University
!  This subroutine solves for the residual for all degrees of freedom
!____________________________________________________________________
!  L. Zhang, 06/24/2004
!  Tulane University
!  Revised the subroutine to array
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine block(xloc, dloc, doloc, p, q, hk, ien, f_fluids,rngface, f_stress)
  use global_constants
  use run_variables
  use fluid_variables
  implicit none

  integer ien(nen,ne)
  real* 8 xloc(nsd,nn)
  real* 8 dloc(ndf,nn),doloc(ndf,nn)
  real* 8 p(ndf,nn),q(ndf,nn),hk(ne)

  real* 8 x(nsd,nen)
  real* 8 d(ndf,nen),d_old(ndf,nen)

  real* 8 eft0,det,effd,effm,effc
  real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
  real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real* 8 f_stress(nsd,nsd,nn)

  real* 8 drt(ndf),drs(ndf)
  real* 8 dr(nsd,ndf)
  real* 8 u,v,w,pp,ug
  real* 8 tau(nsd,nsd)
  real* 8 hg,taum,tauc,vel,ree, taul
  real* 8 res_c,res_a(nsd),res_t(nsd)
  real* 8 prs_c,prs_t(nsd),p_vec(3),prs_cc(nsd)
  real* 8 mu,nu,ro,g(nsd)
  real* 8 tempc(ndf),temp
  real* 8 dtinv,oma,ama
  integer inl, ie, isd, iq, node,jsd
  integer ieface,irng, rngface(neface,ne) !,inface

  real* 8 f_fluids(nsd,nn)
  real* 8 fnode(nsd,nen),fq(nsd)

  real* 8 K(12,12),F(12)
  integer ii,ij,icount,jcount
  real* 8 uold(nsd)
  real* 8 ddd(12)
  real* 8 residual(12)
  real* 8 res_orig(12)
  integer inode,jnode
  real* 8 KK(ndf*nn,ndf*nn)

  KK(:,:) = 0
  dtinv = 1.0/dt
  if(steady) dtinv = 0.0
  oma   = 1.0 - alpha
  ama   = 1.0 - oma
 !=================================================
!f_fluids(:,:)=f_fluids(:,:)/(0.0625/6.0)
!p(1:nsd,1:nn)=p(1:nsd,1:nn)+f_fluids(1:nsd,1:nn)
!=================================================
!===================================================
! f_fluids is actually the FSI force at fluid nodes,
! then we will just subscrib it from p(!:nsd) which is the 
! residuals for the momentum equations. 
! 1 let fnod=0 then fndoe and fq 's effects will be canceled out
! 2 do the subscribition after the elements loop
! Xingshi 09/15/2008
!===================================================
  do ie=1,ne		! loop over elements
     do inl=1,nen	
	     x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
!============================================================================
!		 fnode(1:nsd,inl) = f_fluids(1:nsd,ien(inl,ie))	
                 fnode(:,inl)=0.0
!============================================================================
		 d(1:ndf,inl) =  dloc(1:ndf,ien(inl,ie))
		 d_old(1:ndf,inl) = doloc(1:ndf,ien(inl,ie))
		f_stress(1:nsd,1:nsd,ien(inl,ie)) = 0.0
	 enddo
	ij = 1
	do inl=1,nen
	   do ii=1,ndf
		ddd(ij) = d(ii,inl)
		ij = ij+1
	   end do
	end do
	 hg = hk(ie)


	 
!!!!!!!!!!!!!!!!!!Intialize K and F!!!!!!!!!!!
  K(:,:) = 0.0
  F(:) = 0.0
  residual(:) = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 do iq=1,nquad  ! loop over the quadrature points in each element 
!...  calculate the shape function and the weight at quad point
		if (nsd==2) then
		    if (nen.eq.3) then !calculate shape function at quad point
			   include "sh2d3n.h"
			elseif (nen.eq.4) then
				include "sh2d4n.h"
			endif
		elseif (nsd==3) then
		    if (nen.eq.4) then !calculate shape function at quad point
			   include "sh3d4n.h"
			elseif (nen.eq.8) then
				include "sh3d8n.h"
			endif
		endif
!		stop
!	write(*,*) 'shape functions=',sh(0,1),sh(0,2),sh(0,3),sh(0,4)
	    eft0 = abs(det) * wq(iq) ! calculate the weight at each quad pt
!...  initialize d, dd/dx, dd/dy, dd/dz, dd/dt
		drs(:) = 0.0
		drt(:) = 0.0
		dr(1:nsd,1:ndf)=0.0
		fq(:)=0.0
		uold(:) = 0.0 !!initialize uold
!... calculate vi, dvi/dxj
        do inl=1,nen
		   tempc(1:nsd) = ama*d(1:nsd,inl)+oma*d_old(1:nsd,inl)
		   drs(1:nsd) = drs(1:nsd)+sh(0,inl)*tempc(1:nsd)
		   do isd=1,nsd
			 dr(isd,1:nsd) = dr(isd,1:nsd)+sh(isd,inl)*tempc(1:nsd)
		   enddo
		   fq(:) = fq(:) + sh(0,inl)*fnode(:,inl)        
	    enddo


!... calculate dvi/dt, p, dp/dxi
        do inl=1,nen
		   drt(1:nsd)=drt(1:nsd)+sh(0,inl)*(d(1:nsd,inl)-d_old(1:nsd,inl))*dtinv
		   uold(1:nsd)=uold(1:nsd)+sh(0,inl)*d_old(1:nsd,inl)*dtinv
		   drs(pdf)=drs(pdf)+sh(0,inl)*d(pdf,inl)    		   
		   dr(1:nsd,pdf)=dr(1:nsd,pdf)+sh(1:nsd,inl)*d(pdf,inl)       
	    enddo

!... define u=v1, v=v2, w=v3, pp=p
		if (nsd==2) then
		    u = drs(udf)
			v = drs(vdf)
			w = 0.0
		elseif (nsd==3) then
			u = drs(udf)
		    v = drs(vdf)
		    w = drs(wdf)
		endif
	    pp= drs(pdf)

	    if(stokes) then ! if stokes flow
		   u = 0.0
		   v = 0.0
		   w = 0.0
	    endif

!....  calculate liquid constant and gravity
	    mu = vis_liq  ! liquid viscosity
	    ro = den_liq  ! liquid density
		g  = gravity  ! gravatitional force

	! believe nu is calculated only for turbulent model
		if (nsd==2) then
	    nu = delta(4)*turb_kappa**2*hg**2 * sqrt(2*dr(1,1)**2+(dr(2,1)+dr(1,2))**2 &
	                                            +2*dr(2,2)**2)
		elseif (nsd==3) then
	    nu = delta(4)*turb_kappa**2*hg**2 * sqrt(2*dr(1,1)**2+(dr(2,1)+dr(1,2))**2 &
	                                            +2*dr(2,2)**2+(dr(3,1)+dr(1,3))**2 &
	                                            +2*dr(3,3)**2+(dr(3,2)+dr(2,3))**2)
		endif
	    mu = mu + nu*ro                   

!....  calculate each term in the residual equation
	    res_c = 0.0
		if (nsd==2) then
		  do inl=1,nen
		   	 res_c = res_c+sh(xsd,inl)*d(udf,inl) &
	                    +sh(ysd,inl)*d(vdf,inl)
		  enddo
		elseif (nsd==3) then
		  do inl=1,nen
		     res_c = res_c+sh(xsd,inl)*d(udf,inl) &
	                    +sh(ysd,inl)*d(vdf,inl) &
	                    +sh(zsd,inl)*d(wdf,inl)
		  enddo
		endif

	    do isd = 1, nsd
			if (nsd==2) then
			   res_a(isd)=ro*(drt(isd)+u*dr(1,isd)+v*dr(2,isd)-g(isd))-fq(isd)
			  
			elseif (nsd==3) then
			   res_a(isd)=ro*(drt(isd)+u*dr(1,isd)+v*dr(2,isd)+w*dr(3,isd)-g(isd))-fq(isd)
			endif
	    enddo

		res_t(1:nsd) = dr(1:nsd,pdf) + res_a(1:nsd)

		!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! Mickael 02/01/2005
        ! TAUm, TAUc and TAUl (l for  lsic), See Tezduyar and Sathe, 2003

                                
        vel  = sqrt(u*u+v*v+w*w)  !magnitude of the velocity
                        
        ree  = vel*hg*ro/(2.0*mu)  !Reynolds number
        if(steady.or.(.not.taudt)) then !stablization, taum
			taum = 1.0/sqrt((2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
		else
			taum = 1.0/sqrt((2.0/dt)**2+(2.0*vel/hg)**2+(4.0*mu/hg/hg)**2)
        endif
        tauc = taum/ro
                        
        if (ree.le.3.0) then
			taul = hg*vel*ree/6.0
        else
            taul = hg*vel/2.0
        endif
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!.....   Density optimization

		ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0
	      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c.....   Galerkin Terms (Look at notes)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c....   Calculate stress tau(ij) and pressure
		do isd = 1,nsd
			do jsd = 1,nsd
				tau(isd,jsd) = mu*(dr(isd,jsd) + dr(jsd,isd))
				do inl=1,nen	
					! Compute fluid stress
					f_stress(isd,jsd,ien(inl,ie))= f_stress(isd,jsd,ien(inl,ie)) + tau(isd,jsd)
				enddo
			enddo
		enddo

		prs_t(1:nsd) = res_t(1:nsd)*taum


		!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! Mickael 02/01/2005
        ! TAUl (l for  lsic), See Tezduyar and Sathe, 2003
  
        prs_c = ro*res_c*taul
        !
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	    prs_cc(1:nsd) = res_t(1:nsd)*tauc	      
!.... calculate the residual at each degree of freedom
	    do inl=1,nen ! loop over number of nodes in an element
		 
		   node = ien(inl,ie)
			if (nsd==2) then
			   temp = ro*(u*ph(xsd,inl)+v*ph(ysd,inl))
			elseif (nsd==3) then
			   temp = ro*(u*ph(xsd,inl)+v*ph(ysd,inl)+w*ph(zsd,inl))
			endif
		! Continuty Equation
		   p(pdf,node) = p(pdf,node)-ph(0,inl)*res_c
		   res_orig(3+3*(inl-1))=p(pdf,node)  !!!original residual for p
		! Momentum Equation (Euler Residual)
!===========================================================
! why originally it is minus?
		   p(1:nsd,node) = p(1:nsd,node)-ph(0,inl)*res_a(1:nsd)
!=============================================================
		! Momentum Equation (C: Diffusion terms)
		    if (nsd==2) then
			  do isd=1,nsd
				p(isd,node)=p(isd,node) + ph(isd,inl)*pp -   &
										  ph(1,inl)*tau(1,isd) -  &
										  ph(2,inl)*tau(2,isd)
			  enddo
			elseif (nsd==3) then
			  do isd=1,nsd
				p(isd,node)=p(isd,node) + ph(isd,inl)*pp -   &
										  ph(1,inl)*tau(1,isd) -  &
										  ph(2,inl)*tau(2,isd) -  &
										  ph(3,inl)*tau(3,isd)
			  enddo
			endif
		   res_orig(1+3*(inl-1))=p(1,node)
		   res_orig(2+3*(inl-1))=p(2,node)

!		! Stablization with Tau_moment
!		   if (nsd==2) then
!		   	 p(pdf,node) = p(pdf,node) - ph(xsd,inl)*prs_cc(udf)  &
!	                                   - ph(ysd,inl)*prs_cc(vdf)
!	       elseif (nsd==3) then
!		     p(pdf,node) = p(pdf,node) - ph(xsd,inl)*prs_cc(udf)  &
!	                                   - ph(ysd,inl)*prs_cc(vdf)  &
!	                                   - ph(zsd,inl)*prs_cc(wdf)
!		   endif		! Stablization with Tau_cont    
!		   p(1:nsd,node) = p(1:nsd,node) - prs_t(1:nsd)*temp - ph(1:nsd,inl)*prs_c
	 enddo
!********************************************************
              ! Stiffness matrix and source term
!********************************************************

  do ii = 1,nen
     do ij = 1,nen
!!!!!!!!!!!!!!!!K for x direction momentum equation
	K(1+3*(ii-1),1+3*(ij-1))=K(1+3*(ii-1),1+3*(ij-1))+(ro*sh(0,ii)*sh(0,ij)*dtinv+ro*sh(0,ii)*u*sh(1,ij)+ &
ro*sh(0,ii)*v*sh(2,ij)+mu*(sh(1,ii)*(sh(1,ij)+sh(1,ij))+sh(2,ii)*sh(2,ij)))*eft0
	K(1+3*(ii-1),2+3*(ij-1))=K(1+3*(ii-1),2+3*(ij-1))+mu*sh(2,ii)*sh(1,ij)*eft0
	K(1+3*(ii-1),3+3*(ij-1))=K(1+3*(ii-1),3+3*(ij-1))-sh(1,ii)*sh(0,ij)*eft0 
!!!!!!!!!!!!!!!1K for y direction momentum equation
	K(2+3*(ii-1),1+3*(ij-1))=K(2+3*(ii-1),1+3*(ij-1))+mu*sh(1,ii)*sh(2,ij)*eft0
	K(2+3*(ii-1),2+3*(ij-1))=K(2+3*(ii-1),2+3*(ij-1))+(ro*sh(0,ii)*sh(0,ij)*dtinv+ro*sh(0,ii)*u*sh(1,ij)+ &
ro*sh(0,ii)*v*sh(2,ij)+mu*(sh(2,ii)*(sh(2,ij)+sh(2,ij))+sh(1,ii)*sh(1,ij)))*eft0
	K(2+3*(ii-1),3+3*(ij-1))=K(2+3*(ii-1),3+3*(ij-1))-sh(2,ii)*sh(0,ij)*eft0
!!!!!!!!!!!!!!!!K for continuity equation
	K(3+3*(ii-1),1+3*(ij-1))=K(3+3*(ii-1),1+3*(ij-1))+sh(0,ii)*sh(1,ij)*eft0
	K(3+3*(ii-1),2+3*(ij-1))=K(3+3*(ii-1),2+3*(ij-1))+sh(0,ii)*sh(2,ij)*eft0
	K(3+3*(ii-1),3+3*(ij-1))=K(3+3*(ii-1),3+3*(ij-1))+0
     end do
!!!!!!!!!!!!!!!!Calculate F term!!!!!!!!!!!!!!!!!!!
	F(1+3*(ii-1))=F(1+3*(ii-1))+ro*(sh(0,ii)*(uold(1)+g(1)))*eft0
	F(2+3*(ii-1))=F(2+3*(ii-1))+ro*(sh(0,ii)*(uold(2)+g(2)))*eft0
	F(3+3*(ii-1))=F(3+3*(ii-1))+0
  end do

!********************************************************
	      ! Diagonal Preconditioner
!********************************************************
	    effd = mu*eft0*alpha
	    effm = taum*eft0
	    effc = tauc*eft0

	    do inl=1,nen

		   node=ien(inl,ie)
		   if (nsd==2) then
			 ug = ro*(u*sh(1,inl)+v*sh(2,inl))
		   elseif (nsd==3) then
		     ug = ro*(u*sh(1,inl)+v*sh(2,inl)+w*sh(3,inl))
		   endif
		   temp = alpha*ug + sh(0,inl)*dtinv*ro
		 
		   do isd=1,nsd
			 q(isd,node) = q(isd,node) + (alpha*ro*sh(0,inl)*dr(isd,isd)+temp)*(sh(0,inl)*eft0+ug*effm)
		   enddo
		   if (nsd==2) then
		   	 temp = sh(1,inl)**2 + sh(2,inl)**2
			elseif (nsd==3) then
				temp = sh(1,inl)**2 + sh(2,inl)**2 + sh(3,inl)**2
			endif

		   do isd=1,nsd
		     q(isd,node) = q(isd,node) + (sh(isd,inl)**2+temp)*effd
		   enddo

		   q(pdf,node) = q(pdf,node)+temp*effm
		   q(1:nsd,node) = q(1:nsd,node)+sh(1:nsd,inl)**2*effc
	    enddo

	 enddo ! end of qudrature pts loop

	do ii = 1,12
	   do ij = 1,12
	      residual(ii) = residual(ii)+K(ii,ij)*ddd(ij)
	   end do
	   residual(ii) = residual(ii)-F(ii)
	end do
!write(*,*)'index=',ien(:,ie)
!write(*,*)'diag='
!        do ii=1,12
!	   do ij =1,12
!write(*,*)K(ii,ij)
!	   end do
!        end do
!write(*,*)'K21=',K(2,1)
!write(*,*)'K2,10=',K(2,10)
!write(*,*)'F=',F(:)
!write(*,*)'ddd=',ddd(:)
!write(*,*)'residual = ',residual(:)
!write(*,*)'sh(1,1)=',sh(1,4)
!write(*,*)'sh(2,1)=',sh(2,4)
!write(*,*)'orginal residual=', res_orig(:)
	do ii=1,nen
	   inode=ien(ii,ie)
	   do ij=1,nen
	      jnode =ien(ij,ie)
	      do icount=1,3
	         do jcount=1,3

	             KK(icount+3*(inode-1),jcount+3*(jnode-1))=KK(icount+3*(inode-1),jcount+3*(jnode-1))+ &
			 K(icount+3*(ii-1),jcount+3*(ij-1))
	         end do
	      end do
	   end do
	end do




  enddo ! end of element loop
!write(*,*)'KK(:,:)='
open(11,file='matrixK.txt',status='unknown')
do ii=1,ndf*nn
   do ij=1,ndf*nn
      write(11,*) KK(ii,ij)
   end do
enddo
close(11)
 continue  
continue

  return
end subroutine block

