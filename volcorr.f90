      subroutine volcorr(x_solid,nsd,nen,ien,nn_solid,ne_solid,x_solid_init)
! x_solid: the potision of the solid domain at current time step
! nsd: dimension of the problem
! nen: number of nodes per element of solid mesh
! g: correction of the updated solid points
! ien: connectiviy of the solid mesh
! nn_solid: number of nodes of solid mesh
! ne_solid: number of elements of solid mesh
! e_area: the area of every element
! rhs: (det(F_{ij}))^{-1}-1 right hand side of the equations 


      integer nsd
	  integer nn_solid
	  integer ne_solid
	  integer ien(ne_solid,nen)
      real(8) x_solid(nsd,nn_solid)
	  real(8) g(nsd,nn_solid)

	  real(8) x_solid_init(nsd,nn_solid)
!=============================================================
! shape funciton variable       
	  integer iquad
	  integer nquad
	  real(8) xq(nsd,1)
	  real(8) wq
	  integer iq
	  real(8) sq(0:nsd,nen,1)
	  real(8) sh(0:nsd,nen)
	  real(8) xr(nsd,nsd)
	  real(8) cf(nsd,nsd)
          real(8) sx(3,3)
	  real(8) det
! iquad: number of integration points always 1 here
! nquad: always 1
! sq: the same meaning of in the main program but only 1 in 3nd dimension here
! iq: always 1 here just for convenience
!===============================================================
! variables to calculate rhs & e_area
      real(8) e_area(ne_solid)
	  real(8) rhs(ne_solid)
	  real(8) y(nsd,nen)
	  real(8) rs(nsd)
	  real(8) xj(nsd,nsd)
	  real(8) xji(nsd,nsd)
	  real(8) detJ
	  real(8) toxj(nsd,nsd)
	  real(8) toxji(nsd,nsd)
	  real(8) todet
	  real(8) xto(nsd,nsd)
	  real(8) xot(nsd,nsd)
	  real(8) toc(nsd,nsd)
	  real(8) detf
! det: determinant of Jacobian at current time step
! detf: det(F_{ij}^{-1})
!=======================================
! Seting up left and right hand side variables
      integer ine
	  integer inen
	  integer index
	  integer jne
	  integer jnen
	  integer indey
	  integer inn
	  real(8) x(nsd,nen)
!	  real(8) y(nsd,nen)
	  real(8) crrm(nsd*nn_solid,nsd*nn_solid)
	  real(8) crra(nsd*nn_solid)
	  real(8) s_area(nn_solid)
! crrm: the whole matirx of the left hand side 		
! crra: the whole matrix of the right hand side
! s_area: the sum area of the elements sounding every point
!====================================================
! varialbe to solve the linear system
      integer np
      integer mp
      integer error
      real(8) crra_cg(nsd*nn_solid)
!=============================================================
! First step: calculate the shape funciton in the parent domain, based on nsd & nen      
	  iquad=1
	  iq=1
          xq(:,:)=0.0d0
	  if (nsd==2) then
	     if (nen==3) then
 	!              write(*,*) 'I am here in volcorr 84', xq(:,:)
		     call quad2d3n(iquad,nquad,xq,wq,nsd,1)
         !             write(*,*) 'I am after call quad', xq(:,:)
		 else if (nen==4) then
		     call quad2d4n(iquad,nquad,xq,wq,nsd,1)
		 end if 
		 if(nen==3) then
		  sq(0,1,iq) = xq(1,iq)
		  sq(0,2,iq) = xq(2,iq)
		  sq(0,3,iq) = 1 - xq(1,iq) - xq(2,iq)
        else if (nen==4) then
		  sq(0,1,iq) = (1 - xq(1,iq)) * (1 - xq(2,iq)) / 4
		  sq(0,2,iq) = (1 + xq(1,iq)) * (1 - xq(2,iq)) / 4
		  sq(0,3,iq) = (1 + xq(1,iq)) * (1 + xq(2,iq)) / 4
		  sq(0,4,iq) = (1 - xq(1,iq)) * (1 + xq(2,iq)) / 4
		  
		  sq(1,1,iq) = - (1 - xq(2,iq)) / 4
		  sq(1,2,iq) = + (1 - xq(2,iq)) / 4
		  sq(1,3,iq) = + (1 + xq(2,iq)) / 4
		  sq(1,4,iq) = - (1 + xq(2,iq)) / 4
	
		  sq(2,1,iq) = - (1 - xq(1,iq)) / 4
		  sq(2,2,iq) = - (1 + xq(1,iq)) / 4
		  sq(2,3,iq) = + (1 + xq(1,iq)) / 4
		  sq(2,4,iq) = + (1 - xq(1,iq)) / 4

        endif
	  else if (nsd==3) then
	     if (nen==4) then 
		     call quad3d4n(iquad,nquad,xq,wq,nsd,1)
		 else if  (nen==8) then
		     call quad3d8n(iquad,nquad,xq,wq,nsd,1)
		 end if
	  
		 if(nen.eq.4) then
		  sq(0,1,iq) = xq(1,iq)
		  sq(0,2,iq) = xq(2,iq)
		  sq(0,3,iq) = xq(3,iq)
		  sq(0,4,iq) = 1 - xq(1,iq) - xq(2,iq) - xq(3,iq)
        else
		  sq(0,1,iq) = (1 - xq(1,iq))	&
			   * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,2,iq) = (1 + xq(1,iq))	&
			   * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,3,iq) = (1 + xq(1,iq))	&
			   * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,4,iq) = (1 - xq(1,iq))	&
			   * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(0,5,iq) = (1 - xq(1,iq))	&
			   * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(0,6,iq) = (1 + xq(1,iq))	&
			   * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(0,7,iq) = (1 + xq(1,iq))	&
			   * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(0,8,iq) = (1 - xq(1,iq))	&
			   * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,1,iq) = - (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,2,iq) = + (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,3,iq) = + (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,4,iq) = - (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
		  sq(1,5,iq) = - (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,6,iq) = + (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,7,iq) = + (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(1,8,iq) = - (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,1,iq) = - (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,2,iq) = - (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,3,iq) = + (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,4,iq) = + (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
		  sq(2,5,iq) = - (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,6,iq) = - (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,7,iq) = + (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(2,8,iq) = + (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
		  sq(3,1,iq) = - (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,2,iq) = - (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,3,iq) = - (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
		  sq(3,4,iq) = - (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
		  sq(3,5,iq) = + (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,6,iq) = + (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
		  sq(3,7,iq) = + (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
		  sq(3,8,iq) = + (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
        endif
	  end if

!==================================================================
! Second step: Set up the left and right handside matrix 
      crra(:)=0
	  crrm(:,:)=0
	  s_area(:)=0 

! 2.0 calculate deformation matrix at 1 gau point and area per element
      do ine=1,ne_solid
	     do inen=1,nen
		    y(1:nsd,inen)=x_solid(1:nsd,ien(ine,inen)) ! current position
			x(1:nsd,inen)=x_solid_init(1:nsd,ien(ine,inen)) ! initial position
		 end do
		 rs(1:nsd)=xq(1:nsd,1)
		 call r_element(rs)
		 call r_jacob(y,xj,xji,detJ)
		 call r_jacob(x,toxj,toxji,todet)
		 call r_bdpd_curr(xji)
		 call r_bdpd_init(toxji)
		 call r_stoxc(xto,xot,xj,xji,toxj,toxji,toc,ine)
		 e_area(ine)=wq*detJ
		 call determinant(xot,nsd,nsd,detf)
		 rhs(ine)=detf-1
	  end do
       continue


!loop over element
!==========================================================
! NOTE
! Alougth the system of linear equations is derived from per nodes,
! it is more convenient to calculate the terms in both right and left handside
! matrix from element to elment.
! The formation may be tricky and need to be tested. Also, how I break the coefficients 
! of this linear system should be write down in the notes.
!==========================================================
      
	  
	  do ine=1,ne_solid
! 2.1: get N_{j,i} 
         do inen=1,nen
		    x(:,inen)=x_solid(1:nsd,ien(ine,inen))
		 end do  
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
! 2.2: loop over nsd and elements to set up the left and rigt handside matrix
		do isd=1,nsd
		do inen=1,nen
		   index=ien(ine,inen)+(isd-1)*nn_solid
		   do jsd=1,nsd
		   do jnen=1,nen
		   indey=ien(ine,jnen)+(jsd-1)*nn_solid
           crrm(index,indey)=crrm(index,indey)+sh(isd,inen)*sh(jsd,jnen)
! left handside matrix but missing 1/A_{\alpha} for the diagonal terms

		   end do
		   end do
		   crra(index)=crra(index)+rhs(ine)*sh(isd,inen)

		end do
		end do
		end do
		do ine=1,ne_solid
		   do inen=1,nen  
! right handside matrix
		   s_area(ien(ine,inen))=s_area(ien(ine,inen))+e_area(ine)
! s_area is the area sum of the sounding elements
           end do
		end do

          continue
        do inn=1,nn_solid
		   do isd=1,nsd
		      index=inn+(isd-1)*nn_solid
			  crrm(index,index)=crrm(index,index)+1.0/s_area(inn)
! adding 1/A_{\alpha} to the diagonal terms
		   end do
		end do
! Finish setting up left and right handside matrix
!====================================================================
! Step 3 solve the linear equations & update the correction
! For testing just use gausian
!        open (unit=100, file='a.txt', status='unknown') 
!        open (unit=101, file='b.txt', status='unknown')
!        do ine=1,nsd*nn_solid
!           do inen=1,nsd*nn_solid
!               write(100,*) crrm(ine,inen)
!           end do
!           write(101,*) crra(ine)
!        end do
        np=nsd*nn_solid
        mp=1
        continue
!=======================================================
!We have two options to solve linear system eqautions
! Gausian
!        call gaussj(crrm,np,np,crra,np,mp)
! CG
         error=0
         crra_cg(:)=0.0d0
         call cg_method(np,crrm,crra,crra_cg,error)
		do isd=1,nsd
		   do inn=1,nn_solid
		      index=inn+(isd-1)*nn_solid
		      !g(isd,inn)=crra(index) !For guasian
                      g(isd,inn)=crra_cg(index) !For CG
			  x_solid(isd,inn)=x_solid(isd,inn)+g(isd,inn)
			  end do
		end do
		return 
		end
