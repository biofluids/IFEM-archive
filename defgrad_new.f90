subroutine defgrad_new(xloc,xrefloc,ien,f,jac,finv,jacinv)
! Calculate deformation gradient for ALE by Xingshi Wang
        use fluid_variables, only: nen,ne,nn,nsd,nquad,wq,sq
        implicit none

        integer ien(nen,ne)
        real* 8 xloc(nsd,nn), xrefloc(nsd,nn)
        real* 8 x(nsd,nen)
	real*8  f(nsd,nsd,nquad,ne) ! F
        real* 8 finv(nsd,nsd,nquad,ne) ! F^{-1]
        real* 8 sh(0:nsd,nen), ph(0:nsd,nen)
        real* 8 eft0,det
	real*8  jac(nquad,ne), jacinv(nquad,ne) ! det(F) and det(F^{-1})
	real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)
	integer inl, ie, isd, jsd, ksd, iq, node
	real*8  ftmp(nsd,nsd)
	real*8  ftmp1(nsd,nsd)
	real*8  check(nsd,nsd)
	real*8  jactmp
	integer info
	integer ipiv(nsd)
	real*8  work(nsd)



	! Get shape function 
   do ie=1,ne
	   do inl=1,nen
              do isd=1,nsd 
                 x(isd,inl)=xrefloc(isd,ien(inl,ie))
              enddo
           enddo

      do iq=1,nquad

	if (nsd == 3) then ! 3-D case
              if (nen.eq.4) then
                 include "sh3d4n.h"
              else if (nen.eq.8) then
                 include "sh3d8n.h"
              end if	
	end if

        if (nsd==2) then ! 2-D case
              if (nen.eq.3) then 
                  include "sh2d3n.h"
              elseif (nen.eq.4) then
                  include "sh2d4n.h"
              endif	
	end if

	ftmp(:,:)=0.0d0
	do inl=1,nen
	   node=ien(inl,ie)
		do isd=1,nsd
			do jsd=1,nsd
			   ftmp(isd,jsd)=ftmp(isd,jsd)+sh(jsd,inl)*xloc(isd,node) ! F= d x_i / d xref_j
			end do
		end do
	end do
	
	f(:,:,iq,ie)=ftmp(:,:)
	call determinant(ftmp,nsd,nsd,jactmp) ! det(F)
	jac(iq,ie)=jactmp

!ftmp1(:,:)=ftmp(:,:)

	call DGETRF(nsd,nsd,ftmp,nsd,ipiv,info)
	call DGETRI(nsd,ftmp,nsd,ipiv,work,nsd,info) ! get F^{-1]
	finv(:,:,iq,ie)=ftmp(:,:)
	call determinant(ftmp,nsd,nsd,jactmp) ! det(F^{-1})
	jacinv(iq,ie)=jactmp
! check F * F{-1} = I	
!check(:,:)=0.0d0
!do isd=1,nsd
!	do jsd=1,nsd
!		do ksd=1,nsd
!		check(isd,jsd)=check(isd,jsd)+ftmp1(isd,ksd)*ftmp(ksd,jsd)
!		end do
!	end do
!end do

! write(*,*) 'check', check(:,:)




     end do ! Gausian point loop
   end do ! element loop


return
end
