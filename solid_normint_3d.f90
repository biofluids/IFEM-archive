subroutine solid_normint_3d(x_solid,nsd,nn_solid,ien_sbc,ne_sbc,nen_solid,&
                        ien_solid,ne_solid,solid_norm)
!---------------------------------
! Calculate solid surface normal integral for 3-D cases
! 1 find the points on the boundary surface
! 2 use two vector cross to find the norm of the norm of the surface
! 3 take one interial point and take the dot between this vector and norm 
!   if dot < 0 it is the outward norm and then normalize the norm
! 4 calculate the area of boundary surface: (1) tet, area of triangle; (2) hex, 2 triangles
! 5 take norm*area and distribute it equaly to the each node

implicit none

real(8) x_solid(nsd,nn_solid)
integer nsd
integer nn_solid
integer ien_sbc(ne_sbc,nen_solid+2)
integer ne_sbc
integer nen_solid
integer ien_solid(ne_solid,nen_solid)
integer ne_solid
real(8) solid_norm(nsd,nn_solid)
!------------------------------------
integer count
integer ibs
integer ine
integer ntem
integer nos
real(8), allocatable :: x(:,:)
real(8) p(nsd)
real(8) w
real(8) norm(nsd)
real(8) area
real(8) tot_area
solid_norm(:,:)=0.0d0
tot_area=0.0d0


if (nen_solid == 4) then
allocate (x(nsd,3))
w=1.0/3.0
end if

if (nen_solid == 8) then
allocate (x(nsd,4))
w=1.0/4.0
end if


do ibs=1,ne_sbc
	ine=ien_sbc(ibs,1)
	count=1
		do nos=1,nen_solid
		ntem=ien_solid(ine,nos)
		if (ien_sbc(ibs,nos+1) == 1) then
		x(1:nsd,count)=x_solid(1:nsd,ntem)
		count=count+1
		else
		p(1:nsd)=x_solid(1:nsd,ntem)
		end if
		end do

		if (nen_solid == 4) then
			call outnormal_tet(x,nsd,p,norm,area)
		else if (nen_solid == 8) then
			call outnormal_hex(x,nsd,p,norm,area)
		end if

		tot_area=area+tot_area
                do nos=1,nen_solid
                ntem=ien_solid(ine,nos)
                if (ien_sbc(ibs,nos+1) == 1) then
		solid_norm(:,ntem)=solid_norm(:,ntem)+area*norm(:)*w
		end if
		end do
continue

end do

write(*,*) 'tot surface area', tot_area

deallocate(x)

end 
