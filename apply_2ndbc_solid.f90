subroutine apply_2ndbc_solid(x_solid,nsd,nn_solid,ien_sbc,ne_sbc,nen_solid,&
                        ien_solid,ne_solid,solid_bcforce,solid_stress)
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
real(8) solid_bcforce(nsd,nn_solid)
real(8) solid_stress(nsd*2,nn_solid)
!------------------------------------
integer count
integer ibs
integer ine
integer ntem
integer nos
real(8), allocatable, dimension(:,:) :: x
real(8) p(nsd)
real(8) w
real(8) norm(nsd)
real(8) area
integer isd
integer jsd
real(8) stress_tmp(nsd,nsd)
integer error_id

solid_bcforce(:,:)=0.0d0


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
		! find the local nodes on the 2nd type BC surface
if (count .ne. 4) write(*,*) 'Wrong !!!!!',ibs
		if (nen_solid == 4) then
			call outnormal_tet(x,nsd,p,norm,area)
		else if (nen_solid == 8) then
			call outnormal_hex(x,nsd,p,norm,area)
		end if
		! get the ourward norm

                do nos=1,nen_solid
                ntem=ien_solid(ine,nos)

                if (ien_sbc(ibs,nos+1) == 1 .and. ien_sbc(ibs,nen_solid+2) == -999) then
                stress_tmp(1,1)=solid_stress(1,ntem)
                stress_tmp(2,2)=solid_stress(2,ntem)
                stress_tmp(3,3)=solid_stress(3,ntem)
                stress_tmp(2,3)=solid_stress(4,ntem)
                stress_tmp(3,2)=solid_stress(4,ntem)
                stress_tmp(1,3)=solid_stress(5,ntem)
                stress_tmp(3,1)=solid_stress(5,ntem)
                stress_tmp(1,2)=solid_stress(6,ntem)
                stress_tmp(2,1)=solid_stress(6,ntem)

			do isd=1,nsd
				do jsd=1,nsd
				solid_bcforce(isd,ntem)=solid_bcforce(isd,ntem)+&
				area*stress_tmp(isd,jsd)*norm(jsd)*w
				end do
			end do

		end if
		end do
		! integrate over 2nd BC surface \sigma_{ij} n_{j} w 
continue

end do


deallocate(x)

continue
return
end 
