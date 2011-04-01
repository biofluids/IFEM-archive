subroutine nature_pre_3d(xloc,dloc,ien,rng,bdel,ne_inflow,bdindex,pin,res_bc)
  use global_constants
  use fluid_variables, only: nn,nsd,nen,ne,ndf,neface,vis_liq,den_liq,mapping
implicit none

! input arguments
 real(8) xloc(nsd,nn) ! coordinates of fluid nodes
 real(8) dloc(ndf,nn) ! [u,v,p] for 2-D [u,v,w,p] for 3d
 integer ien(nen,ne)  ! connectivity matrix
 integer rng(neface,ne) ! all boundary information
 integer bdel(ne_inflow) ! element inndex for those on the inflow boundary
 integer ne_inflow ! number of elements on the inflow boundary
 integer bdindex ! set which edge has inflow boundary condition in the whole fluid domain
 real(8) pin ! pressure difference for inflow
 real(8) res_bc(nsd,nn) ! residuale coming from nature B.C. integration
! local variables
integer ie
integer inl
real(8) d(ndf,nen)
integer iface
integer ie_inflow
integer n_edge ! number of nodes per edge
! loop variables
integer ia
integer ix
integer tmp
real(8) p(nsd)
real(8) w
integer isd
real(8) norm(nsd)
real(8) x(nsd,4)
integer li(4) ! local node index who sits on the boundary in each element
real(8) area
real(8) face_area
!=======================================
if (nen == 4) then
n_edge=3 ! set number of nodes per edge NOTE: only for 2-D case right now !!!!
w=1.0/3.0
end if

if (nen == 8) then
n_edge=4
w=1.0/4.0
end if
res_bc(:,:)=0.0d0
li(:)=0
x(:,:)=0.0d0
face_area=0.0d0
!===================================================
  do ie_inflow=1,ne_inflow            ! loop over elements on the inflow boundary
     ie=bdel(ie_inflow)               ! global element index
!===========================================================
     do iface=1,neface
        ! decide the knowns and unknows based on iface
        if (rng(iface,ie) == bdindex) then
		if (nen == 4) then
                        li(1:n_edge)=mapping(iface,1:n_edge,3) ! tet case

!			li(1)=mapping(iface,1,3) ! tet case
!                        li(2)=mapping(iface,2,3) ! tet case
!                        li(3)=mapping(iface,3,3) ! tet case

		end if
		
		if (nen == 8) then
			li(1:n_edge)=mapping(iface,1:n_edge,4) ! hex case
		end if
        end if
     end do ! for iface

	do ia=1,n_edge
		ix=ien(li(ia),ie)
		x(1:nsd,ia)=xloc(1:nsd,ix)
	end do

	do inl=1,nen
		do ia=1,n_edge
			if (inl .ne. li(ia)) tmp=inl
		end do
	end do
	p(1:nsd) = xloc(1:nsd,ien(tmp,ie))

                if (nen == 4) then
                        call outnormal_tet(x(1:nsd,1:n_edge),nsd,p,norm,area)
                else if (nen == 8) then
                        call outnormal_hex(x(1:nsd,1:n_edge),nsd,p,norm,area)
                end if
	norm(:)=-norm(:)
	continue
!=========================
	do ia = 1,n_edge
	ix=ien(li(ia),ie)
		do isd=1,nsd
		res_bc(isd,ix) = res_bc(isd,ix) + pin*area*w*norm(isd) ! force dot with norm
		end do
	end do
!=========================
	face_area=area+face_area
end do
!deallocate(x)
!deallocate(li)


return
end
       


