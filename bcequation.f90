subroutine bcequation(xloc,dloc,ien,ginflow,linflow,nn_inflow,rng,bdel,ne_inflow,bdindex)
! This subroutine is used to solve the boundary equation
! Returen u,v at inflow boundary equation and global error for boundary equation
  use global_constants
!  use run_variables
  use fluid_variables, only: nn,nsd,nen,ne,ndf,neface,vis_liq,den_liq
! input arguments
 real(8) xloc(nsd,nn) ! coordinates of fluid nodes
 real(8) dloc(ndf,nn) ! [u,v,p] for 2-D [u,v,w,p] for 3d
 integer ien(nen,ne)  ! connectivity matrix
 integer ginflow(nn)  ! global index for nodes on inflow  boundary
 integer linflow(nn_inflow) ! boundary equation index for nodes on inflow boundary
 integer nn_inflow    ! # of nodes on inflow boundary
 integer rng(neface,ne) ! all boundary information
 integer bdel(ne_inflow) ! element inndex for those on the inflow boundary
 integer ne_inflow ! number of elements on the inflow boundary
 integer bdindex ! set which edge has inflow boundary condition in the whole fluid domain
!===========================================
! variable for shape function
  real(8) x(nsd,nen)
  real(8) det
  real(8) sh(0:nsd,nen),ph(0:nsd,nen)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)
  real(8) sq(0:nsd,nen,nen) ! This sq is not the same in fluid_variables
			! In fluid_variables, sq is the shape function for integration points
			! Here sq is the shape funtion for the element nodes
!===================================
! Loop variables
integer ie
integer inl
real(8) d(ndf,nen)
integer iface
!===================================
! Boundary equation variables
integer ie_inflow
integer li(2) ! the local index of unknowns only for 2-D case
integer n_edge ! number of nodes per edge
real(8) pin ! pressure difference for inflow
real(8) fx  ! RHS for boundary equation in x
real(8) fy ! RHS for boundary eqaution in y
real(8) a_lhs(nn_inflow*nsd,nn_inflow*nsd) ! LHS matrix A 
real(8) f_rhs(nn_inflow*nsd) ! RHS vector F
integer neq ! counter of # of equation obtained from looping the elements containing boundary edge
integer ia
integer ja
real(8) el_len
integer jx
integer jnl
real(8) inv_a(nn_inflow*nsd,nn_inflow*nsd)
integer indx(nn_inflow*nsd)
!======================================
n_edge=2 ! set number of nodes per edge NOTE: only for 2-D case right now !!!!
pin=10 ! set the pressure difference between inflow and outflow, finally it should be set in input_fluid.in
neq=0 ! set neq to be 0, as beginning
a_lhs(:,:)=0.0d0 ! initial matrix
f_rhs(:)=0.0d0 ! initial vector
!===================================================
  do ie_inflow=1,ne_inflow            ! loop over elements on the inflow boundary
     ie=bdel(ie_inflow)               ! global element index
     do iface=1,neface
     
     	do inl=1,nen       
             x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
	     d(1:ndf,inl) =  dloc(1:ndf,ien(inl,ie))
     	end do
!====================================================================================================
        ! decide the knowns and unknows based on iface
	if (rng(iface,ie) == bdindex) then
		if (nsd==2) then ! only for 2-D case right now
			if (nen==4) then ! Quadratrial case
					if (iface == 1) then ! case I local nodes 1,2 on the bounndary
						li(1)=1
						li(2)=2
						iq=1
					end if
					if (iface == 2) then ! case II local nodes 2,3 on the boundary
                                                li(1)=2
                                                li(2)=3
						iq=2
                                        end if
					if (iface == 3) then ! case III local nodes 3,4 on the boundary
                                                li(1)=3
                                                li(2)=4
						iq=3
                                        end if
					if (iface == 4) then ! case IV local nodes 4,1 on the boundary
                                                li(1)=1
                                                li(2)=4
						iq=4
                                        end if
				! These four cases are decided based on how we define mrng.in

			end if
			if (nen==3) then ! Triangle case
					if (iface == 1) then ! case I local nodes 1,2 on the bounndary
                                                li(1)=1
                                                li(2)=2
						iq=1
                                        end if
                                        if (iface == 2) then ! case II local nodes 2,3 on the boundary
                                                li(1)=2
                                                li(2)=3
						iq=2
                                        end if
                                        if (iface == 3) then ! case III local nodes 3,4 on the boundary
                                                li(1)=1
                                                li(2)=3
						iq=3
                                        end if
			end if
				! These three cases are decided based on how we define mrng.in
		end if
	end if
  end do ! for iface
!======================================================================================================
! Calculate shape function
! iq is the face index of lacal element
call bc_shape2d(sq,nsd,nen,neface) ! get shape function in local sapce
     ! there is no need for more than one integration point
                if (nsd==2) then
                    if (nen.eq.3) then !calculate shape function at quad point^M
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
                endif ! shape function in global space for node iq
!======================================================================================================
! Decide RHS --- set fx and fy
! Note right now we just assume the normal direction is x and tangial direction is y
fx=0.0d0
fy=0.0d0
write(*,*) 'N,x', sh(1,:)
write(*,*) 'N,y', sh(2,:)
do inl=1,nen
	fx=fx+sh(1,inl)*d(1,inl) ! fx= sum{N_i,x * u_i}
end do
do inl=1,n_edge
fx=fx-sh(1,li(inl))*d(1,inl) ! fx=N_3,x * u_3 + N_4,x * u_4 (See notes)
end do
fx=pin/(vis_liq/den_liq)-fx ! fx= P/(nu) - N_3,x * u_3 + N_4,x * u_4 (Also see notes)
! fx has been given till now

do inl=1,nen
	fy=fy+sh(1,inl)*d(2,inl)+sh(2,inl)*d(1,inl) ! fy=sum(N_i,x * v_i)+sum{N_i,y * u_i}
end do
do inl=1,n_edge
fy=fy-sh(1,li(inl))*d(2,li(inl))+sh(2,li(inl))*d(1,li(inl)) ! fy= N_3,y*u_3 + N_4,y * u_4 + N_3,x * v_3 + N_4.x * v_4
end do
fy=-fy
! fy has been given till now
! Assume that only use one integration point in line integration so M_1 = M_2 =1/2 check notes
!=========================================================================================================
! Sort of assembling produce set LHS matrix and RHS matrix using ginflow, mien to set up the linear system of equation
el_len=(x(1,li(1))-x(1,li(2)))**2 + (x(2,li(1))-x(2,li(2)))**2
el_len=sqrt(el_len) ! length of the line element, Only work for 2-D case !!!!
do inl=1,n_edge
ix=ien(li(inl),ie) ! global node index for the node on the boundary
ix=ginflow(ix) ! global node index to local boundary equation index
do jnl=1,n_edge
jx=ien(li(jnl),ie)
jx=ginflow(jx)
a_lhs(ix*2-1,2*jx-1)=a_lhs(ix*2-1,jx*2-1)+sh(1,li(jnl))*el_len ! equation in x direction
a_lhs(ix*2,2*jx-1)=a_lhs(ix*2,2*jx-1)+sh(2,li(jnl))*el_len  ! equation in y direction 
a_lhs(ix*2,2*jx)=a_lhs(ix*2,2*jx)+sh(1,li(jnl))*el_len
end do
f_rhs(ix*2-1)=f_rhs(ix*2-1)+fx*el_len
f_rhs(ix*2)=f_rhs(ix*2)+fy*el_len
end do
! We are done with setting A and F
!=========================================================================================================
!  end do ! for iface loop
end do ! for ie_inflow loop
!==================================================
! For this 1-D boundary equation, essential boundary condtion u=0=v at two ends should be applied
a_lhs(1,:)=0.0d0
a_lhs(1,1)=1.0d0
f_rhs(1)=0.0d0 ! start point for u

a_lhs(2,:)=0.0d0
a_lhs(2,2)=1.0d0
f_rhs(2)=0.0d0 ! start point for v

a_lhs(nn_inflow*2-1,:)=0.0d0
a_lhs(nn_inflow*2-1,nn_inflow*2-1)=1.0d0
f_rhs(nn_inflow*2-1)=0.0d0 ! end point for u

a_lhs(nn_inflow*2,:)=0.0d0
a_lhs(nn_inflow*2,nn_inflow*2)=1.0d0
f_rhs(nn_inflow*2)=0.0d0 ! end point for v

!==================================================

call MIGS(a_lhs,nn_inflow,inv_a,indx)


open(unit=100, file='A', status='unknown')
open(unit=200, file='F', status='unknown')
open(unit=300, file='INVA', status='unknown')
do ia=1,nn_inflow*2
	do ja=1,nn_inflow*2
	write(100,*) a_lhs(ia,ja)
        write(300,*) inv_a(ia,ja)
	end do
        write(200,*) f_rhs(ia)
end do
	
write(*,*) 'After output the linear eqution then stop'

stop

end



