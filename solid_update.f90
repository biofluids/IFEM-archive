subroutine solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel)
  use run_variables, only: tt,dt,its
  use solid_variables
  use r_common
  implicit none

  integer :: klok

  integer,dimension(1:ne_solid,1:nen_solid) :: solid_fem_con   !...connectivity for solid FEM mesh
  real*8,dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  real*8,dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current
  real*8,dimension(1:nsd_solid,1:nn_solid) :: solid_vel         !...velocity
  real*8,dimension(1:nsd_solid,1:nn_solid) :: solid_prevel      !...velocity - previous timestep
  real*8,dimension(1:nsd_solid,1:nn_solid) :: solid_accel       !...acceleration

  integer :: i_solid,i,j,ni,ie,nos,iq
  integer :: ntem,node
  integer :: k

  real*8 :: todet,wp,viter,tot_vol
  real*8 :: avgvel(3),mom(3),rs(3),xj(3,3),xx(3,9)


! Update material point displacement u(t+dt) = u(t) + dt*v(t)
            
!  write(*,*) 'updating the variables for the structure'
  
  if (nrigid == 1) then  !...rigid case: create average velocity and apply to all nodes
     do i_solid=1,n_solid
		mom(1:3) = 0.0   !...momentum
		tot_vol  = 0.0   !...total volume

		element: do ie=(i_solid-1)*ne_solid_1+1,i_solid*ne_solid_1
           do nos=1,nen_solid
              ntem=solid_fem_con(ie,nos)
              xx(1:3,nos)=solid_coor_init(1:3,ntem)
		   enddo
		   gauss_int: do iq = 1,nquad_solid
	          rs(1)=xg(1,iq)
	          rs(2)=xg(2,iq)
	          rs(3)=xg(3,iq)

		     !...calculate determinant
			  do i=1,3
			     do j=1,3
			        xj(i,j)=0.0d0
			        do k=1,nen_solid
			           xj(i,j)=xj(i,j)+r_p(j,k)*xx(i,k)
				    enddo
				 enddo
			  enddo
		      todet = xj(1,1) * (xj(2,2)*xj(3,3) - xj(3,2)*xj(2,3))  &
      		        - xj(2,1) * (xj(1,2)*xj(3,3) - xj(3,2)*xj(1,3))  &
      		        + xj(3,1) * (xj(1,2)*xj(2,3) - xj(2,2)*xj(1,3))
	            
			  wp = wq_solid(iq)

			  tot_vol = tot_vol+wp*todet*density_solid

      		  do ni=1,nen_solid
			     node=solid_fem_con(ie,ni)
			     mom(1:3)=mom(1:3)+wp*todet*density_solid*h(ni)*solid_vel(1:3,node)
			  enddo
		   enddo gauss_int
		enddo element

		avgvel(1:3)=mom(1:3)/tot_vol  !calculate average velocity

		du(1,(i_solid-1)*nn_solid_1+1:i_solid*nn_solid_1)=avgvel(1)
		du(2,(i_solid-1)*nn_solid_1+1:i_solid*nn_solid_1)=avgvel(2)
		du(3,(i_solid-1)*nn_solid_1+1:i_solid*nn_solid_1)=avgvel(3)

		solid_vel(1,(i_solid-1)*nn_solid_1+1:i_solid*nn_solid_1)=avgvel(1)
		solid_vel(2,(i_solid-1)*nn_solid_1+1:i_solid*nn_solid_1)=avgvel(2)
		solid_vel(3,(i_solid-1)*nn_solid_1+1:i_solid*nn_solid_1)=avgvel(3)
     enddo
  endif


 !...acceleration
  solid_accel(1:nsd_solid,1:nn_solid) = (solid_vel(1:nsd_solid,1:nn_solid) - solid_prevel(1:nsd_solid,1:nn_solid))/dt     

 !...save velocity from previous timestep
  solid_prevel(1:nsd_solid,1:nn_solid) = solid_vel(1:nsd_solid,1:nn_solid)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  viter=0.0d0
  do i=1,nn_solid
     du(1,i)=solid_vel(1,i)*dt
     du(2,i)=solid_vel(2,i)*dt
     du(3,i)=solid_vel(3,i)*dt
  enddo

  do i=1,nn_solid
     if (klok .eq. 1) then
        vnorm = vnorm + du(1,i)**2 + du(2,i)**2 + du(3,i)**2
     else
        viter = viter + du(1,i)**2 + du(2,i)**2 + du(3,i)**2
     endif
  enddo
           
  if (klok .eq. 1) then
     vnorm=sqrt(vnorm)
     viter=1.0d0
  else
     viter=sqrt(viter)/vnorm
  endif
  write(*,*) ' norm=',viter

!            do i=1,nnd
!                 predrf(i)=predrf2(i)
!                 predrf(i+nnd)=predrf2(i+nnd)
!			   predrf(i+2*nnd)=predrf2(i+2*nnd)
!                 drf(i)=drf2(i)
!                 drf(i+nnd)=drf2(i+nnd)
!			   drf(i+2*nnd)=drf2(i+2*nnd)
!                 vel_pt(1,i)=0.0d0
!                 vel_pt(2,i)=0.0d0
!			   vel_pt(3,i)=0.0d0
!			enddo


 !...update current position
  solid_coor_curr(1,1:nn_solid) = solid_coor_curr(1,1:nn_solid) + du(1,1:nn_solid)
  solid_coor_curr(2,1:nn_solid) = solid_coor_curr(2,1:nn_solid) + du(2,1:nn_solid)
  solid_coor_curr(3,1:nn_solid) = solid_coor_curr(3,1:nn_solid) + du(3,1:nn_solid)

 !...write to 'vel_time.m' to plot
  write(9500,*) 'vel(',its,')=',solid_vel(1,1),';'
  write(9500,*) 'time(',its,')=',tt,';'

  write(*,*) " solid position updated"

  return
end subroutine solid_update