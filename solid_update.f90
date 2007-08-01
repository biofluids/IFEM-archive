subroutine solid_update(klok,solid_fem_con,solid_coor_init,solid_coor_curr,solid_vel,solid_prevel,solid_accel)
  use run_variables, only: tt,dt,its
  use solid_variables
  use r_common
  implicit none

  integer :: klok
  integer,dimension(1:ne_solid,1:nen_solid) :: solid_fem_con     !...connectivity for solid FEM mesh
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_vel         !...velocity
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_prevel      !...velocity - previous timestep
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_accel       !...acceleration
  integer :: i_solid,i,j,k,ni,ie,nos,iq
  integer :: ntem,node
  real(8) :: todet,wp,viter,tot_vol
  real(8) :: avgvel(nsd_solid),mom(nsd_solid),rs(nsd_solid)
  real(8) :: xj(nsd_solid,nsd_solid),xx(3,9)


! 2. Update material point displacement u(t+dt) = u(t) + dt*v(t)
  if (nrigid == 1) then  !...rigid case: create average velocity and apply to all nodes
     do i_solid=1,n_solid
        mom(1:nsd_solid) = 0.0   !...momentum
        tot_vol  = 0.0   !...total volume
        element: do ie=(i_solid-1)*ne_solid_1+1,i_solid*ne_solid_1
           do nos=1,nen_solid
              ntem=solid_fem_con(ie,nos)
              xx(1:nsd_solid,nos)=solid_coor_init(1:nsd_solid,ntem)
           enddo
           gauss_int: do iq = 1,nquad_solid
              rs(1:nsd_solid) = xq_solid(1:nsd_solid,iq)
             !...calculate determinant
              do i=1,nsd_solid
                 do j=1,nsd_solid
                    xj(i,j)=0.0d0
                    do k=1,nen_solid
                       xj(i,j)=xj(i,j)+r_p(j,k)*xx(i,k)
                    enddo
                 enddo
              enddo
			  if(nsd_solid==3) then
              todet = xj(1,1) * (xj(2,2)*xj(3,3) - xj(3,2)*xj(2,3))  &
                    - xj(2,1) * (xj(1,2)*xj(3,3) - xj(3,2)*xj(1,3))  &
                    + xj(3,1) * (xj(1,2)*xj(2,3) - xj(2,2)*xj(1,3))
			  elseif(nsd_solid==2) then
              todet = xj(1,1) * xj(2,2) - xj(1,2)*xj(2,1)
			  endif

              wp = wq_solid(iq)
              tot_vol = tot_vol+wp*todet*density_solid

              do ni=1,nen_solid
                 node=solid_fem_con(ie,ni)
                 mom(1:nsd_solid)=mom(1:nsd_solid)+wp*todet*density_solid*h(ni)*solid_vel(1:nsd_solid,node)
              enddo
           enddo gauss_int
        enddo element
        
        avgvel(1:nsd_solid)=mom(1:nsd_solid)/tot_vol  !calculate average velocity
		do i=1,nsd_solid
			 solid_vel(i,(i_solid-1)*nn_solid_1+1:i_solid*nn_solid_1)=avgvel(i)
		enddo
        
     enddo
  endif

 !...acceleration
  solid_accel(1:nsd_solid,1:nn_solid) = (solid_vel(1:nsd_solid,1:nn_solid) - solid_prevel(1:nsd_solid,1:nn_solid))/dt     

 !...save velocity from previous timestep
  solid_prevel(1:nsd_solid,1:nn_solid) = solid_vel(1:nsd_solid,1:nn_solid)

  viter=0.0d0
  write(*,*) 'maximum solid velocity is =', maxval(solid_vel(1:nsd_solid,:))
  du(1:nsd_solid,1:nn_solid)=solid_vel(1:nsd_solid,1:nn_solid)*dt

  do i=1,nn_solid
     if (klok .eq. 1) then
	    do j=1,nsd_solid
        vnorm = vnorm + du(j,i)**2
		enddo
     else
	    do j=1,nsd_solid
		viter = viter + du(j,i)**2
		enddo	 
     endif
  enddo
           
  if (klok .eq. 1) then
     vnorm=sqrt(vnorm)
     viter=1.0d0
  else
     viter=sqrt(viter)/vnorm
  endif
  write(*,*) ' norm=',viter

 !...update current position

  do j=1,nsd_solid
	solid_coor_curr(j,1:nn_solid) = solid_coor_curr(j,1:nn_solid) + du(j,1:nn_solid)
  enddo

  write(*,*) " solid position updated"


  return
end subroutine solid_update