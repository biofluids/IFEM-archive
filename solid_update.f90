subroutine solid_update(klok)
  use run_variables, only: tt,dt,its
  use solid_variables
  use r_common
  implicit none

  integer :: klok

  integer :: i_solid,i,j,ni,ie,nos,noj
  integer :: ntem,node
  integer :: k,lx,ly,lz

  real*8 :: todet,wp,viter,tot_vol
  real*8 :: avgvel(3),mom(3),rs(3),xj(3,3),xx(3,9)


! 2. Update material point displacement u(t+dt) = u(t) + dt*v(t)
            
!  write(*,*) 'updating the variables for the structure'
  
  if (nrigid .eq.1) then  !...rigid case: create average velocity and apply to all nodes
     do i_solid=1,n_solid
		mom(1:3)=0.0   ! momentum
		tot_vol=0.0    ! total volume

		do ie=(i_solid-1)*ne_solid_1+1,i_solid*ne_solid_1
           do nos=1,nis
              ntem=solid_fem_con(ie,nos)
              do noj=1,3
                 xx(noj,nos)=solid_coor_init(noj,ntem)
		      enddo
	       enddo
		   do lx=1,nint
	          select case (nis)
	          case (8); rs(1)=xg(lx,nint)
	          case (4); rs(1)=xg_tetra(lx,nint)
	          end select

			  do ly=1,nint
		         select case (nis)
	             case (8); rs(2)=xg(ly,nint)
	             case (4); rs(2)=xg_tetra(ly,nint)
	             end select

			     do lz=1,nint
		            select case (nis)
	                case (8); rs(3)=xg(lz,nint)
	                case (4); rs(3)=xg_tetra(lz,nint)
	                end select

				   !...calculate determinant
				    do i=1,3
				       do j=1,3
					      xj(i,j)=0.0d0
					      do k=1,nis
						     xj(i,j)=xj(i,j)+r_p(j,k)*xx(i,k)
					      enddo
				       enddo
				    enddo
		            todet = xj(1,1) * (xj(2,2)*xj(3,3) - xj(3,2)*xj(2,3))  &
      			          - xj(2,1) * (xj(1,2)*xj(3,3) - xj(3,2)*xj(1,3))  &
      		              + xj(3,1) * (xj(1,2)*xj(2,3) - xj(2,2)*xj(1,3))
	            
				    select case (nis)
	                case (8); wp = wgt(lx,nint)*wgt(ly,nint)*wgt(lz,nint)
	                case (4); wp = wgt_tetra(lx,nint)*wgt_tetra(ly,nint)*wgt_tetra(lz,nint)/6
	                end select

				    tot_vol = tot_vol+wp*todet*density_solid
      			    do ni=1,nis
				        node=solid_fem_con(ie,ni)
			            mom(1:3)=mom(1:3)+wp*todet*density_solid*h(ni)*solid_vel(1:3,node)
				    enddo
			     enddo
			  enddo
		   enddo
		enddo

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

!=======================================================

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