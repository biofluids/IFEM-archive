subroutine r_sreadinit(solid_coor_init,solid_coor_curr)
  use r_common
  use solid_variables
  use run_variables, only: tt
  implicit none

  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current
  integer:: i,T0
  real(8) :: Temp

!===========================================================
!Mickael modified this part
 T0=int(tt)
do i=1,nn_solid
                
	if ((tt >= (T0+0.00)).and.(tt <= (T0+0.35))) then


		if (solid_coor_curr(2,i).gt.0.0) then
			Temp=sqrt((1-((solid_coor_curr(1,i)-10.0)**2)/(10.0**2))*0.0025**2)+solid_coor_curr(2,i)
		else
			Temp=-sqrt((1-((solid_coor_curr(1,i)-10.0)**2)/(10.0**2))*0.0025**2)+solid_coor_curr(2,i)
		endif

        	solid_coor_curr(2,i) = Temp *1.0

	elseif ((tt > (T0+0.35)).and.(tt <= (T0+0.90))) then
		 if (solid_coor_curr(2,i).gt.0.0) then
        		Temp=-sqrt((1-((solid_coor_curr(1,i)-10.0)**2)/(10.0**2))*0.0025**2)+solid_coor_curr(2,i)
        	else
        		Temp=sqrt((1-((solid_coor_curr(1,i)-10.0)**2)/(10.0**2))*0.0025**2)+solid_coor_curr(2,i)
        	endif
        
        	solid_coor_curr(2,i) = Temp *1.0

	elseif ((tt > (T0+0.90)).and.(tt <= (T0+1.0))) then
		 if (solid_coor_curr(2,i).gt.0.0) then
                        Temp=sqrt((1-((solid_coor_curr(1,i)-10.0)**2)/(10.0**2))*0.0025**2)+solid_coor_curr(2,i)
                else
                        Temp=-sqrt((1-((solid_coor_curr(1,i)-10.0)**2)/(10.0**2))*0.0025**2)+solid_coor_curr(2,i)
                endif

                solid_coor_curr(2,i) = Temp *1.0
	endif

enddo


  return
end subroutine r_sreadinit
