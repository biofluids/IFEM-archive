subroutine r_sreadinit(solid_coor_init,solid_coor_curr)
  use r_common
  use solid_variables
  implicit none

	
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_init   !...node position initial
  real(8),dimension(1:nsd_solid,1:nn_solid) :: solid_coor_curr   !...node position current

  integer :: ndum,ndumtest,numBC,i,nn
  real(8) :: x1,y1 !,z1

 ! open(10,file='input_solid_BC_dep.in',status='old',action='read')

 ! read(10,*) numBC
 ! write(*,*) "number of different BC types defined = ",numBC

write(*,*)'nn= ', nn_solid

   write(*,*) solid_coor_init(2,1)

do i=1,nn_solid
 if (solid_coor_init(1,i) .gt. 0.0) then
	solid_coor_curr(1,i) = solid_coor_init(1,i) +0.15 
        solid_coor_curr(2,i) = solid_coor_init(2,i) 

!	write(*,*) solid_coor_curr(2,i)
!	 write(*,*) solid_coor_init(2,i)
elseif (solid_coor_init(1,i) .le. 0.0) then
	solid_coor_curr(1,i) = solid_coor_init(1,i) -0.15
        solid_coor_curr(2,i) = solid_coor_init(2,i) 
endif
enddo













!  do i=1,numBC
!...assign nonlinear initial conditions
! 1000
!  read(10,*) ndum,x1,y1 !,z1
!	write(*,*) ndum, x1, y1
!  ndumtest=ndum-1
!  if (ndumtest >= 0) then
!     if (initdir .eq. 1) then
!        solid_coor_curr(1,ndum) = solid_coor_init(1,ndum) + x1
!        solid_coor_curr(2,ndum) = solid_coor_init(2,ndum) + y1
        !solid_coor_curr(3,ndum) = solid_coor_init(3,ndum) + z1
!	write (*,*) initdir
!	write(*,*) ndumtest
!	write(*,*) solid_coor_curr(1,ndum)


        !xindis(1,ndum)=x1
        !xindis(2,ndum)=y1
        !xindis(3,ndum)=z1
 !    endif
  !   goto 1000
 ! endif
!enddo
!close(10)
  return
end subroutine r_sreadinit
