   program main

    integer,parameter:: MXND = 100000
!    real(8) crd, vel, acc, frc, dt
!    dimensions crd(3*MXND), vel(3*MXND), acc(3*MXND), frc(3,*MXND)
!   -- I believed this should work too
    real(8) crd(3,MXND), vel(3,MXND), acc(3,MXND), frc(3, MXND)
    real(8) dt
!   Initialization:
    integer opt,i
    logical converged
    
    opt = 1
    mxn = MXND
    dt = 0
    converged = .false.
   i=0
!   -- mars returns  nnd and crd
    call mars3D_RPICFD(opt, dt, mxn, nnd, crd, vel, acc, frc)
    vel(1:3,1) =1.0
    vel(1:3,2)=2.0 
    write(*,*) 'crd=',crd(1:3,1)
    write(*,*) 'crd2=',crd(1:3,2)
!   -- loop
    dt =0.5  ! do not need at this point 
    write(*,*) 'starting data exchange'    
      do while (.not.converged)
	i = i+1 
          if  (i <= 2) then 
          opt = 16
!         -- store velocities in vel and accelerations in acc
          call mars3D_RPICFD(opt, dt, mxn, nnd, crd, vel, acc, frc)
!         -- retrieve forces from frc
 	  write(*,*) 'iteration =',i
	  write(*,*) 'force at node 1=', frc(1:3,1)
	  write(*,*) 'force at node 2=', frc(1:3,2)
	  write(*,*) 'force at node 3=', frc(1:3,3)
       else
          converged = .true.
          opt = 32
!         -- store velocities in vel and accelerations in acc
          call mars3D_RPICFD(opt, dt, mxn, nnd, crd, vel, acc, frc)
	  write(*,*) 'force=', frc(1:3,1), frc(1:3,2)
!         -- retrieve forces from frc
       endif
    enddo

!   -- termination (not necessary)
    stp = 12
    call mars3D_RPICFD(opt, dt, mxn, nnd, crd, vel, acc, frc)

end

