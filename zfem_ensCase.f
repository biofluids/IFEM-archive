      subroutine  zfem_ensCase(dt, currentStep,ntsbout)
      implicit real*8 (a-h,o-z)
      include 'main_common'

      dimension time_value(10000)
      integer currentStep

      character*12 file_name
      character*13 mgeo_name
      character*12 pre_name
      character*12 vel_name
      character*12 vor_name
	character*15 stress
	character*15 strain     

      integer ts, fs
      integer file_start_no, file_incre
      integer RESID

      ts = 1
      numbers_of_step = currentStep/ntsbout
      file_start_no = 0
      file_incre = ntsbout


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     write case file 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      open (unit=20, file = 'fem.case', status='unknown')
      write(20, '(A6)') 'FORMAT'
      write(20, '(A20)') 'type:   ensight'
      write(20, *) 
      write(20, '(A8)') 'GEOMETRY'

	file_name = 'fem.geo*****'
c	mgeo_name = 'fem.geo*****'
	pre_name  = 'fem.pre*****'
	vel_name  = 'fem.vel*****'
	stress    = 'fem.stress*****'
	strain    = 'fem.strain*****'
 
 5001 format(A9, 11x, i5, 5x, A13, 1x, A18)
 5002 format(A16, 1x, I2, 1x, A8, 1x, A)

 5100 format(A9, 21x, i5)
 5101 format(A16, 14x, i5)
 5102 format(A22, 8x, i5)
 5103 format(A19, 11x, i5)

      write(20, 5000) 'model:', ts, file_name, 'change_coords_only'
 5000 format(A6, 14x, i5, 5x, A12, 1x, A18) 

      write(20, *) 

      write(20, '(A8)') 'VARIABLE'
      write(20, 5002) 'scalar per node:', ts, 'pressure', pre_name 
c      write(20, 5002) 'scalar per node:', ts, 'vorticit', vor_name 
      write(20, 5002) 'vector per node:', ts, 'velocity', vel_name

      write(20, 5003) ts, stress
      write(20, 5004) ts, strain 
 5003 format('tensor symm per node: ',I2, 1x,'stress ', A)
 5004 format('tensor symm per node: ',I2, 1x,'strain ', A)

      write(20, *)
	write(20, *)

      write(20, '(A4)') 'TIME'
      write(20, 5100) 'time set:', ts
      write(20, 5101) 'number of steps:', numbers_of_step+1
      write(20, 5102) 'filename start number:', file_start_no
      write(20, 5103) 'filename increment:', file_incre

      write(20, *)

      write(20, '(A12)') 'time values:'

	j=0
	do i=1,numbers_of_step
		time_value(i)=dt*ntsbout*i
	enddo

	write(20,5110) 0.0,(time_value(k),k=1,numbers_of_step)
 5110 format(5f14.6 )

      close(20)

      return
      end  
