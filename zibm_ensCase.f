      subroutine  zibm_ensCase(time_value, currentStep,ntsbout)

      include 'main_common'

      dimension time_value(n_step_run)
      integer currentStep

      character*12 file_name
      character*13 mgeo_name
      character*12 pre_name
      character*12 vel_name

      integer ts, fs
      integer file_start_no, file_incre
      integer RESID

      ts = 1
      numbers_of_step = currentStep
      file_start_no = n_step_wr_ib_user_files
      file_incre = n_step_wr_ib_user_files

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     write case file 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      open (unit=20, file = 'ibm.case', status='unknown')
      write(20, '(A6)') 'FORMAT'
	write(20, '(A20)') 'type:   ensight'
c      write(20, '(A20)') 'type:   ensight gold'
      write(20, *) 

      write(20, '(A8)') 'GEOMETRY'
      write(file_name,'(A7, A5)')  'ibm.geo', '*****'
      write(mgeo_name,'(A8, A5)')  'ibm.mgeo', '*****'
      write(pre_name,'(A7, A5)')  'ibm.pre', '*****'
      write(vel_name,'(A7, A5)')  'ibm.vel', '*****'

      write(20, 5000) 'model:', ts, file_name, 'change_coords_only'
      write(20, 5001) 'measured:', ts, mgeo_name, 'change_coords_only'
 
      write(20, *) 
      write(20, '(A8)') 'VARIABLE'
      write(20, 5002) 'scalar per node:', ts, 'pressure', pre_name 
      write(20, 5002) 'vector per node:', ts, 'velocity', vel_name
      write(20, *)

      write(20, '(A4)') 'TIME'
      write(20, 5100) 'time set:', ts
      write(20, 5101) 'number of steps:', numbers_of_step
      write(20, 5102) 'filename start number:', file_start_no
      write(20, 5103) 'filename increment:', file_incre

      write(20, *)

      write(20, '(A12)') 'time values:'
      if(numbers_of_step.le.5) then
         do i= 1, numbers_of_step
            write(20, '(f14.3)') time_value(i)
            write(*, *) 'time_value(i)',time_value(i)
	   enddo
      else
         RESID = mod(numbers_of_step, 5) 
         do i=1, (numbers_of_step-RESID)/5
            write(20, 5110) (time_value(k+5*(i-1)), k=1,5)
	   enddo

         do i=(numbers_of_step-RESID)+1, numbers_of_step
            write(20, '(f14.3)') time_value(i)
	   enddo
 5110    format(f14.3, 1x,f14.3, 1x,f14.3, 1x,f14.3, 1x,f14.3 )
      endif      


 5000 format(A6, 14x, i5, 5x, A12, 1x, A18)
 5001 format(A9, 11x, i5, 5x, A13, 1x, A18)
 5002 format(A16, 1x, I2, 1x, A8, 1x, A)

 5100 format(A9, 21x, i5)
 5101 format(A16, 14x, i5)
 5102 format(A22, 8x, i5)
 5103 format(A19, 11x, i5)

      close(20)

      return
      end  
