      subroutine  zfem_ensCase(time_value, currentStep)
      implicit real*8 (a-h,o-z)
      include 'main_common'

      dimension time_value(10000)
      integer currentStep

      character*12 file_name
      character*13 mgeo_name
      character*12 pre_name
      character*12 vel_name
      character*12 vor_name
      
      character*17 stress_xx
      character*17 stress_yy
      character*17 stress_zz
      character*17 stress_xz
    
      character*17 strain_xx
      character*17 strain_yy
      character*17 strain_zz
      character*17 strain_xz

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
      open (unit=20, file = 'fem.case', status='unknown')
      write(20, '(A6)') 'FORMAT'
      write(20, '(A20)') 'type:   ensight gold'
      write(20, *) 

      write(20, '(A8)') 'GEOMETRY'
      write(file_name,'(A7, A5)')  'fem.geo', '*****'
      write(mgeo_name,'(A8, A5)')  'fem.mgeo', '*****'
      write(pre_name,'(A7, A5)')  'fem.pre', '*****'
      write(vel_name,'(A7, A5)')  'fem.vel', '*****'
      write(vor_name,'(A7, A5)')  'fem.vor', '*****'
      
      write(stress_xx,'(A12, A5)')  'fem.stressXX', '*****'
      write(stress_yy,'(A12, A5)')  'fem.stressYY', '*****'
      write(stress_zz,'(A12, A5)')  'fem.stressZZ', '*****'
      write(stress_xz,'(A12, A5)')  'fem.stressXZ', '*****'
      write(strain_xx,'(A12, A5)')  'fem.strainXX', '*****'
      write(strain_yy,'(A12, A5)')  'fem.strainYY', '*****'
      write(strain_zz,'(A12, A5)')  'fem.strainZZ', '*****'
      write(strain_xz,'(A12, A5)')  'fem.strainXZ', '*****'


      write(20, 5000) 'model:', ts, file_name, 'change_coords_only'
c      write(20, 5001) 'measured:', ts, mgeo_name, 'change_coords_only'
 
      write(20, *) 
      write(20, '(A8)') 'VARIABLE'
      write(20, 5002) 'scalar per node:', ts, 
     $     'pressure', pre_name 
      write(20, 5002) 'scalar per node:', ts, 
     $     'vorticit', vor_name 
      write(20, 5002) 'vector per node:', ts, 
     $     'velocity', vel_name

      write(20, 5003) 'scalar per node:', ts, 
     $     'stress_xx', stress_xx 
      write(20, 5003) 'scalar per node:', ts, 
     $     'stress_yy', stress_yy   
      write(20, 5003) 'scalar per node:', ts, 
     $     'stress_zz', stress_zz   
      write(20, 5003) 'scalar per node:', ts, 
     $     'stress_xz', stress_xz  

      write(20, 5003) 'scalar per node:', ts, 
     $     'strain_xx', strain_xx 
      write(20, 5003) 'scalar per node:', ts, 
     $     'strain_yy', strain_yy 
      write(20, 5003) 'scalar per node:', ts, 
     $     'strain_zz', strain_xx 
      write(20, 5003) 'scalar per node:', ts, 
     $     'strain_xz', strain_xz 

      write(20, *)

c      write(20, 5003) 'tensor symm per node:', ts, 
c     $     'stress', stress_name 
c      write(20, 5003) 'tensor symm per node:', ts, 
c     $     'strain', strain_name 
      write(20, *)

      write(20, '(A4)') 'TIME'
      write(20, 5100) 'time set:', ts
      write(20, 5101) 'number of steps:', numbers_of_step
      write(20, 5102) 'filename start number:', file_start_no
      write(20, 5103) 'filename increment:', file_incre

      write(20, *)

      write(20, '(A12)') 'time values:'
      if(numbers_of_step.le.5) then
         do 1100 i= 1, numbers_of_step
         write(20, '(f14.6)') time_value(i)
 1100    continue
      else
         RESID = mod(numbers_of_step, 5) 
         do 1015 i=1, (numbers_of_step-RESID)/5
            write(20, 5110) (time_value(k+5*(i-1)), k=1,5)
 1015    continue

         do 1020 i=(numbers_of_step-RESID)+1, numbers_of_step
            write(20, '(f14.6)') time_value(i)
 1020    continue
 5110    format(f14.6, 1x,f14.6, 1x,f14.6, 1x,f14.6, 1x,f14.6 )
      endif      


 5000 format(A6, 14x, i5, 5x, A12, 1x, A18)
 5001 format(A9, 11x, i5, 5x, A13, 1x, A18)
 5002 format(A16, 1x, I2, 1x, A8, 1x, A)
 5003 format(A16, 1x, I2, 1x, A9, 1x, A)

 5100 format(A9, 21x, i5)
 5101 format(A16, 14x, i5)
 5102 format(A22, 8x, i5)
 5103 format(A19, 11x, i5)

      close(20)

      return
      end  



      
