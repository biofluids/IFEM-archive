      subroutine inittypeactivations(td, bwriteon)
      implicit real*8 (a-h,o-z)
      
      include 'iba_application_parameters.fh'
      include 'iba_application_variables.fh'
      
      logical bwriteon 
      
      t=0.0d0   
      tperiod=4.0d0*128.0d0
      alpha_star=0.625d0    
      tactive=tperiod/8.0d0
      tau=tperiod/64.0d0
      alpha_solid=0.0d0
      beta=0.0d0
      
      if (bwriteon) then
c         write(6,*) 'initializing force activation: general values'
c         write(6,*) 'alpha_star', alpha_star
c         write(6,*) 'tperiod   ', tperiod
c         write(6,*) 'tactive   ', tactive
c         write(6,*) 'tperiod / tactive ', tperiod / tactive
         call flush(6)
      endif
      return
      end     

