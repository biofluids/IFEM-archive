      subroutine  activatefcus(dnext_pt,dlptlocal_number,
     $     dlptlocal_head,dlptlocal_tail,acttype_con,
     $     stiff0_con,rest0_con,fix_con,stiff,rest)
      implicit real*8 (a-h,o-z)
c     include 'common' 
      
      include 'main_common'      
c      include 'ibd0_implementation_parameters.fh'
c      include 'ibd0_application_parameters.fh'
c      include 'ibd0_automatic_parameters.fh'

      include 'iba_application_parameters.fh'
      include 'iba_application_variables.fh' 

      integer dnext_pt( mn_pt_alloc:mx_pt_alloc )

      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail

      real*8 acttype_con       
      real*8 stiff0_con        
      real*8 rest0_con         
      real*8 fix_con
      dimension acttype_con( mn_point_alloc:mx_point_alloc )
      dimension stiff0_con( mn_point_alloc:mx_point_alloc )
      dimension rest0_con( mn_point_alloc:mx_point_alloc )
      dimension fix_con(mn_point_alloc:mx_point_alloc )

      real*8 stiff( mn_pt_alloc:mx_pt_alloc )
      real*8 rest( mn_pt_alloc:mx_pt_alloc )

      integer    n
      integer    ipt

      logical   flagtmp

      if (dlptlocal_number .eq. 0) then
         return
      endif
      ipt=dlptlocal_head
      
      do 100 n=1, dlptlocal_number
         
         flagtmp =(((tmin1 .lt. (tperiod*
     $        fix_con(ipt))).and.((tperiod*
     $        fix_con(ipt)) .le. tmax1)) .or.
     $        ((tmin2.lt. (tperiod*fix_con(ipt)))
     $        .and. ((tperiod*fix_con(ipt)) 
     $        .le. tmax2)) .and.(acttype_con(ipt) .eq. i_active))
c         if (flagtmp) then
c            rest(ipt)=alpha_solid*rest0_con(ipt)
c         else
c            rest(ipt)=rest0_con(ipt)
c         endif
         rest(ipt)=rest0_con(ipt)
         stiff(ipt)=stiff0_con(ipt)
         ipt=dnext_pt(ipt)
 100  continue
      return
      end 



