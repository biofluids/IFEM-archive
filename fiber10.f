      subroutine  ib1_copypointcoordtoconnection(dnext_pt,
     $     dlptlocal_number,dlptlocal_head,dlptlocal_tail,
     $     coord_pt,coord1stpt_con)
      implicit real*8 (a-h,o-z)
c     include 'common' 

      include 'main_common'            
c      include 'ibd0_implementation_parameters.fh'
c      include 'ibd0_application_parameters.fh'
c      include 'ibd0_automatic_parameters.fh'

      integer dnext_pt( mn_point_alloc:mx_point_alloc )
      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail

      real*8 coord_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      real*8 coord1stpt_con( ix:iz,mn_point_alloc:mx_point_alloc )

      integer    icon 
      integer    ipt  

      if (dlptlocal_number  .eq. 0) then
         return
      endif

      do 100 icon = dlptlocal_head, dlptlocal_tail
         ipt = icon
         coord1stpt_con(ix, icon) = coord_pt(ix, ipt)
         coord1stpt_con(iy, icon) = coord_pt(iy, ipt)
         coord1stpt_con(iz, icon) = coord_pt(iz, ipt)
 100  continue
      return
      end

