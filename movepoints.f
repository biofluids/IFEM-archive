      subroutine  movepoints(mnelalloc, mxelalloc,
     $     listnext, list_number, list_head, list_tail,
     $     vel_pt, coord_pt,acttype_con,fix_con)
      
      implicit real*8 (a-h,o-z)
c     include 'common' 
      include 'r_common' 
      include 'main_common'            
c      include 'ibd0_implementation_parameters.fh'
c      include 'ibd0_application_parameters.fh'
c      include 'ibd0_automatic_parameters.fh'

      dimension listnext( mnelalloc:mxelalloc )

      dimension vel_pt( ix:iz, mnelalloc:mxelalloc )
      dimension coord_pt( ix:iz, mnelalloc:mxelalloc )
      dimension acttype_con(mn_point_alloc:mx_point_alloc )
      dimension fix_con(mn_point_alloc:mx_point_alloc )

      if (list_number .eq. 0) then
         return
      endif
      
      ipt = list_head

      do 100 n = 1, list_number
         
         icon = ipt

         if (fix_con(icon) .eq. -1.0 ) then
            
            coord_pt(ix,ipt) = coord_pt(ix,ipt)
            coord_pt(iy,ipt) = coord_pt(iy,ipt)
            coord_pt(iz,ipt) = coord_pt(iz,ipt)
            
         elseif (fix_con(icon) .eq. -2.0) then
            
            coord_pt(ix,ipt) = coord_pt(ix,ipt)
            coord_pt(iy,ipt) = coord_pt(iy,ipt)
            coord_pt(iz,ipt) = coord_pt(iz,ipt)

c            coord_pt(iy,ipt) = coord_pt(iy,ipt)+ vel_pt(iy,ipt) 
c            coord_pt(iz,ipt) = coord_pt(iz,ipt)+ vel_pt(iz,ipt)
c        
c         elseif (fix_con(icon) .eq. 1.0) then
c            
c            coord_pt(ix,ipt) = coord_pt(ix,ipt)
c            coord_pt(iy,ipt) = coord_pt(iy,ipt)
c            coord_pt(iz,ipt) = coord_pt(iz,ipt) + vel_pt(iz,ipt)
   
         else            
            coord_pt(ix,ipt) = coord_pt(ix,ipt) + vel_pt(ix,ipt)
            coord_pt(iy,ipt) = coord_pt(iy,ipt) + vel_pt(iy,ipt)
            coord_pt(iz,ipt) = coord_pt(iz,ipt) + vel_pt(iz,ipt)
            
         endif
         
         ipt = listnext(ipt)
         
 100  continue
      
      return
      end 
      



