      subroutine  calcaccel(klok,mnelalloc,mxelalloc,
     $     listnext, list_number, list_head, list_tail,
     $     vel_pt,prevel_pt,accel_pt,dt)
      
      implicit real*8 (a-h,o-z)
c     include 'common' 
      include 'main_common' 
      
      dimension listnext( mnelalloc:mxelalloc )
      
      dimension vel_pt( ix:iz, mnelalloc:mxelalloc )
      dimension accel_pt( ix:iz, mnelalloc:mxelalloc )
      dimension prevel_pt( ix:iz, mnelalloc:mxelalloc )

      if (list_number .eq. 0) then
         return
      endif
      
      ipt = list_head
      
      do 100 n = 1, list_number

!         accel_pt(ix,ipt) = (vel_pt(ix,ipt)-prevel_pt(ix,ipt))/xfactor     
!         accel_pt(iy,ipt) = (vel_pt(iy,ipt)-prevel_pt(iy,ipt))/xfactor     
!         accel_pt(iz,ipt) = (vel_pt(iz,ipt)-prevel_pt(iz,ipt))/xfactor
         accel_pt(ix,ipt) = (vel_pt(ix,ipt)-prevel_pt(ix,ipt))/dt     
         accel_pt(iy,ipt) = (vel_pt(iy,ipt)-prevel_pt(iy,ipt))/dt     
         accel_pt(iz,ipt) = (vel_pt(iz,ipt)-prevel_pt(iz,ipt))/dt


         prevel_pt(ix,ipt) = vel_pt(ix,ipt)
         prevel_pt(iy,ipt) = vel_pt(iy,ipt)
         prevel_pt(iz,ipt) = vel_pt(iz,ipt)
        
         ipt = listnext(ipt)
         
 100  continue

      return
      end 
      
