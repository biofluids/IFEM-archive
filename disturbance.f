      subroutine  disturbance(u,v,w,istep,mnelalloc,mxelalloc,
     $     listnext, list_number, list_head, list_tail,
     $     dmnac1,dmnac2,dmnac3,dmxac1,dmxac2,dmxac3,
     $     dmnlc1,dmnlc2,dmnlc3,dmxlc1,dmxlc2,dmxlc3,
     $     vel_pt, coord_pt, force_pt, accel_pt)
      
      implicit real*8 (a-h,o-z)
      include 'r_common' 
      include 'main_common'            

      integer    dmnlc1,dmxlc1
      integer    dmnlc2,dmxlc2
      integer    dmnlc3,dmxlc3

      integer    mnelalloc
      integer    mxelalloc

      integer    listnext( mnelalloc:mxelalloc )

      integer    list_number
      integer    list_head
      integer    list_tail  

      integer    dmnac1
      integer    dmnac2
      integer    dmnac3
      integer    dmxac1
      integer    dmxac2
      integer    dmxac3
      real        u 
      real        v 
      real        w 
      dimension   u(dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      dimension   v(dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      dimension   w(dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)

      real       vel_pt( ix:iz, mnelalloc:mxelalloc )
      real       coord_pt( ix:iz, mnelalloc:mxelalloc )
      real       force_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      real       accel_pt( ix:iz, mnelalloc:mxelalloc )

      integer    ipt
      integer    n

      if (list_number .eq. 0) then
         return
      endif
      
c++++++++
c this is to check 3 fiber points case
c++++++++
      if (n_ibmfem .ne. 1) then
c         force_pt(iz,2)= force_pt(iz,2)+0.1*dsin(pi*(istep-1)/2)
c         force_pt(iz,5)= force_pt(iz,5)+0.1*dsin(pi*(istep-1)/2)
c         force_pt(iz,8)= force_pt(iz,8)+0.1*dsin(pi*(istep-1)/2)
c         force_pt(iz,11)= force_pt(iz,11)+0.1*dsin(pi*(istep-1)/2)
      endif

      dmass = 0.1d2
      write(*,*) dmass,accel_pt(iz,3)

      if (n_ibmfem .ne. 1) then
         force_pt(iz,3)= force_pt(iz,3) - dmass*accel_pt(iz,3)
         force_pt(ix,3)= force_pt(ix,3) - dmass*accel_pt(ix,3)
      endif
         
      return
      end 





