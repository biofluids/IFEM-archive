      subroutine  zmatlab(
     $     klok, td,
     $     ncloud_run,  mncloud_run, mxcloud_run, nmark_run,
     $     nmk_cloud,   mnmk_cloud,  mxmk_cloud,
     $     num_fiber, num_point,
     $     dnext_mk,
     $     dlfreemk_number, dlfreemk_head, dlfreemk_tail,
     $     dlmklocal_number, dlmklocal_head, dlmklocal_tail,
     $     mkpin,
     $     coord_mk,
     $     dnext_pt,
     $     dlptlocal_number, dlptlocal_head, dlptlocal_tail,
     $     pt_iptexp,
     $     coord_pt, attrib_fcu,
     $     force_con,force_pt,
     $     vel_pt,accel_pt,f1,f2,f3,tem,
     $     u,  v,  w, p,
     $     mnlatwr1,  mxlatwr1,
     $     mnlatwr2,  mxlatwr2,
     $     mnlatwr3,  mxlatwr3,
     $     dmnac1,    dmnac2,   dmnac3,
     $     dmxac1,    dmxac2,   dmxac3,
     $     dmnlc1,     dmnlc2,    dmnlc3,
     $     dmxlc1,     dmxlc2,    dmxlc3)

      implicit real*8 (a-h,o-z)
      include 'main_common'

      dimension nmk_cloud(mn_cloud_alloc:mx_cloud_alloc)
      dimension mnmk_cloud(mn_cloud_alloc:mx_cloud_alloc)
      dimension mxmk_cloud(mn_cloud_alloc:mx_cloud_alloc)
      
      integer dnext_mk 
      dimension dnext_mk( mn_marker_alloc:mx_marker_alloc )
      
      integer dlfreemk_number
      integer dlfreemk_head  
      integer dlfreemk_tail  

      integer dlmklocal_number
      integer dlmklocal_head  
      integer dlmklocal_tail  
      
      dimension mkpin( mn_marker_alloc:mx_marker_alloc )
      dimension coord_mk( ix:iz, mn_marker_alloc:mx_marker_alloc)
      
      integer dnext_pt( mn_pt_alloc:mx_pt_alloc )

      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail
      
      integer pt_iptexp( mn_pt_alloc:mx_pt_alloc )
      
      dimension coord_pt( ix:iz, mn_pt_alloc:mx_pt_alloc )
      dimension attrib_fcu( mn_pt_alloc:mx_pt_alloc,
     $     mn_attr_fcu:mx_attr_fcu )
      
      dimension force_con(ix:iz,mn_point_alloc:mx_point_alloc)
      dimension force_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      dimension vel_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      dimension accel_pt( ix:iz, mn_point_alloc:mx_point_alloc )

      integer dmnac1
      integer dmnac2
      integer dmnac3
      integer dmxac1
      integer dmxac2
      integer dmxac3

      integer dmnlc1
      integer dmnlc2
      integer dmnlc3
      integer dmxlc1
      integer dmxlc2
      integer dmxlc3
      
      dimension  f1( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      dimension  f2( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      dimension  f3( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      
      dimension  u( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  v( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  w( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  p( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      

      character*5  fileroot
      character*8 name1
      character*8 name2
      character*8 name3
      character*8 name4
      character*8 name5
      character*8 name6
      character*8 name7
      character*8 name8
      character*8 name9
      character*8 name10
      character*8 name11
      
 200  format(a5   )
 201  format(a4,i1)
 202  format(a3,i2)
 203  format(a2,i3)
 204  format(a1,i4)
 205  format(   i5)
      
      if (klok .eq. 0) then
         write(fileroot, 200) '00000'
      elseif (klok .lt. 10) then
         write(fileroot, 201) '0000',klok
      elseif (klok .lt. 100) then
         write(fileroot, 202) '000' ,klok
      elseif (klok .lt. 1000) then
         write(fileroot, 203) '00'  ,klok
      elseif (klok .lt. 10000) then
         write(fileroot, 204) '0'   ,klok
      elseif (klok .lt. 100000) then    
         write(fileroot, 205)   klok
      else
         write(0,*) 'klok .ge. 100000: modify subroutine createfileroot'
         call exit(1)
      endif
c++++++++
c     matlab output
c++++++++
      write(name1,'(A3,A5)') 'u1_', fileroot
      write(name2,'(A3,A5)') 'u2_', fileroot
      write(name3,'(A3,A5)') 'u3_', fileroot
      write(name4,'(A3,A5)') 'u4_', fileroot
      write(name5,'(A3,A5)') 'u5_', fileroot
      write(name6,'(A3,A5)') 'p1_', fileroot
      write(name7,'(A3,A5)') 'p2_', fileroot
      write(name8,'(A3,A5)') 'p3_', fileroot
      write(name9,'(A3,A5)') 'p4_', fileroot
      write(name10,'(A3,A5)') 'p5_', fileroot
      write(name11,'(A3,A5)') 'pp_', fileroot

      do 30 k = mnk, mxk
         
         write(9500 , 1200) name1, k-mnk+1, u(3,0,k)*unit_velocity
         write(9500 , 1201) name2, k-mnk+1, u(15,0,k)*unit_velocity   
         write(9500 , 1202) name3, k-mnk+1, u(31,0,k)*unit_velocity   
         write(9500 , 1203) name4, k-mnk+1, u(47,0,k)*unit_velocity   
         write(9500 , 1204) name5, k-mnk+1, u(63,0,k)*unit_velocity   
         
         write(9500 , 1300) name6, k-mnk+1, p(3,0,k)*unit_pressure   
         write(9500 , 1301) name7, k-mnk+1, p(15,0,k)*unit_pressure   
         write(9500 , 1302) name8, k-mnk+1, p(31,0,k)*unit_pressure   
         write(9500 , 1303) name9, k-mnk+1, p(47,0,k)*unit_pressure   
         write(9500 , 1304) name10, k-mnk+1, p(63,0,k)*unit_pressure   
         
 30   continue
      
 1200 format(A8, '(',i4,')=',f14.9,';')
 1201 format(A8, '(',i4,')=',f14.9,';')
 1202 format(A8,'(',i4,')=',f14.9,';')
 1203 format(A8,'(',i4,')=',f14.9,';')
 1204 format(A8,'(',i4,')=',f14.9,';')
 1300 format(A8,'(',i4,')=',f14.9,';')
 1301 format(A8,'(',i4,')=',f14.9,';')
 1302 format(A8,'(',i4,')=',f14.9,';')
 1303 format(A8,'(',i4,')=',f14.9,';')
 1304 format(A8,'(',i4,')=',f14.9,';')
      
      
c     horizontal pressure line
c      write(*,*) 'horizontal line'
      do 40 i = dmnlc1+3, dmxlc1,2
         write(9500 , 1440) name11, (i+1)/2-1, p(i,0,32)*unit_pressure   
 40   continue

 1440 format(A8,'(',i4,')=',f14.9,';')
      
      
      return
      end  
