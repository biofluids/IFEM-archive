      subroutine  zibm_tec(
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
      integer dnext_mk( mn_marker_alloc:mx_marker_alloc )
      
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
      dimension  vort_i( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      dimension  vort_j( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      dimension  vort_k( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      
      dimension num_point(max_point)

      character*5  fileroot
      character*12 name_file1
      character*13 name_file2
      character*13 name_file3
      character*9 name_file4
      character*8 name_file5

      ntvs=n_step_wr_ib_user_files

      if (n_step_wr_ib_user_files .eq. 1) then
         njump=1
      else
         njump=2
      endif

      ifileunit = 15
      radp = 0.05d0

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

      write(name_file4,'(A2, A5, A2)')  'up', fileroot, '.m'
      write(name_file5,'(A1, A5, A2)')  's', fileroot, '.m'


c modified on Nov. 5, 2002, change to 3D

c++++++++
c fluid domain
c++++++++

      write(9001,*) 'Title="ibm method, tecplot output"'
      write(9001,*) 'Variables="X","Y","Z","U","V","W","P",
     $     "vor_i","vor_j","vor_k"'

      write(9001,1000) (dmxlc1-dmnlc1)/2+1,  
     $     (dmxlc2-dmnlc2)/2+1, 
     $     (dmxlc3-dmnlc3)/2+1
 1000 format('Zone I=', i5, ', J=',i5, ', K=',i5, ',F=Point')

c++++++++
c     calculate vorticity
c++++++++

      do k = dmnlc3, dmxlc3
       do j = dmnlc2, dmxlc2
         do i= dmnlc1, dmxlc1
            d_u_z = (u(i,j,k+1) - u(i,j,k-1))*unit_velocity
            d_u_y = (u(i,j+1,k+1) - u(i,j-1,k))*unit_velocity
            deltaz = 2.0d0*unit_length
            d_v_x = (v(i+1,j,k) - v(i-1,j,k))*unit_velocity
            d_v_z = (v(i,j,k+1) - v(i,j,k-1))*unit_velocity
            deltax = 2.0d0*unit_length
            d_w_x = (w(i+1,j,k) - w(i-1,j,k))*unit_velocity
            d_w_y = (w(i,j+1,k) - w(i,j-1,k))*unit_velocity
            deltax = 2.0d0*unit_length
            vort_i(i,j,k) = d_w_y/deltay - d_v_z/deltaz
            vort_j(i,j,k) = d_u_z/deltaz - d_w_x/deltax
            vort_k(i,j,k) = d_v_x/deltax - d_u_y/deltay
         enddo
      enddo
      enddo

      do k= dmnlc3,dmxlc3,2
         do j= dmnlc2,dmxlc2,2
            do i= dmnlc1,dmxlc1,2
            write(9001 ,1001) i*unit_length,
     $           j*unit_length,
     $           k*unit_length,
     $           u(i,j,k)*unit_velocity,
     $           v(i,j,k)*unit_velocity,
     $           w(i,j,k)*unit_velocity,
     $           p(i,j,k)*unit_pressure,
     $           vort_i(i,j,k)*unit_pressure,
     $           vort_j(i,j,k)*unit_pressure,
     $           vort_k(i,j,k)*unit_pressure
 1001       format(10(e14.6,1x))
         enddo
      enddo
      enddo
      
c++++++++
c     fiber positions
c++++++++

      ipt = dlptlocal_head

      if (klok .eq. 1) then
         do 131 i=1,num_fiber
            if (num_point(i) .eq. 1) then
               write(9001 ,1117) coord_pt(ix,ipt)*unit_length,
     $              coord_pt(iy,ipt)*unit_length,
     $              coord_pt(iz,ipt)*unit_length,
     $              (klok-1)/ntvs+1
 1117          format('GEOMETRY X=', e14.6, ',Y=', e14.6, ',Z=',e14.6,
     $              'T = CIRCLE',
     $              ',C=BLACK, FC=BLACK, ZN=',i5,',CS=GRID')
               write(9001,*) radp
               ipt=dnext_pt(ipt)
            else
               write(9001 ,1110) (klok-1)/ntvs+1
 1110          format('GEOMETRY X=0, Y=0 Z=0, M=GRID, C=blue,', 
     $              'T=LINE, ZN=',i5,',F=Point')
               write(9001, *) 1
               write(9001 ,*) num_point(i)
               do 810 j=1,num_point(i)
                  write(9001 ,1112) coord_pt(ix,ipt)*unit_length,
     $                 coord_pt(ix,ipt)*unit_length,
     $                 coord_pt(iz,ipt)*unit_length
 1112             format(3(e14.6,1x))
                  ipt=dnext_pt(ipt)
 810           continue
            endif
 131     continue
      else
         do 132 i=1,num_fiber
            if (num_point(i) .eq. 1) then
               write(9001 ,1117) coord_pt(ix,ipt)*unit_length,
     $              coord_pt(iy,ipt)*unit_length,
     $              coord_pt(iz,ipt)*unit_length,
     $              (klok-1)/ntvs+njump
               write(9001,*) radp
               ipt=dnext_pt(ipt)
            else
               write(9001 ,1110) (klok-1)/ntvs+njump
               write(9001, *) 1
               write(9001 ,*) num_point(i)
               do 910 j=1,num_point(i)
                  write(9001 ,1112) coord_pt(ix,ipt)*unit_length,
     $                 coord_pt(iy,ipt)*unit_length,
     $                 coord_pt(iz,ipt)*unit_length
                  ipt=dnext_pt(ipt)
 910           continue
            endif
 132     continue
      endif

      return
      end  
