      subroutine  zibm_tec(
     $     klok, td,
     $     ncloud_run,  mncloud_run, mxcloud_run,
     $	 num_fiber, num_point,
     $     dnext_pt,
     $     dlptlocal_number, dlptlocal_head, dlptlocal_tail,
     $     pt_iptexp,
     $     coord_pt, attrib_fcu,
     $     force_con,force_pt,
     $     vel_pt,accel_pt,f1,f2,f3,tem,
     $     d,xn,ien)

      implicit real*8 (a-h,o-z)
      include 'main_common'
	include 'global.h'

      
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


c      dimension  f1( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
c      dimension  f2( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
c      dimension  f3( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)

	real*8  d(ndf,nn),xn(nsd,nn)
c      dimension  vort_i( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
c      dimension  vort_j( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
c      dimension  vort_k( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      
      dimension num_point(max_point)

      character*5  fileroot
      character*12 name_file1
      character*13 name_file2
      character*13 name_file3
      character*9 name_file4
      character*8 name_file5

	integer ien(nen,ne)

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
      write(9001,*) 'Variables="X","Y","Z","U","V","W","P"'

c      write(9001,*) 'Variables="X","Y","Z","U","V","W","P",
c     $     "vor_i","vor_j","vor_k"'


c     LUCY CHECK THIS HEADER FORMAT!!!!!!
      write(9000,777) klok, nn, ne
 777  format('Zone T=', i3, 
     $       ', N=', i5, ', E=', i5, 'ET=BRICK, F=FEPOINT')

c 1000 format('Zone I=', i5, ', J=',i5, ', K=',i5, ',F=Point')

c++++++++
c     calculate vorticity
c++++++++

c      do k = dmnlc3, dmxlc3
c       do j = dmnlc2, dmxlc2
c         do i= dmnlc1, dmxlc1
c            d_u_z = (u(i,j,k+1) - u(i,j,k-1))*unit_velocity
c            d_u_y = (u(i,j+1,k+1) - u(i,j-1,k))*unit_velocity
c            deltaz = 2.0d0*unit_length
c            d_v_x = (v(i+1,j,k) - v(i-1,j,k))*unit_velocity
c            d_v_z = (v(i,j,k+1) - v(i,j,k-1))*unit_velocity
c            deltax = 2.0d0*unit_length
c            d_w_x = (w(i+1,j,k) - w(i-1,j,k))*unit_velocity
c            d_w_y = (w(i,j+1,k) - w(i,j-1,k))*unit_velocity
c            deltax = 2.0d0*unit_length
c            vort_i(i,j,k) = d_w_y/deltay - d_v_z/deltaz
c            vort_j(i,j,k) = d_u_z/deltaz - d_w_x/deltax
c            vort_k(i,j,k) = d_v_x/deltax - d_u_y/deltay
c         enddo
c      enddo
c      enddo

	do i=1,nn
		write(9000,601) xn(1,nn),xn(2,nn),xn(3,nn),
     +		d(1,nn),d(2,nn),d(3,nn),d(4,nn)
	enddo
 601	format(7(e14.6,1x))
	do i=1,ne
		write(9000,602) (ien(j,i),j=1,nen)
	enddo
602	format(8i10)      
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
