      subroutine  zfem_tec(
     $     klok, td,
     $     ncloud_run,  mncloud_run, mxcloud_run, nmark_run,
     $     nmk_cloud,   mnmk_cloud,  mxmk_cloud,
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
     $     mnlatwr1,  nskipwr1,  mxlatwr1,
     $     mnlatwr2,  nskipwr2,  mxlatwr2,
     $     mnlatwr3,  nskipwr3,  mxlatwr3,
     $     dmnac1,    dmnac2,   dmnac3,
     $     dmxac1,    dmxac2,   dmxac3,
     $     dmnlc1,     dmnlc2,    dmnlc3,
     $     dmxlc1,     dmxlc2,    dmxlc3)
      
      implicit real*8 (a-h,o-z)
      include 'r_common'
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

      
      dimension coord_pt( ix:iz, mn_pt_alloc:mx_pt_alloc )
      dimension attrib_fcu( mn_pt_alloc:mx_pt_alloc,
     $     mn_attr_fcu:mx_attr_fcu )
      
      dimension force_con(ix:iz,mn_point_alloc:mx_point_alloc)
      dimension force_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      dimension vel_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      dimension accel_pt( ix:iz, mn_point_alloc:mx_point_alloc )

      dimension  vort( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      dimension  f1( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      dimension  f2( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      dimension  f3( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      

      dimension  u( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  v( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  w( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )
      dimension  p( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3 )


      klok1= klok/n_step_wr_ib_user_files
      nn=klok1/jump_frame

      if (mod(klok1,jump_frame) .eq. 0) then

c++++++++
c     write fluid domain
c++++++++
         write(9000,*) 'Title="fem method tecplot output"'
         write(9000,*) 'Variables="X","Z","U","W","P","Cir"'
         write(9000,777) (dmxlc1-dmnlc1-3-2*nbou)/2+1, 
     $        (mxk-mnk)/2+1
 777     format('Zone I=', i3, 
     $        ', J=', i3, ', F=Point')
         j=0


c++++++++
c     calculate vorticity
c++++++++
         
         do 400 k = mnk, mxk
            do 403 i= dmnlc1+3, dmxlc1
               deltau = (u(i,j,k+1) - u(i,j,k-1))*unit_velocity
               deltaz = 2.0d0*unit_length
               deltaw = (w(i+1,j,k) - w(i-1,j,k))*unit_velocity
               deltax = 2.0d0*unit_length
               vort(i,j,k) = deltau/deltaz - deltaw/deltax
 403        continue
 400     continue


         do 501 k = mnk, mxk,2
            do 503 i= dmnlc1+3+nbou, dmxlc1-nbou, 2
               write(9000,601) i*unit_length,
     $              k*unit_length,
     $              u(i,j,k)*unit_velocity,
     $              w(i,j,k)*unit_velocity,
     $              p(i,j,k)*unit_pressure,
     $              vort(i,j,k)
 601           format(6(e14.6,1x))
 503        continue
 501     continue
         
c++++++++     
c     write marker point
c++++++++     
c         imk = dlmklocal_head         
c         do 70 nmk = 1, nmk_cloud(1)
c            write(9000,889) nn         
c 889        format('GEOMETRY X=0, Y=0, M=GRID, C=blue,', 
c     $           'T=LINE, ZN=',i4,',F=Point')
c            
c            write(9000,*) 1 
c            write(9000,*) 2
c            write(9000,602) coord_mk(ix,nmk)*unit_length,
c     $           coord_mk(iz,nmk)*unit_length 
c            write(9000,602) coord_mk(ix,nmk)*unit_length*1.01d0,
c     $           coord_mk(iz,nmk)*unit_length*1.01d0
c            
c 602        format(2(e14.6,1x))
c            imk=dnext_mk(imk)
c 70      continue
         
c++++++++     
c     output to tecplot header
c++++++++          
         if (n_ibmfem .eq. 1) then
            write(9001,891) 
 891        format('TITLE="Nonlinear Rubber Program"')
            write(9001,892)
 892        format('VARIABLES=X,Z,DISPX,DISPZ,STRESXX,STRESSZZ,',
     $           'STRESSXZ,STRESSYY,PRE,STRAINXX,STRAINZZ,STRAINXZ,',
     $           'STRAINYY')
            call r_print
            call r_main2(nn)
         endif
      endif
      
      return
      end  

