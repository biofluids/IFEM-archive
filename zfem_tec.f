      subroutine  zfem_tec(
     $     klok, td,
     $     ncloud_run,  mncloud_run, mxcloud_run,
     $     dnext_pt,
     $     dlptlocal_number, dlptlocal_head, dlptlocal_tail,
     $     pt_iptexp,
     $     coord_pt, attrib_fcu,
     $     force_con,force_pt,
     $     vel_pt,accel_pt,f1,f2,f3,tem,
     $     d,xn,ien)
      
      implicit real*8 (a-h,o-z)
      include 'r_common'
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

c      dimension  vort( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
c      dimension  f1( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
c      dimension  f2( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
c      dimension  f3( dmnac1:dmxac1, dmnac2:dmxac2, dmnac3:dmxac3)
      
	real* 8 d(ndf,nn)
	real* 8 xn(nsd,nn)
	integer ien(nen,ne)

      klok1= klok/n_step_wr_ib_user_files
	jump_frame=1
      nn_klok=klok1/jump_frame

      if (mod(klok1,jump_frame) .eq. 0) then

c++++++++
c     write fluid domain output in 'tecf.dat'
c++++++++
      write(9000,*) 'Title="fem method tecplot output"'
c         write(9000,*) 'Variables="X","Z","U","W","P","Cir"'
	write(9000,*) 'Variables="X" "Y" "Z" "U" "V" "W" "P"'

      write(9000,777) klok, nn, ne
 777  format(' Zone T="',i3,  
     $       '", N=', i5, ', E=', i5, ', F=FEPOINT, ET=BRICK')
      j=0

c++++++++
c     calculate vorticity
c++++++++
         
c         do 400 k = mnk, mxk
c            do 403 i= dmnlc1+3, dmxlc1
c               deltau = (u(i,j,k+1) - u(i,j,k-1))*unit_velocity
c               deltaz = 2.0d0*unit_length
c               deltaw = (w(i+1,j,k) - w(i-1,j,k))*unit_velocity
c               deltax = 2.0d0*unit_length
c               vort(i,j,k) = deltau/deltaz - deltaw/deltax
c 403        continue
c 400     continue


	do i=1,nn
		write(9000,601) xn(1,i),xn(2,i),xn(3,i),
     +		d(1,i),d(2,i),d(3,i),d(4,i)
	enddo
	do i=1,ne
		write(9000,602) (ien(j,i),j=1,nen)
	enddo
 601	format(7(e14.6,1x))
602	format(8i10)
         
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
            call r_main2(nn_klok,klok)
         endif
      endif
      
      return
      end  

