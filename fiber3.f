      subroutine  calcpointforces(dilocaldomain,dtid,
     $     dmasksharecon,dnext_pt,dlptlocal_number,
     $     dlptlocal_head,dlptlocal_tail,dnext_con,
     $     dlconshare_number,dlconshare_head,
     $     dlconshare_tail,idomain_pt,idomain_ptcon,
     $     coord_pt,coord1stpt_con,acttype_con,stiff0_con,
     $     rest0_con,fix_con,force_fcu,
     $     stiffness_con,restlength_con,force_pt,accel_pt,vel_pt)
      implicit real*8 (a-h,o-z)
      include 'r_common'

      include 'main_common'      
c      include 'ibd0_implementation_parameters.fh'
c      include 'ibd0_application_parameters.fh'
c      include 'ibd0_automatic_parameters.fh'

      integer dilocaldomain
      integer dtid( mn_de:mx_de )

      integer dnext_con( mn_point_alloc:mx_point_alloc, i_lo:i_hi )
      integer dlconshare_number( null_de:mx_de, i_lo:i_hi )
      integer dlconshare_head( null_de:mx_de, i_lo:i_hi )
      integer dlconshare_tail( null_de:mx_de, i_lo:i_hi )

      logical dmasksharecon( mn_de:mx_de )

      integer dnext_pt( mn_point_alloc:mx_point_alloc )

      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail

      real*8 idomain_pt( mn_pt_alloc:mx_pt_alloc )
      real*8 idomain_ptcon( mn_pt_alloc:mx_pt_alloc, i_lo:i_hi )

      real*8 coord_pt( ix:iz, mn_point_alloc:mx_point_alloc )

      real*8 acttype_con        
      real*8 stiff0_con         
      real*8 rest0_con          
      real*8 fix_con 

      dimension acttype_con(mn_point_alloc:mx_point_alloc )
      dimension stiff0_con(mn_point_alloc:mx_point_alloc )
      dimension rest0_con(mn_point_alloc:mx_point_alloc )
      dimension fix_con(mn_point_alloc:mx_point_alloc )

      real*8 coord1stpt_con(ix:iz,mn_point_alloc:mx_point_alloc)
      real*8 stiffness_con(mn_point_alloc:mx_point_alloc)
      real*8 restlength_con(mn_point_alloc:mx_point_alloc)
      real*8 force_fcu(ix:iz,mn_point_alloc:mx_point_alloc,i_lo:i_hi)
      real*8 force_pt(ix:iz, mn_point_alloc:mx_point_alloc)
      real*8 accel_pt( ix:iz, mn_point_alloc:mx_point_alloc )
      dimension vel_pt( ix:iz, mn_point_alloc:mx_point_alloc )

      integer    nlongcons
      integer    n
      integer    icon
      integer    ipt
      integer    ierror

      call activatefcus(dnext_pt,dlptlocal_number,dlptlocal_head,
     $     dlptlocal_tail,acttype_con,stiff0_con,rest0_con,
     $     fix_con,stiffness_con,restlength_con)

      call ib1_copypointcoordtoconnection(dnext_pt,
     $     dlptlocal_number,dlptlocal_head,dlptlocal_tail,
     $     coord_pt,coord1stpt_con)

      
      if (n_ibmfem .eq. 0) then
         call iba_calcconnectionforces(dnext_pt,dlptlocal_number,
     $        dlptlocal_head,dlptlocal_tail,coord1stpt_con,
     $        coord1stpt_con(ix, mn_point_alloc+1),acttype_con,
     $        stiffness_con,restlength_con,
     $        force_fcu(ix,mn_point_alloc, i_hi),nlongcons)
         if (nlongcons .gt. 0) then
            write(6,*) 'the number of long links = ', nlongcons
            call exit(1)
         endif
      endif

c++++
c Oct. 22
         call iba_calcconnectionforces(dnext_pt,dlptlocal_number,
     $        dlptlocal_head,dlptlocal_tail,coord1stpt_con,
     $        coord1stpt_con(ix, mn_point_alloc+1),acttype_con,
     $        stiffness_con,restlength_con,
     $        force_fcu(ix,mn_point_alloc, i_hi),nlongcons)
         if (nlongcons .gt. 0) then
            write(6,*) 'the number of long links = ', nlongcons
            call exit(1)
         endif
c++++
! in fiber11.f
      call ibg_femcalcptforce(dnext_pt,
     $     dlptlocal_number,dlptlocal_head,dlptlocal_tail,
     $     acttype_con,fix_con,coord_pt,
     $     force_fcu(ix,mn_point_alloc,i_hi),force_pt,
     $     accel_pt,vel_pt)

      return
      end 
