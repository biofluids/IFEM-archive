      subroutine  iba_calcconnectionforces(dnext_pt,
     $     dlptlocal_number,dlptlocal_head,dlptlocal_tail,
     $     coord1stpt_con,coord2ndpt_con,acttype_con,stiff_con,
     $     rest_con,force_con,nlongcons)
      implicit real*8 (a-h,o-z)
c      include 'common' 
      include 'main_common'            
c      include 'ibd0_implementation_parameters.fh'
c      include 'ibd0_application_parameters.fh'
c      include 'ibd0_automatic_parameters.fh'

      integer dnext_pt(mn_point_alloc:mx_point_alloc )
      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail

      dimension coord1stpt_con(ix:iz,mn_point_alloc:mx_point_alloc )
      dimension coord2ndpt_con(ix:iz,mn_point_alloc:mx_point_alloc )

      dimension stiff_con(mn_point_alloc:mx_point_alloc )
      dimension rest_con(mn_point_alloc:mx_point_alloc )
      dimension acttype_con(mn_point_alloc:mx_point_alloc )
      dimension force_con(ix:iz,mn_point_alloc:mx_point_alloc )

      logical rflag

      dimension unitvector(n_dim_space)
      dimension force(n_dim_space)
      
      nlongcons = 0

      if (dlptlocal_number .eq. 0) then
         return
      endif

c+++++++
c         do 100 icon = dlptlocal_head,dlptlocal_tail - 1
      if (n_ibmfem .eq. 0) then
         ntemp = dlptlocal_head
      else ! fem+ibm
         ntemp = nptfem+1
      endif
         do 100 icon = ntemp,dlptlocal_tail - 1
c++++++
         if (acttype_con(icon) .eq. -1.0) then
            force(1) = 0.0d0
            force(2) = 0.0d0
            force(3) = 0.0d0
         else
            r=sqrt((coord2ndpt_con(ix,icon)-
     $           coord1stpt_con(ix,icon))**2+
     $           (coord2ndpt_con(iy,icon)-
     $           coord1stpt_con(iy,icon))**2
     $           + (coord2ndpt_con(iz,icon) - 
     $           coord1stpt_con(iz,icon))**2)

            ratio = r-rest_con(icon)

            rflag= r .lt. rest_con(icon)
            
            if (rflag) then
               tension = 0.0d0 
            else
               tension = stiff_con(icon)*ratio 
c               write(*,*) .0d0/unit_force*unit_length
c               tension=70.0d0/unit_force*unit_length*ratio
            endif

            if (r .ne. 0.0d0) then
               unitvector(1)=(coord2ndpt_con(ix,icon)-
     $              coord1stpt_con(ix,icon) )/r
               
               unitvector(2)=(coord2ndpt_con(iy,icon)-
     $              coord1stpt_con(iy,icon) )/r
               
               unitvector(3)=(coord2ndpt_con(iz,icon)-
     $              coord1stpt_con(iz,icon) )/r
            endif
            
            force(1) =   tension * unitvector(1)
            force(2) =   tension * unitvector(2)
            force(3) =   tension * unitvector(3)
            
         endif
      
         force_con(ix, icon) = force(1)
         force_con(iy, icon) = force(2)
         force_con(iz, icon) = force(3)

c         write(*,*) force_con(1, 4524)
c         write(*,*) icon, force_con(1, 7006)
c         write(*,*) icon, force_con(1, 2663)
 100  continue
      return
      end






























