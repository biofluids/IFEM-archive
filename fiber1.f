      subroutine  calc_ptfuturedomain(dnext_pt,
     $     dlptlocal_number, dlptlocal_head, 
     $     dlptlocal_tail,coord_pt,idomain_pt) 
      implicit real*8 (a-h,o-z)
c     include 'common' 
      
      include 'main_common'      
c      include 'ibd0_implementation_parameters.fh'
c      include 'ibd0_application_parameters.fh'
c      include 'ibd0_automatic_parameters.fh'

      integer dnext_pt( mn_pt_alloc:mx_pt_alloc )
      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail
             
      dimension coord_pt(ix:iz,mn_pt_alloc:mx_pt_alloc)
             
      dimension idomain_pt(mn_pt_alloc:mx_pt_alloc)

      integer calcperbinxyz

      if (dlptlocal_number .eq. 0) then
         return
      endif
      ipt = dlptlocal_head
      
      do 100 n = 1, dlptlocal_number
         
         idomain_pt(ipt) = mn_de + calcperbinxyz(
     $        coord_pt(ix,ipt),coord_pt(iy,ipt),coord_pt(iz,ipt),
     $        n_lc1,n_lc2,n_lc3,n_de1,n_de2,n_de3)
         
         ipt = dnext_pt(ipt)
         
 100  continue
      return
      end

