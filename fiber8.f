      subroutine  create_conexchlists(dnext_pt,dlptlocal_number, 
     $     dlptlocal_head,dlptlocal_tail,idomain_ptcon,
     $     dnext_con,dlconshare_number,dlconshare_head,
     $     dlconshare_tail)
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

      dimension idomain_ptcon( mn_pt_alloc:mx_pt_alloc, i_lo:i_hi )

      integer dnext_con( mn_pt_alloc:mx_pt_alloc, i_lo:i_hi )
      integer dlconshare_number( null_de:mx_de, i_lo:i_hi )
      integer dlconshare_head( null_de:mx_de, i_lo:i_hi )
      integer dlconshare_tail( null_de:mx_de, i_lo:i_hi )

      do 100 idircon = i_lo, i_hi
	! in link6.f
        call splitlistkey(mn_pt_alloc,mx_pt_alloc,null_pt,
     $        dnext_pt,dlptlocal_number,dlptlocal_head,
     $        dlptlocal_tail,idomain_ptcon(mn_pt_alloc,idircon), 
     $        null_de,mx_de,dnext_con(mn_pt_alloc,idircon),
     $        dlconshare_number(null_de,idircon),
     $        dlconshare_head(null_de,idircon),
     $        dlconshare_tail(null_de,idircon))
 100  continue
      return
      end
