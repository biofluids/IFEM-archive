      subroutine  init_domainaccounting(imnel,imxel,inullel,
     $     dnext_pt,dlfreeel_number,dlfreeel_head,dlfreeel_tail,
     $     dlellocal_number,dlellocal_head,dlellocal_tail,
     $     dlelsend_number,dlelsend_head,dlelsend_tail, 
     $     dlelrecv_number,dlelrecv_head,dlelrecv_tail)
      implicit real*8 (a-h,o-z)

      include 'main_common'      
c      include 'ibd0_implementation_parameters.fh'
c      include 'ibd0_application_parameters.fh'
c      include 'ibd0_automatic_parameters.fh'

      integer dnext_pt( imnel:imxel )

      integer dlfreeel_number
      integer dlfreeel_head
      integer dlfreeel_tail

      integer dlellocal_number
      integer dlellocal_head
      integer dlellocal_tail

      integer dlelsend_number( mn_de:mx_de )
      integer dlelsend_head( mn_de:mx_de )
      integer dlelsend_tail( mn_de:mx_de )

      integer dlelrecv_number( mn_de:mx_de )
      integer dlelrecv_head( mn_de:mx_de )
      integer dlelrecv_tail( mn_de:mx_de )

      call setfulllist(imnel,imxel,inullel,dnext_pt,
     $     dlfreeel_number,dlfreeel_head,dlfreeel_tail)
      call setemptylist(imnel,imxel,inullel,dlellocal_number,
     $     dlellocal_head,dlellocal_tail)

      do 100 idomain = mn_de, mx_de

         call setemptylist(imnel,imxel,inullel,
     $        dlelsend_number(idomain),
     $        dlelsend_head(idomain),dlelsend_tail(idomain)) 
         
         call setemptylist(imnel,imxel,inullel,
     $        dlelrecv_number(idomain),
     $        dlelrecv_head(idomain),dlelrecv_tail(idomain))
 100  continue
      return
      end
