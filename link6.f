      subroutine  splitlistkey(imnel, imxel, inullel, 
     $     list_next,listmain_number,listmain_head,listmain_tail,
     $     key, imnkey, imxkey,listskey_next,   
     $     listskey_number, listskey_head, listskey_tail)

      implicit real*8 (a-h,o-z)     

      dimension    list_next( imnel:imxel )
      dimension    key( imnel:imxel )
      dimension    listskey_number( imnkey:imxkey )
      dimension    listskey_head( imnkey:imxkey )
      dimension    listskey_tail( imnkey:imxkey )
      dimension    listskey_next( imnel-1:imxel )

      do 100 i = imnkey, imxkey
         call setemptylist(imnel, imxel, inullel,
     $        listskey_number(i), 
     $        listskey_head(i), 
     $        listskey_tail(i)) 
 100  continue
      
      if (listmain_number .eq. 0) then
         return                 ! ret !
      endif
      iel = listmain_head

      do 200 n = 1, listmain_number
         
         if (listskey_number(key(iel)) .gt. 0) then
            listskey_next( listskey_tail(key(iel))) = iel
         else
            listskey_head(           key(iel) ) = iel
         endif
         
         listskey_tail( key(iel)) = iel
         listskey_number(key(iel)) = listskey_number(key(iel)) + 1
         
         iel = list_next(iel)
         
 200  continue


      do 300 i = imnkey, imxkey
         listskey_next(listskey_tail(i)) = inullel
 300  continue
      

      return
      end
