      subroutine  setfulllist(imnel,imxel,inullel, 
     $     list_next,list_number,list_head,list_tail)

      implicit real*8 (a-h,o-z)

      dimension    list_next( imnel:imxel )

      list_number = 1 + (imxel - imnel)
      list_head  = imnel
      list_tail  = imxel

      do 100 i = imnel, imxel - 1
         list_next(i)  = i + 1
 100  continue

      list_next(imxel) = inullel

      return
      end
