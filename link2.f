      subroutine  setemptylist(imnel,imxel,inullel,
     $     list_number,list_head,list_tail)

      implicit real*8 (a-h,o-z)

      list_number = 0
      list_head  = inullel
      list_tail  = list_head

      return
      end

