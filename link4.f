      subroutine  alloclist(
     $     nelalloc,imnel, imxel,inullel,list_next, 
     $     listfree_number,listfree_head,listfree_tail,
     $     list_number, list_head, list_tail)

      implicit real*8 (a-h,o-z)

      dimension   list_next( imnel:imxel )

      call setemptylist(imnel, imxel, inullel, 
     $     list_number, list_head, list_tail)
      
      if (nelalloc .eq. 0) then
         return
      end if
      
      nelfree = listfree_number
      
      if (nelalloc .gt. nelfree) then
         call exit()  
      end if

      iel = listfree_head

      if (nelalloc .gt. 1) then
        do 100 n = 1, nelalloc-1
          iel= list_next(iel)
 100   continue
      endif

      list_number   = nelalloc
      list_head     = listfree_head
      list_tail     = iel

      listfree_number = listfree_number - list_number
      listfree_head   = list_next(list_tail)

      list_next( list_tail) = inullel

      if (listfree_number .eq. 0) then
         listfree_tail = inullel
      endif

      
      return
      end
