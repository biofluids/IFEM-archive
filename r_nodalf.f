      subroutine r_nodalf
      implicit real*8 (a-h,o-z)
      include 'r_common'
      do 10 i=1,numfn
         ni=(ndirfn(i)-1)*nnd+nodefn(i)
         predrf(ni)=predrf(ni)+fnodo(nodefn(i),ndirfn(i))
   10 continue
      return
      end
