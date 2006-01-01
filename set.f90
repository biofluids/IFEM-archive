!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setid(d,id,mn)
  use fluid_variables
  implicit none
  
  integer :: mn,id(mn,nn)
  real(8)  d(mn,nn)

  integer :: inc,idf

  do inc=1,nn
     do idf=1,mn 
        if (id(idf,inc).eq.0) d(idf,inc) = 0.0                
     enddo
  enddo

  return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setd(d,f,id,mn)
  use fluid_variables
  implicit none

  integer :: mn,id(mn,nn)
  real(8)  d(mn,nn), f(mn,nn)

  integer :: inc,idf

  do inc=1,nn
     do idf=1,mn
        if (id(idf,inc).eq.0) d(idf,inc) = f(idf,inc)         
     enddo
  enddo

  return
end
