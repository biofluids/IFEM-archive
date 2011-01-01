subroutine setsolid_id(d,id,mn)
  use solid_variables, only : nn_solid
  implicit none
  
  integer :: mn,id(mn,nn_solid)
  real(8)  d(mn,nn_solid)
  integer :: inc,idf

  do inc=1,nn_solid
     do idf=1,mn 
        if (id(idf,inc).eq.0) d(idf,inc) = 0.0                
     enddo
  enddo
  return
end

