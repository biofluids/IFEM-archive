subroutine facemap
  use fluid_variables, only: mapping
  implicit none
 !  triangle 
  mapping(:,1,1) = (/1, 2, 3, 0, 0, 0/)
  mapping(:,2,1) = (/2, 3, 1, 0, 0, 0/)
  mapping(:,3,1) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,4,1) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,5,1) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,6,1) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,7,1) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,8,1) = (/0, 0, 0, 0, 0, 0/)
 !  quadrilateral
  mapping(:,1,2) = (/1, 2, 3, 4, 0, 0/)
  mapping(:,2,2) = (/2, 3, 4, 1, 0, 0/)
  mapping(:,3,2) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,4,2) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,5,2) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,6,2) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,7,2) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,8,2) = (/0, 0, 0, 0, 0, 0/)
 !  tetrahedron
  mapping(:,1,3) = (/3, 1, 2, 3, 0, 0/)
  mapping(:,2,3) = (/2, 2, 3, 1, 0, 0/)
  mapping(:,3,3) = (/1, 4, 4, 4, 0, 0/)
  mapping(:,4,3) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,5,3) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,6,3) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,7,3) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,8,3) = (/0, 0, 0, 0, 0, 0/)
 !  hexahedron  
  mapping(:,1,4) = (/1, 1, 2, 3, 4, 5/)
  mapping(:,2,4) = (/4, 2, 3, 4, 1, 6/)
  mapping(:,3,4) = (/3, 6, 7, 8, 5, 7/)
  mapping(:,4,4) = (/2, 5, 6, 7, 8, 8/)
  mapping(:,5,4) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,6,4) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,7,4) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,8,4) = (/0, 0, 0, 0, 0, 0/)
 !  space-time triangle 
  mapping(:,1,5) = (/1, 2, 3, 0, 0, 0/)
  mapping(:,2,5) = (/2, 3, 1, 0, 0, 0/)
  mapping(:,3,5) = (/5, 6, 4, 0, 0, 0/)
  mapping(:,4,5) = (/4, 5, 6, 0, 0, 0/)
  mapping(:,5,5) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,6,5) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,7,5) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,8,5) = (/0, 0, 0, 0, 0, 0/)
 !  space-time quadrilateral 
  mapping(:,1,6) = (/1, 2, 3, 4, 0, 0/)
  mapping(:,2,6) = (/2, 3, 4, 1, 0, 0/)
  mapping(:,3,6) = (/6, 7, 8, 5, 0, 0/)
  mapping(:,4,6) = (/5, 6, 7, 8, 0, 0/)
  mapping(:,5,6) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,6,6) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,7,6) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,8,6) = (/0, 0, 0, 0, 0, 0/)
 !  space-time tetrahedron  
  mapping(:,1,7) = (/3, 1, 2, 3, 0, 0/)
  mapping(:,2,7) = (/2, 2, 3, 1, 0, 0/)
  mapping(:,3,7) = (/1, 4, 4, 4, 0, 0/)
  mapping(:,4,7) = (/7, 5, 6, 7, 0, 0/)
  mapping(:,5,7) = (/6, 6, 7, 5, 0, 0/)
  mapping(:,6,7) = (/5, 8, 8, 8, 0, 0/)
  mapping(:,7,7) = (/0, 0, 0, 0, 0, 0/)
  mapping(:,8,7) = (/0, 0, 0, 0, 0, 0/)
 !  space-time hexahedron  
  mapping(:,1,7) = (/1, 1, 2, 3, 4, 5/)
  mapping(:,2,7) = (/4, 2, 3, 4, 1, 6/)
  mapping(:,3,7) = (/3, 6, 7, 8, 5, 7/)
  mapping(:,4,7) = (/2, 5, 6, 7, 8, 8/)
  mapping(:,5,7) = (/9, 9,10,11,12,13/)
  mapping(:,6,7) = (/12,10,11,12, 9,14/)
  mapping(:,7,7) = (/11,14,15,16,13,15/)
  mapping(:,8,7) = (/10,13,14,15,16,16/)
  return
end subroutine facemap
