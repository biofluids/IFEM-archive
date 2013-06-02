subroutine getinf_rkpm(inf,ninf,x,xna,adist,nn,nsd,maxconn)
  implicit none

  integer :: ninf,nn,nsd,maxconn
  real(8) x(nsd), xna(nsd,nn), adist(nsd,nn)
  real(8) r(nsd)
  integer inf(maxconn)
  integer i
! comments by Jubiao Yang on 03/14/13
! xna - coord. of all fluid nodes (nsd,nn)
! x -   coord. of a solid node (nsd)


!cccccccccccccccccc 
!   x = the coordinate of the point to be calculated for
!   xna = the coordinate of all points
!   r = the distance between point x and other points in the system
!   inf = a collection of all the influence points
!   ninf = total number of influence points
!   adist = the radial distance of the influence domain
!cccccccccccccccccc
 !  open(unit=400, file='interface_RKPM.dat',status='unknown')
  ninf = 0
  do i = 1,nn
     r(1:nsd) = x(1:nsd) - xna(1:nsd,i)
         if (nsd==3) then
                if ((abs(r(1))<=2*adist(1,i)).and.(abs(r(2))<=2*adist(2,i)).and.(abs(r(3))<=2*adist(3,i))) then
                        ninf = ninf + 1
                        inf(ninf) = i
                endif
         elseif (nsd==2) then
                if ((abs(r(1))<=2*adist(1,i)).and.(abs(r(2))<=2*adist(2,i))) then
                        ninf = ninf + 1
                        inf(ninf) = i
                     !   write(400,*) i  
                     !   write(*,*) '****I am in getinf ****'      
                endif
         endif
   enddo
  !write(*,*) 'ninf=  ', inf(1)
  if (ninf > maxconn) then
     write (*,*) "Too many influence nodes!"
     write (*,*) ninf
  elseif (ninf.lt.4) then
     write (*,*) "Not enough influence nodes!"
     write (*,*) ninf
     write (*,*) x(1), x(2)     ! added by Jubiao Yang on 03/14/2013 to see what solid node cannot find enough inf nodes
  endif

  return
end subroutine getinf_rkpm

