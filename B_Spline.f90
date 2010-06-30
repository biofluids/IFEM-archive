!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!5th order B_Spline function!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine B_Spline(dx,hs,nsd,S)

  integer nsd
  real(8) dx(nsd)
  real(8) hs
  real(8) S

  integer i
  real(8) x

  S=1.0
  do i=1,nsd
     x=dx(i)/hs
     if ((x.ge.0.0).and.(x.le.1.0))then
	S=S*((3-x)**5-6*(2-x)**5+15*(1-x)**5)/120
     else if ((x.gt.1.0).and.(x.le.2.0))then
	S=S*((3-x)**5-6*(2-x)**5)/120
     else if ((x.gt.2.0).and.(x.le.3.0))then
	S=S*(3-x)**5/120
     else if (x.gt.3.0) then
	S=S*0.0
     end if
  end do

end subroutine B_Spline








