!  calculate 0order derivative of 5th order B Spline fuction!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine B_Spline_0order(dx,hs,S)

  real(8) dx
  real(8) hs
  real(8) S
  integer Btype
  real(8) x
  real(8) hsg

  S=1.0
  Btype=5
  hsg=1.0*hs
   if(Btype==5) then
     x=abs(dx)/hsg
     if ((x.ge.0.0).and.(x.le.1.0))then
        S=((3.0-x)**5.0-6.0*(2.0-x)**5.0+15.0*(1.0-x)**5.0)/120.0
     else if ((x.gt.1.0).and.(x.le.2.0))then
        S=((3.0-x)**5.0-6.0*(2.0-x)**5.0)/120.0
     else if ((x.gt.2.0).and.(x.le.3.0))then
        S=(3.0-x)**5.0/120.0
     else if (x.gt.3.0) then
        S=0.0
     end if
   else if(Btype==3)then
     x=abs(dx)/hsg
     if ((x.ge.0.0).and.(x.le.1.0)) then
	S=2.0/3.0-x**2.0+0.5*x**3.0
     else if ((x.gt.1.0).and.(x.le.2.0))then
	S=1.0/6.0*(2.0-x)**3.0
     else if (x.gt.2.0) then
	S=0.0
     end if
   end if

end subroutine B_Spline_0order
