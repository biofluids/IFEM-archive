!!calculate 2nd order derivative of B_Spline function
!!!!!!!!!!!!!!!!!!!!!!

subroutine B_Spline_2order(dx,hsg,S)

  real(8) dx
  real(8) hs
  real(8) S
  real(8) x
  integer Btype
  real(8) hsg
  Btype=5
  hs=1.0*hsg
  x=abs(dx)/hs
if(Btype==5) then
  if ((x.ge.0.0).and.(x.le.1.0)) then
	S=(20.0/(hs**2.0)*(3.0-x)**3.0-120.0/(hs**2.0)*(2.0-x)**3.0+300.0/(hs**2.0)*(1.0-x)**3.0)/120.0
  else if ((x.gt.1.0).and.(x.le.2.0)) then
	S=(20.0/(hs**2.0)*(3.0-x)**3.0-120.0/(hs**2.0)*(2.0-x)**3.0)/120.0
  else if ((x.gt.2.0).and.(x.le.3.0)) then
	S=(20.0/(hs**2.0)*(3.0-x)**3.0)/120.0
  else if (x.gt.3.0) then
	S=0.0
  end if
else if(Btype==3)then
   if ((x.ge.0.0).and.(x.le.1.0)) then
	S=-2.0/hs**2.0+3.0/hs**2.0*x
   else if((x.gt.1.0).and.(x.le.2.0))then
	S=1.0/hs**2.0*(2.0-x)
   else if(x.gt.2.0) then
	S=0.0
   end if
end if
!  S=S*hs**2
end subroutine B_Spline_2order
