!!calculate 1st order derivative of B_Spline function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine B_Spline_1order(dx,hsg,S)

  real(8) dx
  real(8) hs
  real(8) S
  real(8) x
  integer Btype
  real(8) hsg
  S=1.0
  Btype=5
  hs=1.0*hsg
  x=abs(dx)/hs
if(Btype==5) then
  if ((x.ge.0.0).and.(x.le.1.0)) then
	S=(-5.0/hs*(3.0-x)**4.0+30.0/hs*(2.0-x)**4.0-75.0/hs*(1.0-x)**4.0)/120.0
  else if ((x.gt.1.0).and.(x.le.2.0)) then
	S=(-5.0/hs*(3.0-x)**4.0+30.0/hs*(2.0-x)**4.0)/120.0
  else if ((x.gt.2.0).and.(x.le.3.0)) then
	S=(-5.0/hs*(3.0-x)**4.0)/120.0
  else if (x.gt.3.0) then
	S=0.0
  end if
elseif(Btype==3) then
  if ((x.ge.0.0).and.(x.le.1.0))then
	S=-2.0/hs*x+1.5/hs*x**2.0
  else if ((x.gt.1.0).and.(x.le.2.0))then
	S=-0.5/hs*(2.0-x)**2.0
  else if(x.gt.2.0) then
	S=0.0
  end if
end if

  if((dx*abs(dx)).lt.0.0) then
    S=-S
  end if  ! derivative of absolute value
end subroutine B_Spline_1order
