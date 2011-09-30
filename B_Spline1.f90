!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!5th order B_Spline function!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine B_Spline1(dx,hs,nsd,S,INFO)

  integer nsd
  real(8) dx(nsd)
  real(8) hs
  real(8) S
  integer Btype
  integer INFO

  integer i
  real(8) x
  real(8) hsg

  S=1.0
  Btype=5
  hsg=1.0*hs
!write(*,*)'hsg=',hsg
  INFO=1
  if(maxval(dx(:)).gt.3.0*hsg) then
    INFO=0
    goto 100
  end if
  if(Btype==5) then
  do i=1,nsd
     x=dx(i)/hsg
     if ((x.ge.0.0).and.(x.le.1.0))then
	S=S*((3.0-x)**5-6.0*(2.0-x)**5+15.0*(1.0-x)**5)/120.0
     else if ((x.gt.1.0).and.(x.le.2.0))then
	S=S*((3.0-x)**5-6.0*(2.0-x)**5)/120.0
     else if ((x.gt.2.0).and.(x.le.3.0))then
	S=S*(3.0-x)**5/120.0
     else if (x.gt.3.0) then
	S=S*0.0
	INFO=0
	goto 100
     end if
  end do
  else if(Btype==3) then
  do i=1,nsd
     x=dx(i)/hsg
     if((x.ge.0.0).and.(x.le.1.0))then
	S=S*(2.0/3.0-(x**2.0)+0.5*(x**3.0))
     else if((x.gt.1.0).and.(x.le.2.0))then
	S=S*(1.0/6.0*((2.0-x)**3.0))
     else if(x.gt.2.0) then
	S=S*0.0
     end if
   end do
  end if
100 continue
end subroutine B_Spline1








