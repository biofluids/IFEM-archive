!=================================================

!calculate Bspline function using distance
!compute 00,dx,dy,dz,xx,yy,zz,xy,xz,yz
!00 for no diff, dx=dx,dy=dy,dz=dz,xx=d^2/dx^2,xy=d^2/dxy
!b spline includes 3 terms,fa,fb,fc
!L is the length ((x-xp)^2+(y-yp)^2)^0.5
!derivatives are calcualted using chain theorem
!==================================================

subroutine B_Spline_inter(xloc,xploc,nsd,hs,dtype,S)

  real(8) xloc(nsd),xploc(nsd)
  integer nsd
  real(8) hs,h
  character(2) dtype
  real(8) S
  real(8) L
  real(8) x,y,xp,yp
  real(8) fa,fb,fc,fa_L,fb_L,fc_L,L_x,L_y,L_z   
  h=hs
  S=0.0
  if(nsd==2) then
    x=xloc(1)
    y=xloc(2)
    xp=xploc(1)
    yp=xploc(2)
    L=sqrt((x-xp)**2.0+(y-yp)**2.0)
    fa=(3.0-L/h)**5.0
    fb=6.0*(2.0-L/h)**5.0
    fc=15.0*(1.0-L/h)**5.0
    fa_L=-(5.0*(L/h - 3.0)**4.0)/h
    fb_L=-(30.0*(L/h - 2.0)**4.0)/h
    fc_L=-(75.0*(L/h - 1)**4.0)/h
    L_x=(2.0*x - 2.0*xp)/(2.0*((x - xp)**2.0 + (y - yp)**2.0)**(0.5))
    L_y=(2.0*y - 2.0*yp)/(2.0*((x - xp)**2.0 + (y - yp)**2.0)**(0.5))
    fa_LL=-(20.0*(L/h - 3.0)**3.0)/h**2.0
    fb_LL=-(120.0*(L/h - 2.0)**3.0)/h**2.0
    fc_LL=-(300.0*(L/h - 1.0)**3.0)/h**2.0
    L_xx=1.0/((x - xp)**2.0 + (y - yp)**2.0)**(0.5) - (2.0*x - 2.0*xp)**2.0/(4.0*((x - xp)**2.0 + (y - yp)**2.0)**(1.5))
    L_yy=1.0/((x - xp)**2.0 + (y - yp)**2.0)**(0.5) - (2.0*y - 2.0*yp)**2.0/(4.0*((x - xp)**2.0 + (y - yp)**2.0)**(1.5))
    L_xy=-((2.0*x - 2.0*xp)*(2.0*y - 2.0*yp))/(4.0*((x - xp)**2.0 + (y - yp)**2.0)**(1.5))

!...no derivative
    if(dtype=='00') then
      if((L.ge.0.0) .and. (L.lt.hs)) then
	S=1.0/120.0*(fa-fb+fc)
      else if((L.ge.hs).and.(L.lt.2.0*hs)) then
	S=1.0/120.0*(fa-fb)
      else if((L.ge.2.0*hs) .and. (L.lt.3.0*hs)) then
	S=1.0/120.0*fa
      else if(L.ge.3.0*hs) then
	S=0.0
      end if
!...dS/dx
    else if(dtype=='dx') then
      if((L.ge.0.0) .and. (L.lt.hs)) then
        S=1/120*(fa_L*L_x-fb_L*L_x+fc_L*L_x)
      else if((L.ge.hs).and.(L.lt.2*hs)) then
        S=1/120*(fa_L*L_x-fb_L*Lx)
      else if((L.ge.2*hs) .and. (L.lt.3*hs)) then
        S=1/120*(fa_L*L_x)
      else if(L.ge.3*hs) then
        S=0.0
      end if
!...dS/dy
    else if(dtype=='dy') then
      if((L.ge.0.0) .and. (L.lt.hs)) then
        S=1/120*(fa_L*L_y-fb_L*L_y+fc_L*L_y)
      else if((L.ge.hs).and.(L.lt.2*hs)) then
        S=1/120*(fa_L*L_y-fb_L*L_y)
      else if((L.ge.2*hs) .and. (L.lt.3*hs)) then
        S=1/120*(fa_L*L_y)
      else if(L.ge.3*hs) then
        S=0.0
      end if
!...d^2 S/dx^2
    else if(dtype=='xx') then
      if((L.ge.0.0) .and. (L.lt.hs)) then
        S=1/120*((fa_LL+fb_LL+fc_LL)*L_x**2+(fa_L+fb_L+fc_L)*L_xx)
      else if((L.ge.hs).and.(L.lt.2*hs)) then
        S=1/120*((fa_LL+fb_LL)*L_x**2+(fa_L+fb_L)*L_xx)
      else if((L.ge.2*hs) .and. (L.lt.3*hs)) then
        S=1/120*(fa_LL*L_x**2+fa_L*L_xx)
      else if(L.ge.3*hs) then
        S=0.0
      end if
!...d^2 S/dy^2
    else if(dtype=='yy') then
      if((L.ge.0.0) .and. (L.lt.hs)) then
        S=1/120*((fa_LL+fb_LL+fc_LL)*L_y**2+(fa_L+fb_L+fc_L)*L_yy)
      else if((L.ge.hs).and.(L.lt.2*hs)) then
        S=1/120*((fa_LL+fb_LL)*L_y**2+(fa_L+fb_L)*L_y)
      else if((L.ge.2*hs) .and. (L.lt.3*hs)) then
        S=1/120*(fa_LL*L_y**2+fa_L*L_yy)
      else if(L.ge.3*hs) then
        S=0.0
      end if
!...d^2 S/dxy
    else if(dtype=='xy') then
      if((L.ge.0.0) .and. (L.lt.hs)) then
        S=1/120*((fa_LL+fb_LL+fc_LL)*L_x*L_y+(fa_L+fb_L+fc_L)*L_xy)
      else if((L.ge.hs).and.(L.lt.2*hs)) then
        S=1/120*((fa_LL+fb_LL)*L_x*L_y+(fa_L+fb_L)*L_xy)
      else if((L.ge.2*hs) .and. (L.lt.3*hs)) then
        S=1/120*(fa_LL*L_x*L_y+fa_L*L_xy)
      else if(L.ge.3*hs) then
        S=0.0
      end if
    end if

  end if

end subroutine B_Spline_inter










      
