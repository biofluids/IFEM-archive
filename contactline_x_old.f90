!used to construct contact line

subroutine contactline_x(x_inter,x,norm_inter,conele,ne_conele,x_center,hg,&
			I_fluid_center,corr_Ip,ien)

  use fluid_variables, only:nsd,nn,ne,nen
  use interface_variables
  use allocate_variables
  real(8) x_inter(nsd,maxmatrix),x(nsd,nn),norm_inter(nsd,maxmatrix)
  real(8) x_center(nsd,ne),hg(ne),I_fluid_center(ne),corr_Ip(maxmatrix)
  integer ien(nen,ne)
  real(8) x_inter_temp(nsd,maxmatrix)
  integer i,j,icount,jcount,nn_inter_temp,nn_con
  real(8) x_fix(nsd),norm_fix(nsd),norm_wall(nsd),x_wall(nsd)
  real(8) thelta,pi
  real(8) y_limit
  real(8) yl,yr
  integer nn_con_regen
  real(8) x_con_regen(nsd,maxmatrix)
  integer conele(1000,2),ne_conele(2)

  integer inl,nxcan,nit,flag,flag_n
  real(8) xlocan(nsd),err_p,delta,xcan(nsd,5000)
  real(8) II,dI(nsd),ddI(3*(nsd-1)),norm_p(nsd),curv_p

  real(8) temp, xlocan_temp(nsd)
  real(8) min_y,max_y
  nn_con_regen=0
  pi=3.1415926
  thelta=pi/2.0+45*pi/180.0
  y_limit=-0.5+2.5*max_hg
  nn_inter_temp=0
  nn_con=0

  min_y=minval(x_inter(2,1:nn_inter))
  max_y=maxval(x_inter(2,1:nn_inter))
!  x_wall(2,:)=-0.5
if((min_y.lt.-0.5+2.0*max_hg).and.(max_y.gt.0.5-2.0*max_hg)) then
  do i=1,nn_inter
     if(abs(x_inter(2,i)).lt.abs(y_limit)) then
       nn_inter_temp=nn_inter_temp+1
       x_inter_temp(:,nn_inter_temp)=x_inter(:,i)
     end if
  end do
end if
if((min_y.lt.-0.5+2.0*max_hg).and.(max_y.lt.0.5-2.0*max_hg)) then
  do i=1,nn_inter
     if(x_inter(2,i).gt.y_limit) then
       nn_inter_temp=nn_inter_temp+1
       x_inter_temp(:,nn_inter_temp)=x_inter(:,i)
     end if
  end do
end if

if((min_y.gt.-0.5+2.0*max_hg).and.(max_y.gt.0.5-2.0*max_hg)) then
  do i=1,nn_inter
     if(x_inter(2,i).lt.0.5-2.0*max_hg) then
       nn_inter_temp=nn_inter_temp+1
       x_inter_temp(:,nn_inter_temp)=x_inter(:,i)
     end if
  end do
end if



  nn_con_regen=0
do flag_n=1,2

  if(flag_n==1) then
    if(minval(x_inter(2,1:nn_inter)).gt.-0.5+2.0*max_hg) goto 300
    y_limit=-0.5+2.5*max_hg
  else
    if(maxval(x_inter(2,1:nn_inter)).lt.0.5-2.0*max_hg) goto 300
    y_limit=0.5-2.5*max_hg
  end if
  nxcan=0
  do i=1,ne_conele(flag_n)
     do j=1,ne_inter
        if(conele(i,flag_n)==inter_ele(j)) then
!          nxcan=nxcan+1
	  xlocan(:)=0.0
	  do inl=1,nen
	     xlocan(1)=xlocan(1)+1.0/real(nen)*x(1,ien(inl,conele(i,flag_n)))
	  end do
	  xlocan(2)=y_limit
	  nit=1
	  err_p=999.0
	  delta=0.0
	  xlocan_temp(1:nsd)=xlocan(1:nsd)
	  do while((nit.le.5).and.(err_p.gt.1.0e-8))
	     xlocan(1)=xlocan(1)+delta
	     call get_indicator_derivative_1st(xlocan,x_inter,x_center,hg,&
	          I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
	     delta=(0.5-II)/dI(1)
	     err_p=abs(0.5-II)
	     nit=nit+1
	  end do
	  if(err_p.gt.1.0e-8) goto 200
	  temp=(xlocan_temp(1)-xlocan(1))**2+(xlocan_temp(2)-xlocan(2))**2
	  temp=sqrt(temp)
	  if(temp.gt.max_hg/2.0) goto 200
	  nxcan=nxcan+1
          xcan(:,nxcan)=xlocan(:)
	  flag=1
	  do icount=1,nxcan-1
	     if(abs(xcan(1,nxcan)-xcan(1,icount)).lt.max_hg) then
	       flag=0
	     end if
	  end do
	  if(flag==1) then
	    if(norm_p(1).lt.0) then
	      norm_wall(1)=-sin(thelta)
	    else
	      norm_wall(1)=sin(thelta)
	    end if
	    if(flag_n==1) then
	       norm_wall(2)=-cos(thelta)
	    else
	       norm_wall(2)=cos(thelta)
	    end if
	    x_fix(:)=xlocan(:)
	    norm_fix(:)=norm_p(:)
	    if(flag_n==1) then
              x_wall(2)=-0.5
	    else
	      x_wall(2)=0.5
	    end if
	    call solve_contact(x_fix,x_wall,norm_fix,norm_wall,x_con_regen,nn_con_regen)
	  end if
	end if
200 continue
     end do
  end do
300 continue
end do
  nn_inter=nn_inter_temp
  do i=1,nn_inter_temp
     x_inter(:,i)=x_inter_temp(:,i)
  end do
  do i=1,nn_con_regen
     if(  (abs(x_con_regen(1,i)).lt.1.0) .and. (abs(x_con_regen(2,i)).lt.0.5))then
     nn_inter=nn_inter+1
     x_inter(:,nn_inter)=x_con_regen(:,i)
     end if
  end do

end subroutine contactline_x
