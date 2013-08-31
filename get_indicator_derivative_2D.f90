!================================================
!used to calculate the derivative of the indicator for a point
!=================================================

subroutine get_indicator_derivative_2D(x,xp,x_center,I_fluid_center,corr_Ip, &
				  II,dI,ddI,&
				  norm_p, curv_p)

  use interface_variables
  use fluid_variables,only:nsd,ne,nn
  use allocate_variables, only:nn_center_domain,center_domain
  use mpi_variables

  real(8) x(nsd),xp(nsd,maxmatrix),x_center(nsd,nn_center)
  real(8) I_fluid_center(nn_center),corr_Ip(maxmatrix)
!  real(8) sum_0d(2),sum_1d(2,nsd),sum_2d(2,3*(nsd-1))
  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) norm_p(nsd),curv_p

  integer i,j,isd,jsd
  real(8) hsg,S(nsd),Sp(nsd),Spp(nsd),temp,dx
  real(8) W0,W1(nsd),W2(3*(nsd-1))!phi; dphi/dx,dphi/dy ;d^2phi/dx^2, d^2phi/dy^2, d^phi/dxdy
  integer flag, jcount
  real(8) curv_A,curv_B
!====================
  real(8) M(nsd+1,nsd+1,nsd*4-2), Mvec(nsd+1,nsd+1,nsd*4-2) !M,Mx,My,Mxx,Myy,Mxy
  real(8) B(nsd+1,nsd*4-2)       !B,Bx,By,Bxx,Byy,Bxy
  real(8) P(nsd+1)
  real(8) vec(nsd+1,nsd*4-2)     !vec,vec_x,vec_y,vec_xx,vec_yy,vec_xy
  integer IP(nsd+1), INFO
  real(8) Phi(nsd*4-2,nn_center_domain)           !Phi,Phix,Phiy,Phixx,Phiyy,Phixy
  real(8) Phi_inter(nsd*4-2,nn_inter)
  integer icount
  real(8) dtemp(nsd*4-2)
  real(8) Mtemp(nsd+1,nsd+1)
!=======================

  integer index_center(nn_center_domain),index_inter(nn_inter)
!used to calculate the 1st&2nd derivative of indicator
  
  II=0.0
  dI(:)=0.0
  ddI(:)=0.0
  index_center(:)=1
  index_inter(:)=1
  do flag=1,2

     if(flag==1) then
	jcount=nn_center_domain
     else if(flag==2) then
	jcount=nn_inter
     end if

  do j=1,jcount
     do isd=1,nsd
	do jsd=1,nsd  !loop over jsd,if jsd=isd, do the derivative

	  if(flag==1) then
	   dx=x(jsd)-x_center(jsd,center_domain(j))
	   if(abs(dx).gt.3.0*hsp) index_center(j)=0
	  else if(flag==2) then
	   dx=x(jsd)-xp(jsd,j)
	   if(abs(dx).gt.3.0*hsp) index_inter(j)=0
	  end if

	   if(jsd==isd) then
		call B_Spline_1order(dx,hsp,Sp(jsd))
		call B_Spline_2order(dx,hsp,Spp(jsd))
	   else
		call B_Spline_0order(dx,hsp,Sp(jsd))
		Spp(jsd)=Sp(jsd)
	   end if
	end do

	W1(isd)=1.0
	W2(isd)=1.0
	do jsd=1,nsd
	   W1(isd)=W1(isd)*Sp(jsd)
	   W2(isd)=W2(isd)*Spp(jsd)
	end do   !calculate S1x, S1xx
     end do

     if(nsd==2) then
	do isd=1,nsd
	   if(flag==1) then
		dx=x(isd)-x_center(isd,center_domain(j))
	   else if(flag==2) then
		dx=x(isd)-xp(isd,j)
	   end if
	   call B_Spline_1order(dx,hsp,Sp(isd))
	   call B_Spline_0order(dx,hsp,Spp(isd))
	end do

	W0=Spp(1)*Spp(2)
	W2(3)=Sp(1)*Sp(2)
     end if

     if(nsd==3) then
       do isd=1,nsd
	   if(flag==1) then
                dx=x(isd)-x_center(isd,center_domain(j))
           else if(flag==2) then
                dx=x(isd)-xp(isd,j)
           end if
           call B_Spline_1order(dx,hsp,Sp(isd))
           call B_Spline_0order(dx,hsp,Spp(isd))
       end do
 
       W0=Spp(1)*Spp(2)*Spp(3)
       W2(4)=Sp(1)*Sp(2)*Spp(3)
       W2(5)=Sp(1)*Spp(2)*Sp(3)
       W2(6)=Spp(1)*Sp(2)*Sp(3)
     end if


     if(flag==1) then  !for center points
        if(nsd==2) then
	Phi(1,j)=W0
	Phi(2:3,j)=W1(1:nsd)
	Phi(4:5,j)=W2(1:nsd)
	Phi(6,j)=W2(3)
	else if(nsd==3) then
        Phi(1,j)=W0
	Phi(2:4,j)=W1(1:nsd)
	Phi(5:7,j)=W2(1:nsd)
	Phi(8:10,j)=W2(4:6)
	end if

!	II=II+W0*I_fluid_center(center_domain(j))
!        dI(1:nsd)=dI(1:nsd)+W1(1:nsd)*I_fluid_center(center_domain(j))
!        ddI(1:nsd)=ddI(1:nsd)+W2(1:nsd)*I_fluid_center(center_domain(j))
!        ddI(3)=ddI(3)+W2(3)*I_fluid_center(center_domain(j))


     else if(flag==2) then ! for inter points

	if(nsd==2) then
	Phi_inter(1,j)=W0
	Phi_inter(2:3,j)=W1(1:nsd)
	Phi_inter(4:5,j)=W2(1:nsd)
	Phi_inter(6,j)=W2(3)
	else if(nsd==3) then
	Phi_inter(1,j)=W0
	Phi_inter(2:4,j)=W1(1:nsd)
	Phi_inter(5:7,j)=W2(1:nsd)
	Phi_inter(8:10,j)=W2(4:6)
	end if

!	II=II+W0*corr_Ip(j)
!	dI(1:nsd)=dI(1:nsd)+W1(1:nsd)*corr_Ip(j)
!	ddI(1:nsd)=ddI(1:nsd)+W2(1:nsd)*corr_Ip(j)
!	ddI(3)=ddI(3)+W2(3)*corr_Ip(j)

     end if
   end do! end of j loop
  end do ! end of flag loop

!==========================================
!used for nonuniform mesh with correction
!==========================================
  M(:,:,:)=0.0
  Mvec(:,:,:)=0.0
  P(:)=0.0
  P(1)=1.0
  B(:,:)=0.0
  vec(:,:)=0.0
  do j=1,nn_center_domain
   if(index_center(j)==1) then
!=====construct vec vecx vecy vecxx vecyy vecxy==========
     vec(1,1)=1.0
     vec(2:nsd+1,1)=x(1:nsd)-x_center(1:nsd,center_domain(j))
     vec(2,2)=1.0
     vec(3,3)=1.0
     hsg=c_w(center_domain(j))
!=========================================================
     do icount=1,nsd+1
	do jcount=1,nsd+1
	   Mvec(icount,jcount,1)=vec(icount,1)*vec(jcount,1)

	   Mvec(icount,jcount,2)=vec(icount,2)*vec(jcount,1)+ &
						       vec(icount,1)*vec(jcount,2)

	   Mvec(icount,jcount,3)=vec(icount,3)*vec(jcount,1)+ &
						       vec(icount,1)*vec(jcount,3)

	   Mvec(icount,jcount,4)=vec(icount,4)*vec(jcount,1)+ &
						       2*vec(icount,2)*vec(jcount,2)+ &
						       vec(icount,1)*vec(jcount,4)

	   Mvec(icount,jcount,5)=vec(icount,5)*vec(jcount,1)+ &
						       2*vec(icount,3)*vec(jcount,3)+ &
						       vec(icount,1)*vec(jcount,5)

	   Mvec(icount,jcount,6)=vec(icount,6)*vec(jcount,1)+ &
						       vec(icount,2)*vec(jcount,3)+ &
						       vec(icount,3)*vec(jcount,2)+ &
						       vec(icount,1)*vec(jcount,6)
	end do
    end do
    do icount=1,nsd+1
	do jcount=1,nsd+1

	   M(icount,jcount,1)=M(icount,jcount,1)+Mvec(icount,jcount,1)*Phi(1,j)/(hsp**nsd)*(hsg**nsd) !M

	   M(icount,jcount,2)=M(icount,jcount,2)+(Mvec(icount,jcount,2)*Phi(1,j)+ &
						  Mvec(icount,jcount,1)*Phi(2,j))/(hsp**nsd)*(hsg**nsd)

	   M(icount,jcount,3)=M(icount,jcount,3)+(Mvec(icount,jcount,3)*Phi(1,j)+ &
						  Mvec(icount,jcount,1)*Phi(3,j))/(hsp**nsd)*(hsg**nsd)

	   M(icount,jcount,4)=M(icount,jcount,4)+(Mvec(icount,jcount,4)*Phi(1,j)+ &
						  2*Mvec(icount,jcount,2)*Phi(2,j)+ &
						  Mvec(icount,jcount,1)*Phi(4,j))/(hsp**nsd)*(hsg**nsd)

	   M(icount,jcount,5)=M(icount,jcount,5)+(Mvec(icount,jcount,5)*Phi(1,j)+ &
						  2*Mvec(icount,jcount,3)*Phi(3,j)+ &
						  Mvec(icount,jcount,1)*Phi(5,j))/(hsp**nsd)*(hsg**nsd)

	   M(icount,jcount,6)=M(icount,jcount,6)+(Mvec(icount,jcount,6)*Phi(1,j)+ &
						  Mvec(icount,jcount,2)*Phi(3,j)+ &
						  Mvec(icount,jcount,3)*Phi(2,j)+ &
						  Mvec(icount,jcount,1)*Phi(6,j))/(hsp**nsd)*(hsg**nsd)

	end do
    end do
   end if  !!!!!!!
  end do
  P(:)=0.0
  P(1)=1.0
  Mtemp(1:nsd+1,1:nsd+1)=M(1:nsd+1,1:nsd+1,1)
if(Mtemp(1,1).lt.1.0e-6) then
  II=999.0
  goto 123
end if
  call DGESV(nsd+1,1,Mtemp,nsd+1,IP,P,nsd+1,INFO)
  B(:,1)=P(:)  ! get B
!if(myid==0)write(*,*)'B=',B(:,1)
!stop
  P(:)=0.0
  do icount=1,nsd+1
     do jcount=1,nsd+1
	P(icount)=P(icount)-M(icount,jcount,2)*B(jcount,1)
     end do
  end do
  Mtemp(1:nsd+1,1:nsd+1)=M(1:nsd+1,1:nsd+1,1)
  call DGESV(nsd+1,1,Mtemp,nsd+1,IP,P,nsd+1,INFO)
  B(:,2)=P(:) ! get Bx

  P(:)=0.0
  do icount=1,nsd+1
     do jcount=1,nsd+1
	P(icount)=P(icount)-M(icount,jcount,3)*B(jcount,1)
     end do
  end do
  Mtemp(1:nsd+1,1:nsd+1)=M(1:nsd+1,1:nsd+1,1)
  call DGESV(nsd+1,1,Mtemp,nsd+1,IP,P,nsd+1,INFO)
  B(:,3)=P(:) ! get By

  P(:)=0.0
  do icount=1,nsd+1
     do jcount=1,nsd+1
	P(icount)=P(icount)-M(icount,jcount,4)*B(jcount,1)- &
			    2*M(icount,jcount,2)*B(jcount,2)
     end do
  end do
  Mtemp(1:nsd+1,1:nsd+1)=M(1:nsd+1,1:nsd+1,1)
  call DGESV(nsd+1,1,Mtemp,nsd+1,IP,P,nsd+1,INFO)
  B(:,4)=P(:) ! get Bxx

  P(:)=0.0
  do icount=1,nsd+1
     do jcount=1,nsd+1
	P(icount)=P(icount)-M(icount,jcount,5)*B(jcount,1)- &
			    2*M(icount,jcount,3)*B(jcount,3)
     end do
  end do
  Mtemp(1:nsd+1,1:nsd+1)=M(1:nsd+1,1:nsd+1,1)
  call DGESV(nsd+1,1,Mtemp,nsd+1,IP,P,nsd+1,INFO)
  B(:,5)=P(:) ! get Byy

  P(:)=0.0
  do icount=1,nsd+1
     do jcount=1,nsd+1
	P(icount)=P(icount)-M(icount,jcount,6)*B(jcount,1)- &
			    M(icount,jcount,2)*B(jcount,3)- &
			    M(icount,jcount,3)*B(jcount,2)
     end do
  end do
  Mtemp(1:nsd+1,1:nsd+1)=M(1:nsd+1,1:nsd+1,1)
  call DGESV(nsd+1,1,Mtemp,nsd+1,IP,P,nsd+1,INFO)
  B(:,6)=P(:) ! get Bxy

  vec(:,:)=0.0
  debugA=0.0
  debugB=0.0
  do j=1,nn_center_domain
   if(index_center(j)==1) then
!====construct vec vecx vecy vecxx vecyy vecxy==========
     vec(1,1)=1.0
     vec(2:nsd+1,1)=x(1:nsd)-x_center(1:nsd,center_domain(j))
     vec(2,2)=1.0
     vec(3,3)=1.0
     dtemp(:)=0.0
     hsg=c_w(center_domain(j))
     do icount=1,nsd+1
	dtemp(1)=dtemp(1)+vec(icount,1)*B(icount,1)
	dtemp(2)=dtemp(2)+vec(icount,2)*B(icount,1)+vec(icount,1)*B(icount,2)
	dtemp(3)=dtemp(3)+vec(icount,3)*B(icount,1)+vec(icount,1)*B(icount,3)
	dtemp(4)=dtemp(4)+vec(icount,4)*B(icount,1)+ &
			  2*vec(icount,2)*B(icount,2)+vec(icount,1)*B(icount,4)
	dtemp(5)=dtemp(5)+vec(icount,5)*B(icount,1)+ &
			  2*vec(icount,3)*B(icount,3)+vec(icount,1)*B(icount,5)
	dtemp(6)=dtemp(6)+vec(icount,6)*B(icount,1)+vec(icount,2)*B(icount,3)+ &
			  vec(icount,3)*B(icount,2)+vec(icount,1)*B(icount,6)
     end do

     W0=dtemp(1)*Phi(1,j)/(hsp**nsd)*(hsg**nsd)
     W1(1)=(dtemp(2)*Phi(1,j)+dtemp(1)*Phi(2,j))/(hsp**nsd)*(hsg**nsd)
     W1(2)=(dtemp(3)*Phi(1,j)+dtemp(1)*Phi(3,j))/(hsp**nsd)*(hsg**nsd)
     W2(1)=(dtemp(4)*Phi(1,j)+2*dtemp(2)*Phi(2,j)+dtemp(1)*Phi(4,j))/(hsp**nsd)*(hsg**nsd)
     W2(2)=(dtemp(5)*Phi(1,j)+2*dtemp(3)*Phi(3,j)+dtemp(1)*Phi(5,j))/(hsp**nsd)*(hsg**nsd)
     W2(3)=(dtemp(6)*Phi(1,j)+dtemp(2)*Phi(3,j)+dtemp(3)*Phi(2,j)+dtemp(1)*Phi(6,j))/(hsp**nsd)*(hsg**nsd)

     II=II+W0*I_fluid_center(center_domain(j))      
     dI(1:nsd)=dI(1:nsd)+W1(1:nsd)*I_fluid_center(center_domain(j))
     ddI(1:nsd)=ddI(1:nsd)+W2(1:nsd)*I_fluid_center(center_domain(j))
     ddI(3)=ddI(3)+W2(3)*I_fluid_center(center_domain(j))
   end if !!!!
  end do
!===========================================================
  do j=1,nn_inter
   if(index_inter(j)==1) then
     vec(1,1)=1.0
     vec(2:nsd+1,1)=x(1:nsd)-xp(1:nsd,j)
     vec(2,2)=1.0
     vec(3,3)=1.0
     dtemp(:)=0.0
     do icount=1,nsd+1
        dtemp(1)=dtemp(1)+vec(icount,1)*B(icount,1)
        dtemp(2)=dtemp(2)+vec(icount,2)*B(icount,1)+vec(icount,1)*B(icount,2)
        dtemp(3)=dtemp(3)+vec(icount,3)*B(icount,1)+vec(icount,1)*B(icount,3)
        dtemp(4)=dtemp(4)+vec(icount,4)*B(icount,1)+ &
                          2*vec(icount,2)*B(icount,2)+vec(icount,1)*B(icount,4)
        dtemp(5)=dtemp(5)+vec(icount,5)*B(icount,1)+ &
                          2*vec(icount,3)*B(icount,3)+vec(icount,1)*B(icount,5)
        dtemp(6)=dtemp(6)+vec(icount,6)*B(icount,1)+vec(icount,2)*B(icount,3)+ &
                          vec(icount,3)*B(icount,2)+vec(icount,1)*B(icount,6)
     end do

     W0=dtemp(1)*Phi_inter(1,j)!/(hsp**nsd)*(hsg**nsd)
     W1(1)=(dtemp(2)*Phi_inter(1,j)+dtemp(1)*Phi_inter(2,j))!/(hsp**nsd)*(hsg**nsd)
     W1(2)=(dtemp(3)*Phi_inter(1,j)+dtemp(1)*Phi_inter(3,j))!/(hsp**nsd)*(hsg**nsd)
     W2(1)=(dtemp(4)*Phi_inter(1,j)+2*dtemp(2)*Phi_inter(2,j)+dtemp(1)*Phi_inter(4,j))!/(hsp**nsd)*(hsg**nsd)
     W2(2)=(dtemp(5)*Phi_inter(1,j)+2*dtemp(3)*Phi_inter(3,j)+dtemp(1)*Phi_inter(5,j))!/(hsp**nsd)*(hsg**nsd)
     W2(3)=(dtemp(6)*Phi_inter(1,j)+dtemp(2)*Phi_inter(3,j)+dtemp(3)*Phi_inter(2,j)+dtemp(1)*Phi_inter(6,j))!/(hsp**nsd)*(hsg**nsd)
       II=II+W0*corr_Ip(j)
       dI(1:nsd)=dI(1:nsd)+W1(1:nsd)*corr_Ip(j)
       ddI(1:nsd)=ddI(1:nsd)+W2(1:nsd)*corr_Ip(j)
       ddI(3)=ddI(3)+W2(3)*corr_Ip(j)
   end if
  end do


  temp=sqrt(dI(1)**2+dI(2)**2)
  norm_p(1:nsd)=-dI(1:nsd)/temp

     curv_A=ddI(1)/temp+dI(1)*(-0.5)/(temp**3)*(2*dI(1)*ddI(1)+2*dI(2)*ddI(3))
     curv_B=ddI(2)/temp+dI(2)*(-0.5)/(temp**3)*(2*dI(1)*ddI(3)+2*dI(2)*ddI(2))
     curv_p=curv_A+curv_B
123 continue
end subroutine get_indicator_derivative_2D
























