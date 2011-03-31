!!calculate curv using B_Spline function!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_curv_Bspline(xp,xg,Ig_ini,infdomain_inter,corr_Ip,curv_inter,hg,norm_inter)

  use interface_variables
  use fluid_variables,only:nsd,ne,nn
  use mpi_variables
  real(8) xp(nsd,maxmatrix),xg(nsd,ne)
  real(8) Ig_ini(ne)
  integer infdomain_inter(maxmatrix)
  real(8) corr_Ip(maxmatrix)
  real(8) curv_inter(maxmatrix)
  real(8) hg(ne)
  real(8) norm_inter(nsd,maxmatrix)

  integer i,j,isd,jsd
  real(8) hs,Sp(nsd),temp1(nsd),dx,temp2,temp,curv_x,curv_y
  real(8) dI(nsd,nn_inter)    ! 1st order derivative
  real(8) dI2(nsd,nn_inter)   ! 2nd order derivative
  real(8) dxy(nn_inter)       ! dI/dxdy

  integer gcount(maxmatrix),pcount(maxmatrix)
  real(8) length(nsd)
  integer flag,num

  dI(1:nsd,1:nn_inter)=-norm_inter(1:nsd,1:nn_inter) ! give dI/dx,dI/dy,dI/dz at xp
  dI2(:,:) = 0.0
  dxy(:)   = 0.0

  curv_inter(:) = 0.0
  do i=1,nn_inter
!     hs=hg(infdomain_inter(i))
     num=0
     do j=1,ne
	hs=hg(j)
	do isd=1,nsd
	   do jsd=1,nsd
	      dx=xp(jsd,i)-xg(jsd,j)
	      if(jsd==isd) then
		call B_Spline_2order(dx,hs,Sp(jsd))
	      else
		call B_Spline_0order(dx,hs,Sp(jsd))
	      end if
	   end do
	   temp2=1.0
	   do jsd=1,nsd
	      temp2=temp2*Sp(jsd)
	   end do
	   dI2(isd,i)=dI2(isd,i)+Ig_ini(j)*temp2
	end do  ! get dI2/dx2
	do isd=1,nsd
	   dx=xp(isd,i)-xg(isd,j)
	   call B_Spline_1order(dx,hs,temp1(isd))
	end do
	dxy(i)=dxy(i)+Ig_ini(j)*temp1(1)*temp1(2) ! get dI/dxy
!====================================================================	
	length(:)=abs(xp(:,i)-xg(:,j))
	flag=1
	do isd=1,nsd
	   if(length(isd).gt.3*hs) then
	     flag=0
	   end if
	end do
	if(flag==1) then
	   num=num+1
	end if
!=====================================================================
     end do
     gcount(i)=num
  end do
!write(*,*)'curv_inter=',curv_inter(1:nn_inter)
  do i=1,nn_inter
!     hs=hg(infdomain_inter(i))
     num=0
     do j=1,nn_inter
	hs=hg(infdomain_inter(j))/1.0
	do isd=1,nsd
	   do jsd=1,nsd
	      dx=xp(jsd,i)-xp(jsd,j)
	      if(jsd==isd) then
	        call B_Spline_2order(dx,hs,Sp(jsd))
	      else
	        call B_Spline_0order(dx,hs,Sp(jsd))
	      end if
	   end do
	   temp2=1.0
	   do jsd=1,nsd
	      temp2=temp2*Sp(jsd)
	   end do
	   dI2(isd,i)=dI2(isd,i)+corr_Ip(j)*temp2
	end do
	do isd=1,nsd
	   dx=xp(isd,i)-xp(isd,j)
	   call B_Spline_1order(dx,hs,temp1(isd))
	end do
	dxy(i)=dxy(i)+corr_Ip(j)*temp1(1)*temp1(2)
!======================================================
        length(:)=abs(xp(:,i)-xp(:,j))
        flag=1
        do isd=1,nsd
           if(length(isd).gt.3*hs) then
             flag=0
           end if
        end do
        if(flag==1) then
           num=num+1
        end if
!========================================================
     end do
     pcount(i)=num
  end do

  do i=1,nn_inter
     temp=dI(1,i)**2+dI(2,i)**2
     temp=sqrt(temp)
     curv_x=dI2(1,i)/temp+dI(1,i)*(-0.5)/(temp**3)*(2*dI(1,i)*dI2(1,i)+2*dI(2,i)*dxy(i))
     curv_y=dI2(2,i)/temp+dI(2,i)*(-0.5)/(temp**3)*(2*dI(1,i)*dxy(i)+2*dI(2,i)*dI2(2,i))
     curv_inter(i)=curv_x+curv_y
  end do
if(myid==0)then
write(*,*)'curv_max=',maxval(abs(curv_inter(1:nn_inter)))
end if
!open(111,file='bigcurvp.dat',status='unknown')
!  write(111,'(a12)')'#Version 1.0'
!  write(111,'(a21)')'#EnSight Point Format'

!  do i=1,nn_inter
!write(*,*)i,curv_inter(i),'# of g',gcount(i),'# of p',pcount(i)
!  if(abs(curv_inter(i)).gt.10.0) then
!     write(111,999)xp(1,i),xp(2,i),0.0
!  end if
!  end do
!999 format(f14.10,f14.10,f14.10)
!close(111)
end subroutine get_curv_Bspline































