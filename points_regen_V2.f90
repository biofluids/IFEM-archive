!======================================================
!     Regenerate interfacial points
!======================================================

subroutine points_regen(x,x_inter,x_center,x_inter_regen,nn_inter_regen, &
			I_fluid_center,corr_Ip, hg, infdomain,ien)

  use fluid_variables, only:nsd,nn,ne,nen
  use interface_variables
  use allocate_variables, only:regen_ele_loc,nn_regen_ele_loc,center_domain,nn_center_domain

  real(8) x(nsd,nn), x_inter(nsd,maxmatrix),x_center(nsd,ne)
  real(8) x_inter_regen(nsd,maxmatrix)
  integer nn_inter_regen
  real(8) Ic_inter
  real(8) I_fluid_center(ne),corr_Ip(maxmatrix)
  real(8) hg(ne)
  integer infdomain(maxmatrix),ien(nen,ne)

  integer i,j,isd,inl,node,ie,icount,jcount
  real(8) x_fluid(nsd,nen)   !global coordinates of elements' nodes
  integer nn_sub             !# of candidate points per fluid element edg
  integer nn_ele             !# of candidate points per fluid element
  real(8) x_loc_can(nsd,5000)!local corodinates of candidate points
  real(8) xlocan(nsd), xlocan_temp(nsd)

  real(8) sh(nen) !shape function
  real(8) x_glo(nsd)         !global coordinate for a candidate point
  real(8) distance, hs
  real(8) dx(nsd),Sp
  integer nn_cr
  real(8) x_inter_regen_loc(nsd,5000) !coordiates of local regen points

  real(8) Sp_sum_0d(2), Sp_sum_1d(2,nsd),Sp_sum_2d(2,3*(nsd-1))
  real(8) II,dI(nsd),ddI(3*(nsd-1))

  integer nit
  real(8) err_p,temp,delta(nsd)

  integer d_flag
  real(8) d_range, A
  integer curv_flag
  real(8) pi

  real(8) x_ele(nsd,nen)
  integer inf
  integer nn_local_new
  real(8) dcurv,curv_can
  nn_cr=20
  nn_inter_regen=0
  x_inter_regen(:,:)=0.0

  do ie=1,ne_inter
!  do ie=32,32
     curv_flag=0
     do inl=1,nen
	node=ien(inl,inter_ele(ie))
	x_ele(1:nsd,inl)=x(1:nsd,node)
     end do   !for each element, find the global coordiates for each node

    do inl=1,nen
       x_fluid(1,inl)=mean_coor(1,ie)+1.1*hg(inter_ele(ie))* &
		(mean_coor(3,ie)*cos(pi/4.0+real(inl-1)*pi/2.0)- &
		 mean_coor(4,ie)*sin(pi/4.0+real(inl-1)*pi/2.0))
       x_fluid(2,inl)=mean_coor(2,ie)+1.1*hg(inter_ele(ie))* &
		(mean_coor(3,ie)*sin(pi/4.0+real(inl-1)*pi/2.0)+ &
		 mean_coor(4,ie)*cos(pi/4.0+real(inl-1)*pi/2.0))
    end do
!a    x_fluid(1:nsd,1:nen)=x_ele(1:nsd,1:nen)


     nn_sub=35
     nn_local=0
     do while((nn_local.le.nn_cr).and.(nn_sub.le.35)) !begin regeneration loop
	nn_local=0  !reset nn_local
	x_inter_regen_loc(:,:)=0.0  !reset local coordinates

	if((nsd==2) .and. (nen==4)) then  !2d4n
	  do i=1,nn_sub
	     do j=1,nn_sub
		x_loc_can(1,nn_sub*(i-1)+j)=2.0/nn_sub*i-1.0-1.0/nn_sub
		x_loc_can(2,nn_sub*(i-1)+j)=2.0/nn_sub*j-1.0-1.0/nn_sub
	     end do
	  end do
	  nn_ele=nn_sub**2
	end if       !assign local coordinates for candidate points(2d4n)

	do i=1,nn_ele  !begin candidate points loop
	   if(nsd==2 .and. nen==4) then
             sh(1)=0.25*(1-x_loc_can(1,i))*(1-x_loc_can(2,i))
             sh(2)=0.25*(1+x_loc_can(1,i))*(1-x_loc_can(2,i))
             sh(3)=0.25*(1+x_loc_can(1,i))*(1+x_loc_can(2,i))
             sh(4)=0.25*(1-x_loc_can(1,i))*(1+x_loc_can(2,i))
           end if  !calculate the shape function

	   xlocan(:)=0.0
	   do inl=1,nen
		xlocan(1:nsd)=xlocan(1:nsd)+sh(inl)*x_fluid(1:nsd,inl)
	   end do
!	   d_flag=0
!	   do icount=1,nn_inter
!	      d_range=sqrt((xlocan(1)-x_inter(1,icount))**2+(xlocan(2)-x_inter(2,icount))**2)
!	      if(d_range.lt.hg(inter_ele(ie))/2.0) then
!		d_flag=1
!	      end if
!	   end do
!	   if(d_flag==0) then
!	     goto 200
!	   end if
	
!=========================================================================
!    point projection
	   nit=1
	   err_p=999.0
	   delta(1:nsd)=0.0
	   xlocan_temp(1:nsd)=xlocan(1:nsd)
	 do while((nit.le.10).and.(err_p.gt.1.0e-9))
	   xlocan(1:nsd)=xlocan(1:nsd)+delta(1:nsd)
!	   call get_weight_derivative(xlocan,x_inter,x_center,infdomain,hg,&
!					Sp_sum_0d,Sp_sum_1d,Sp_sum_2d)
	   call get_indicator_derivative(xlocan,x_inter,x_center,infdomain,hg,&
		I_fluid_center,corr_Ip,Sp_sum_0d,Sp_sum_1d,Sp_sum_2d,&
		II,dI,ddI)
!	   if((abs(II-Ic_inter).gt.0.15*Ic_inter).and.(abs(II-Ic_inter).lt.0.01*Ic_inter)) then
	   if(((II-Ic_inter).gt.0.2*Ic_inter).or.((Ic_inter-II).gt.0.2*Ic_inter)) then
	      goto 200
	   end if
	   if(nsd==2) then
	      temp=(dI(1)**2+dI(2)**2)
	      if((temp.lt.0.0001).or.(abs(dI(1)).lt.0.0001))then
		goto 200
	      end if
!	      temp=(dI(1)**2+dI(2)**2)/dI(1)
	      delta(1)=(Ic_inter-II)*dI(1)/temp
	      delta(2)=delta(1)*dI(2)/dI(1)
	   end if
!	   err_p=abs(II-Ic_inter)
	   err_p=max(abs(II-Ic_inter),abs(delta(1)/maxval(hg(:))),abs(delta(2)/maxval(hg(:))))
	   nit=nit+1
	 end do   
!===========================================================================
!write(*,*)err_p,abs(delta(1)/maxval(hg(:))),abs(delta(2)/maxval(hg(:)))
	call get_curv_num(xlocan(1:nsd),x_inter,x_center,hg,infdomain,&
			I_fluid_center,corr_Ip,&
			dcurv,curv_can,dI)
	if(dcurv.gt.10.0) then
	  err_p=999.0
	end if
!	temp=sqrt(dI(1)**2+dI(2)**2)
!	A=ddI(1)/temp+dI(1)*(-0.5)/(temp**3)*(2*dI(1)*ddI(1)+2*dI(2)*ddI(3))+ &
!	  ddI(2)/temp+dI(2)*(-0.5)/(temp**3)*(2*dI(1)*ddI(3)+2*dI(2)*ddI(2))

!	if(abs(A).gt.100.0) then
!	  goto 200
!	end if

        if(err_p.lt.1.0e-9) then
	distance=sqrt((xlocan(1)-xlocan_temp(1))**2+(xlocan(2)-xlocan_temp(2))**2)
	if(distance.le.hg(inter_ele(ie))/3.0) then

	  d_flag=0
	  do icount=1,nn_local
	     temp=(x_inter_regen_loc(1,icount)-xlocan(1))**2+(x_inter_regen_loc(2,icount)-xlocan(2))**2
	     temp=sqrt(temp)
	     if(temp.lt.hg(inter_ele(ie))/real(nn_cr)) then
	       d_flag=1
	       goto 300
	     end if
	  end do

300 continue
	  d_flag=0
	  if(d_flag==0)then 
	    nn_local=nn_local+1
	    x_inter_regen_loc(1:nsd,nn_local)=xlocan(1:nsd)
	  end if
	end if
	end if
200 continue

	end do ! end of candidate point loop
	nn_sub=nn_sub+4
     end do   ! end of do while

     nn_local_new=0
!     do icount=1,nn_local
!        call point_in_ele(x_inter_regen_loc(1:nsd,icount),x_ele,nsd,nen,inf)
!	if(inf==1) then
!	   nn_local_new=nn_local_new+1
!           nn_inter_regen=nn_inter_regen+1
!           x_inter_regen(1:nsd,nn_inter_regen)= &
!			x_inter_regen_loc(1:nsd,icount)
!
!           regen_point_flag(nn_inter_regen)=inter_ele(ie)
!	end if
!     end do
    nn_local_new=0
    do icount=1,nn_local
	inf=0
	do jcount=1,nn_inter_regen
	   temp=(x_inter_regen(1,jcount)-x_inter_regen_loc(1,icount))**2+&
		(x_inter_regen(2,jcount)-x_inter_regen_loc(2,icount))**2
	   temp=sqrt(temp)
	   if(temp.lt.hg(inter_ele(ie))/real(nn_cr)) then
	     inf=1
	     goto 111
	   end if
	end do
	if(inf==0) then
!	call point_in_ele(x_inter_regen_loc(1:nsd,icount),x_ele,nsd,nen,inf)
!	if(inf==1)  then
	    nn_local_new=nn_local_new+1
!	end if
	    nn_inter_regen=nn_inter_regen+1
	    x_inter_regen(1:nsd,nn_inter_regen)= &
		x_inter_regen_loc(1:nsd,icount)
!if(abs(x_inter_regen(1,nn_inter_regen)).lt.1.0e-6) then
!   write(*,*)ie,nn_inter_regen
!end if
	  regen_point_flag(nn_inter_regen)=inter_ele(ie)
	end if
111 continue
    end do




     if(curv_flag==1.or.nn_local==0.or.nn_local_new.lt.3) then
	nn_regen_ele=nn_regen_ele+1
	regen_ele(nn_regen_ele)=inter_ele(ie)
     end if
!if(nn_local.gt.0) then
write(*,*)'ie=',ie,'nn_local',nn_local,'nnloc_new',nn_local_new
!end if

  if(inter_ele(ie)==717) then
  open(336,file='regen717.dat',status='unknown')
  write(336,'(a12)')'#Version 1.0'
  write(336,'(a21)')'#EnSight Point Format'
  do j=1,nn_local_new
     write(336,127)x_inter_regen(1,nn_inter_regen-nn_local_new+j),x_inter_regen(2,nn_inter_regen-nn_local_new+j),0.0
  end do
  close(336)
127 format(f14.10,f14.10,f14.10)
  end if

  end do ! end of loop over interfacial elements

end subroutine points_regen














     
  
