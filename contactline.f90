
subroutine contactline(x_inter,x,x_center,hg,I_fluid_center,corr_Ip,ien,vel_fluid,vol_nn)

  use fluid_variables, only:nsd,nn,ne,nen
  use interface_variables
  use allocate_variables
  use mpi_variables

  real(8) x_inter(nsd,maxmatrix),x(nsd,nn)
  real(8) x_center(nsd,ne),hg(ne),I_fluid_center(ne),corr_Ip(maxmatrix)
  integer ien(nen,ne)
  real(8) vel_fluid(nsd,nn),vol_nn(nn)

  real(8) pi,theta,min_y,max_y,center_y,bot_y
  real(8) gridy
  real(8) x_fix(nsd),norm_fix(nsd),norm_wall(nsd),x_wall(nsd)
  integer nn_con_regen
  real(8) x_con_regen(nsd,maxmatrix)
  integer nn_inter_outer
  real(8) x_inter_outer(nsd,maxmatrix)

  integer i,j,icount,jcount,inl,nit

  real(8) xlocan(nsd),xlocan_temp(nsd),err_p,delta,xcan(nsd,100)
  real(8) II,dI(nsd),ddI(3*(nsd-1)),norm_p(nsd),curv_p
  real(8) temp,ca !capillary number

  integer loc_can
  real(8) x_loc_can(nsd,50)
  pi=3.14159
  theta=pi/2.0+93*pi/180.0
!=============================
! assign parameters
  gridy=max_hg
  bot_y=0.0
  min_y=0.5*gridy
  max_y=1.5*gridy
  center_y=0.5*(min_y+max_y)

  nn_inter_outer=0
  x_inter_outer(:,:)=0.0
  nn_con_regen=0
  x_con_regen(:,:)=0.0
  do i=1,nn_inter
     if(x_inter(2,i).gt.center_y) then
       nn_inter_outer=nn_inter_outer+1
       x_inter_outer(:,nn_inter_outer)=x_inter(:,i)
     end if
  end do


  loc_can=0
  do i=1,ne_inter
     xlocan(:)=0.0
     do inl=1,nen
        xlocan(:)=xlocan(:)+1.0/real(nen)*x(:,ien(inl,inter_ele(i)))
     end do

     if((xlocan(2).gt.min_y-0.5*gridy).and.(xlocan(2).lt.max_y-0.5*gridy)) then
!       if(myid==0) write(*,*)'find the element'
       xlocan(2)=center_y-gridy/30.0
       
       nit=1
       err_p=999.0
       delta=0.0
       xlocan_temp(:)=xlocan(:)

       do while((nit.lt.5).and.(err_p.gt.1.0e-8))
          xlocan(1)=xlocan(1)+delta
	  call get_indicator_derivative_1st(xlocan,x_inter,x_center,hg,&
	       I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
	  delta=(0.5-II)/dI(1)
	  err_p=abs(0.5-II)
	  nit=nit+1
	  if(abs(dI(1)).lt.1.0e-5) then
	    err_p=999.0
	    nit=999
	  end if
       end do
if(myid==0) then
  if(inter_ele(i)==175) then
     write(*,*)'temp=',temp,'err_p=',err_p,'max_hg/1.5=',max_hg/1.5
  end if
end if

       if(err_p.gt.1.0e-6) goto 200
       temp=(xlocan_temp(1)-xlocan(1))**2+(xlocan_temp(2)-xlocan(2))**2
       temp=sqrt(temp)
       if(temp.gt.max_hg/1.5) goto 200

if(myid==0) then
  if(inter_ele(i)==175) then
     write(*,*)'temp=',temp,'err_p=',err_p,'max_hg/1.5=',max_hg/1.5
  end if
end if


       
       do j=1,loc_can
          temp=(xlocan(1)-x_loc_can(1,j))**2+(xlocan(2)-x_loc_can(2,j))**2
          temp=sqrt(temp)
          if(temp.lt.max_hg/2.0) goto 200
       end do
       loc_can=loc_can+1
       x_loc_can(:,loc_can)=xlocan(:)
        


       call find_ca(xlocan,x,vel_fluid,vol_nn,ca,norm_p,theta) 

       if(norm_p(1).lt.0.0) then
         norm_wall(1)=cos(theta+pi/2.0)
	 norm_wall(2)=sin(theta+pi/2.0)
       else
         norm_wall(1)=cos(pi/2.0-theta)
	 norm_wall(2)=sin(pi/2.0-theta)
       end if
       x_fix(:)=xlocan(:)
       norm_fix(:)=norm_p(:)
       x_wall(2)=bot_y
       call solve_contact(x_fix,x_wall,norm_fix,norm_wall,x_con_regen,nn_con_regen)
      
      end if !! end of regen
200 continue      
  end do ! end of loop over inter elements

  do i=1,nn_inter_outer
     x_inter(:,i)=x_inter_outer(:,i)
  end do

  nn_inter=nn_inter_outer

  do i=1,nn_con_regen
     if(x_con_regen(2,i).gt.bot_y) then
       nn_inter=nn_inter+1
       x_inter(:,nn_inter)=x_con_regen(:,i)
     end if
  end do

end subroutine contactline



















