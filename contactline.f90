
subroutine contactline(x_inter,x,x_center,hg,I_fluid_center,corr_Ip,ien,vel_fluid,vol_nn,&
			infdomain_inter,nn_con_ele,con_ele,norm_con_ele,flag_contact,index_bcnode,I_fluid,ne_regen,ele_regen)

  use fluid_variables, only:nsd,nn,ne,nen
  use interface_variables
  use allocate_variables
  use mpi_variables
  use run_variables,only:its
  include 'mpif.h'

  real(8) x_inter(nsd,maxmatrix),x(nsd,nn)
  real(8) x_center(nsd,ne),hg(ne),I_fluid_center(ne),corr_Ip(maxmatrix),I_fluid(nn)
  integer ien(nen,ne)
  real(8) vel_fluid(nsd,nn),vol_nn(nn)
  integer infdomain_inter(maxmatrix)
  integer nn_con_ele,con_ele(nn_con_ele),index_bcnode(2,nn_con_ele)
  real(8) norm_con_ele(nsd,nn_con_ele)

  real(8) pi,theta,min_y,max_y,center_y,bot_y
  real(8) gridy
  real(8) x_fix(nsd),norm_fix(nsd),norm_wall(nsd),x_wall(nsd)
  integer nn_con_regen,nn_org_regen,nn_bottom ! number of points for contactline/ for original line formed by regeneration
  real(8) x_con_regen(nsd,maxmatrix),x_org_regen(nsd,60),x_bottom(nsd,60),x_con_regen_temp(nsd,maxmatrix)
  integer nn_inter_outer
  real(8) x_inter_outer(nsd,maxmatrix)

  integer i,j,icount,jcount,inl,nit,ie

  real(8) xlocan(nsd),xlocan_temp(nsd),err_p,delta(nsd),xcan(nsd,100)
  real(8) II,dI(nsd),ddI(3*(nsd-1)),norm_p(nsd),curv_p
  real(8) temp,ca !capillary number

  integer loc_can
  real(8) x_loc_can(nsd,5000)
  integer flag,index_con_ele(nn_con_ele)
  real(8) tang(nsd) ! tangential direction
  real(8) norm_e(nsd) ! normal direction for element
  integer temp_n
  real(8) rot_x(nsd) !temp x for rotation
  integer flag_contact,flag_ad_re,flag_contact_temp
  real(8) ratio,anglet
  integer ne_regen,ele_regen(100),ele_regen_temp(100)


  integer nn_loc,base,top,loc_index
  integer lower, local_nn(ncpus),local_nn_temp(ncpus)
  integer flag_ca



  pi=3.14159
  theta=static_angle*pi/180.0
!=============================
! assign parameters
if(myid==0)write(*,*)'begin contact line'

  if(nn_con_ele.le.ncpus) then
    if(myid+1.le.nn_con_ele) then
      nn_loc=1
    else
      nn_loc=0
    end if
  else
     base=floor(real(nn_con_ele)/real(ncpus))
     top=nn_con_ele-base*ncpus
     if(myid+1.le.top) then
        nn_loc=base+1
     else
        nn_loc=base
     end if
  end if




  gridy=0.9*max_hg

  nn_inter_outer=0
  x_inter_outer(:,:)=0.0
  nn_con_regen=0
  x_con_regen(:,:)=0.0
  index_con_ele(:)=0
  flag_contact=0
  flag_contact_temp=0
  ne_regen=0
  loc_can=0
  ne_regen=0
  ele_regen(:)=0
  x_org_regen(:,:)=0.0
  x_bottom(:,:)=0.0
  nn_bottom=0
  nn_org_regen=0
!  do ie=749,749!1,nn_con_ele!71,71!nn_con_ele
  do loc_index=1,nn_loc
     ie=myid+1+(loc_index-1)*ncpus
     if(  (I_fluid(index_bcnode(1,ie))-0.5)*(I_fluid(index_bcnode(2,ie))-0.5).gt.0.0  ) goto 200 ! only perform within the ele with inter points
     i=con_ele(ie)   !this is the element id
      flag_contact_temp=1
     ratio=(I_fluid(index_bcnode(1,ie))-0.5)/(0.5-I_fluid(index_bcnode(2,ie)))

     xlocan(:)=0.0
     xlocan(:)=1.0/(ratio+1)*x(:,index_bcnode(1,ie))+ratio/(ratio+1)*x(:,index_bcnode(2,ie))
     norm_e(:)=norm_con_ele(:,ie)
     tang(1)=norm_e(2)
     tang(2)=-norm_e(1)       


       nit=1
       err_p=999.0
       delta(:)=0.0
       xlocan_temp(:)=xlocan(:)
       do while((nit.lt.5).and.(err_p.gt.1.0e-8))
          xlocan(:)=xlocan(:)+delta(:)

	  call get_indicator_derivative_2D_1st(xlocan,x_inter,x_center,hg,&
	       I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
	  delta(1)=(0.5-II)*dI(1)/(dI(1)**2+dI(2)**2)
	  delta(2)=(0.5-II)*dI(2)/(dI(1)**2+dI(2)**2)

	  err_p=abs(0.5-II)
	  nit=nit+1
       end do
!write(*,*)'========================================='
!write(*,*)'err_p=',err_p,'xlcan=',xlocan(:),'xc=',xlocan_temp(:)
!write(*,*)'ratio=',ratio
       if(err_p.gt.1.0e-6) then
          write(*,*)'something wrong in contact line'
	  goto 200
	end if
       temp=(xlocan_temp(1)-xlocan(1))**2+(xlocan_temp(2)-xlocan(2))**2
       temp=sqrt(temp)
       if(temp.gt.max_hg/1.5) then
	  write(*,*)'tow long distance in contact line'
	  goto 200
       end if

       anglet=acos(norm_p(1)*norm_e(1)+norm_p(2)*norm_e(2))/pi*180
write(*,*)'diff_angle=',anglet-static_angle, myid
if(abs(anglet-static_angle).lt.ad_re_angle) goto 200
       call regen_contact(xlocan,x,x_inter,x_center,hg,I_fluid_center,corr_Ip,nn_org_regen,x_org_regen, &
         		  gridy,norm_e,tang,nn_con_ele,con_ele,ne_regen,ele_regen,ien,norm_p,nn_bottom,x_bottom)
	call find_ca(xlocan,x,vel_fluid,vol_nn,ca,norm_p,theta,anglet,flag_ca)


       do j=1,nn_org_regen
          nn_con_regen=nn_con_regen+1
          x_con_regen(:,nn_con_regen)=x_org_regen(:,j)
       end do

if(flag_ca==0) then
  do j=1,nn_bottom
     nn_con_regen=nn_con_regen+1
     x_con_regen(:,nn_con_regen)=x_bottom(:,j)
  end do
goto 200
end if



!       call find_ca(xlocan,x,vel_fluid,vol_nn,ca,norm_p,theta) 
      x_fix(2)=norm_e(1)*xlocan(1)+norm_e(2)*xlocan(2)
      x_fix(1)=tang(1)*xlocan(1)+tang(2)*xlocan(2)
      norm_fix(2)=norm_e(1)*norm_p(1)+norm_e(2)*norm_p(2)
      norm_fix(1)=tang(1)*norm_p(1)+tang(2)*norm_p(2)
      x_wall(2)=x_fix(2)-gridy


       if(norm_fix(1).lt.0.0) then
         norm_wall(1)=cos(theta+pi/2.0)
	 norm_wall(2)=sin(theta+pi/2.0)
       else
         norm_wall(1)=cos(pi/2.0-theta)
	 norm_wall(2)=sin(pi/2.0-theta)
       end if
       temp_n=nn_con_regen
       call solve_contact(x_fix,x_wall,norm_fix,norm_wall,x_con_regen,nn_con_regen)
       do j=temp_n+1,nn_con_regen
	  rot_x(1)=norm_e(1)*x_con_regen(2,j)+tang(1)*x_con_regen(1,j)
	  rot_x(2)=norm_e(2)*x_con_regen(2,j)+tang(2)*x_con_regen(1,j)
	  x_con_regen(:,j)=rot_x(:)
       end do

200 continue      

  end do ! end of loop over inter elements


777 continue
flag_contact=0
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(flag_contact_temp,flag_contact,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)
if(flag_contact==0) goto 999


local_nn_temp(:)=0
local_nn(:)=0
local_nn_temp(myid+1)=nn_con_regen
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(local_nn_temp(1),local_nn(1),ncpus,mpi_integer,mpi_sum,mpi_comm_world,ierror)
lower=0
do icount=1,myid
   lower=lower+local_nn(icount)
end do
x_con_regen_temp(:,:)=0.0
do icount=1,nn_con_regen
   x_con_regen_temp(:,lower+icount)=x_con_regen(:,icount)
end do
x_con_regen(:,:)=0.0
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(x_con_regen_temp(1,1),x_con_regen(1,1),nsd*maxmatrix,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
nn_con_regen=0
do icount=1,ncpus
   nn_con_regen=nn_con_regen+local_nn(icount)
end do

local_nn_temp(:)=0
local_nn(:)=0
local_nn_temp(myid+1)=ne_regen
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(local_nn_temp(1),local_nn(1),ncpus,mpi_integer,mpi_sum,mpi_comm_world,ierror)
lower=0
do icount=1,myid
   lower=lower+local_nn(icount)
end do
ele_regen_temp(:)=0
do icount=1,ne_regen
   ele_regen_temp(lower+icount)=ele_regen(icount)
end do
ele_regen(:)=0
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(ele_regen_temp(1),ele_regen(1),100,mpi_integer,mpi_sum,mpi_comm_world,ierror)
ne_regen=0
do icount=1,ncpus
   ne_regen=ne_regen+local_nn(icount)
end do






  nn_inter_outer=0
  do i=1,nn_inter
     do j=1,ne_regen
        if(infdomain_inter(i)==ele_regen(j)) goto 666
     end do
     nn_inter_outer=nn_inter_outer+1
     x_inter_outer(:,nn_inter_outer)=x_inter(:,i)
666 continue
  end do


  do i=1,nn_inter_outer
     x_inter(:,i)=x_inter_outer(:,i)
  end do

  nn_inter=nn_inter_outer
!nn_inter=0

  do i=1,nn_con_regen
!     if(x_con_regen(2,i).gt.bot_y) then

       nn_inter=nn_inter+1
       x_inter(:,nn_inter)=x_con_regen(:,i)

!     end if
  end do
999 continue

end subroutine contactline



















