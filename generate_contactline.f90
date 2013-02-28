

subroutine generate_contactline(x_inter_regen,nn_inter_regen,x,x_inter,x_center,hg,I_fluid_center, &
			corr_Ip,ien,nn_con_ele,con_ele,norm_con_ele,vel_fluid,&
				vol_nn,infdomain_inter)

  use fluid_variables,only:nsd,ne,nn,nen,ndf
  use interface_variables
!  use allocate_variables
  use mpi_variables
  include 'mpif.h'

  real(8) x_inter_regen(nsd,maxmatrix),x(nsd,nn)
  real(8) x_inter(nsd,maxmatrix),x_center(nsd,nn_center),hg(ne)
  integer nn_inter_regen
  real(8) I_fluid_center(nn_center),corr_Ip(maxmatrix)
  integer ien(nen,ne),nn_con_ele,con_ele(nn_con_ele)
  real(8) norm_con_ele(nsd,nn_con_ele)
  real(8) vel_fluid(nsd,nn),vol_nn(nn)
  integer infdomain_inter(maxmatrix)


  integer nn_loc,base,top,loc_index

  integer i,j,icount,jcount!,nit

  integer finf,index_con

  real(8) xc(nsd),norm_c(nsd),norm_e(nsd),xr(nsd),norm_r(nsd),curv_c
  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) xlocan(nsd),xlocan_temp(nsd),err_p!,delta(nsd)
  
  integer nn_contact_proc
  real(8) x_contact_proc(nsd,maxmatrix)
  real(8) norm_g(nsd),curv_g !not using ,just something dummy

  real(8) interval,tangx(nsd),tangy(nsd),modmod,temp

  integer nitt,flag
  real(8) deltaa(nsd)

  integer NumRegen
  integer IndexRef
  real(8) x_cor_regen(nsd,30) !30=maxconn

  integer ne_regen
  integer ele_regen(nn_con_ele)

  real(8) anglet,thelta,ca
  integer flag_ca

  real(8) x_fix(2),x_wall(2),norm_fix(2),norm_wall(2)
  integer nn_2d
  real(8) x_2d(2,30)

  integer local_nn(ncpus),local_nn_temp(ncpus)
  integer nn_contact_temp,lower
  real(8) x_contact_temp(nsd,maxmatrix)
  integer ne_ele_regen_temp,ele_regen_temp(nn_con_ele)
  real(8) vec_contact(nsd)


  NumRegen=0
  IndexRef=0
  x_cor_regen(:,:)=0.0
  ele_regen(:)=0
  ne_regen=0
  nn_contact_proc=0
  x_contact_proc(:,:)=0.0
  maxconn=30
  interval=hg_sp/5.0

  if(nn_inter_regen.le.ncpus) then
    if(myid+1.le.nn_inter_regen) then
      nn_loc=1
    else 
      nn_loc=0
    end if
  else
    base=floor(real(nn_inter_regen)/real(ncpus))
    top=nn_inter_regen-base*ncpus
    if(myid+1.le.top) then
	nn_loc=base+1
    else
	nn_loc=base
    end if
  end if

  icount=0

  do loc_index=1,nn_loc
     i=myid+1+(loc_index-1)*ncpus  ! id
     xc(:)=x_inter_regen(:,i)
     call getinf_el_3d_con(finf,xc,x,nn,nsd,ne,nen,ien,maxconn,nn_con_ele,con_ele,index_con)
     norm_e(:)=norm_con_ele(:,index_con)
     call get_indicator_derivative_3D(xc,x_inter,x_center,hg,I_fluid_center,corr_Ip,II,dI,ddI,norm_c,curv_c)
     xlocan(:)=xc(:)
     icount=1
     x_cor_regen(:,:)=0.0
     x_cor_regen(:,icount)=xc(:)
     NumRegen=0
     IndexRef=0
     do j=1,12
        nitt=1
	err_p=999.0
	deltaa(:)=0.0
	xlocan(1:nsd)=xlocan(1:nsd)+interval*norm_e(1:nsd)
	xlocan_temp(:)=xlocan(:)
	do while((nitt.le.5).and.(err_p.gt.1.0e-6))
	   temp=sqrt(deltaa(1)**2+deltaa(2)**2+deltaa(3)**2)
	   if(temp.gt.hg_sp) deltaa(:)=deltaa(:)/temp*hg_sp
	   xlocan(:)=xlocan(:)+deltaa(:)
	   call get_indicator_derivative_3D_1st(xlocan,x_inter,x_center,hg,I_fluid_center,corr_Ip,II,dI,ddI,norm_g,curv_g)
           tangx(1)=norm_e(2)*norm_g(3)-norm_e(3)*norm_g(2)
           tangx(2)=norm_e(3)*norm_g(1)-norm_e(1)*norm_g(3)
           tangx(3)=norm_e(1)*norm_g(2)-norm_e(2)*norm_g(1)

           tangy(1)=tangx(2)*norm_e(3)-tangx(3)*norm_e(2)
           tangy(2)=tangx(3)*norm_e(1)-tangx(1)*norm_e(3)
           tangy(3)=tangx(1)*norm_e(2)-tangx(2)*norm_e(1)

	   modmod=sqrt(tangy(1)**2+tangy(2)**2+tangy(3)**2)
	   tangy(:)=tangy(:)/modmod
	   temp=dI(1)*tangy(1)+dI(2)*tangy(2)+dI(3)*tangy(3)
	  if(abs(temp).lt.1.0e-6)write(*,*)'something wrong in generate contactline'
	   deltaa(:)=(0.5-II)*tangy(:)/temp
	   nitt=nitt+1
	   err_p=abs(II-0.5)
	end do
	if(err_p.lt.1.0e-6) then
	  call getinf_el_3d_con(finf,xlocan,x,nn,nsd,ne,nen,ien,maxconn,nn_con_ele,con_ele,index_con)
!	  if(finf==0) goto 200
	  
	  icount=icount+1
	  x_cor_regen(:,icount)=xlocan(:)
	  NumRegen=icount
	  if( (xlocan(1)-xc(1))*norm_e(1)+(xlocan(2)-xc(2))*norm_e(2)+(xlocan(3)-xc(3))*norm_e(3).lt.1.9*hg_sp) then
	      IndexRef=icount
	      norm_r(:)=norm_g(:)
	      vec_contact(:)=tangy(:) !tangential direction of wall, used to cal contact vel
	      xr(:)=xlocan(:)
	  end if

	  flag=1
	  do jcount=1,ne_regen
	     if(finf==ele_regen(jcount))flag=0
	  end do
	  if(flag==1) then
	    ne_regen=ne_regen+1
	    ele_regen(ne_regen)=finf
	  end if

	else
	  xlocan(:)=xlocan_temp(:) !if not generate, perform for the next point
	end if


      end do ! end of loop over j


200 continue

!           call get_indicator_derivative_3D(x_cor_regen(:,3),x_inter,x_center,hg,I_fluid_center,corr_Ip,II,dI,ddI,norm_g,curv_g)

      anglet=acos(norm_c(1)*norm_e(1)+norm_c(2)*norm_e(2)+norm_c(3)*norm_e(3))/3.14159*180.0
write(*,*)anglet,curv_c
      if(abs(anglet-static_angle).lt.ad_re_angle) then
	 do jcount=1,NumRegen
	    nn_contact_proc=nn_contact_proc+1
	    x_contact_proc(:,nn_contact_proc)=x_cor_regen(:,jcount)
	 end do
      else
!	 write(*,*)'vec_contact=',vec_contact(:)
	 call find_ca(xr,x,vel_fluid,vol_nn,ca,vec_contact,thelta,anglet,flag_ca)
	 if(flag_ca==-1) then
	    do jcount=1,NumRegen
	       nn_contact_proc=nn_contact_proc+1
	       x_contact_proc(:,nn_contact_proc)=x_cor_regen(:,jcount)
	    end do
	  else
	    do jcount=IndexRef+1,NumRegen
	       nn_contact_proc=nn_contact_proc+1
	       x_contact_proc(:,nn_contact_proc)=x_cor_regen(:,jcount)
	    end do

            tangy(1)=norm_e(2)*norm_r(3)-norm_e(3)*norm_r(2)
	    tangy(2)=norm_e(3)*norm_r(1)-norm_e(1)*norm_r(3)
	    tangy(3)=norm_e(1)*norm_r(2)-norm_e(2)*norm_r(1)

	    modmod=sqrt(tangy(1)**2+tangy(2)**2+tangy(3)**2)
	    tangy(:)=tangy(:)/modmod

            tangx(1)=tangy(2)*norm_e(3)-tangy(3)*norm_e(2)
            tangx(2)=tangy(3)*norm_e(1)-tangy(1)*norm_e(3)
            tangx(3)=tangy(1)*norm_e(2)-tangy(2)*norm_e(1)
	    modmod=sqrt(tangx(1)**2+tangx(2)**2+tangx(3)**2)
	    tangx(:)=tangx(:)/modmod

	    temp=(xr(1)-xc(1))*norm_e(1)+(xr(2)-xc(2))*norm_e(2)+(xr(3)-xc(3))*norm_e(3)
!	    write(*,*)'distcance=',temp,'ca=',ca,'thelta=',thelta/3.14159*180.0
	    thelta=3.14159/2.0*120.0/90.0
	    x_fix(1)=0.0
	    x_fix(2)=0.0
	    x_wall(1)=0.0
	    x_wall(2)=-temp!0.5*max_hg
	    norm_fix(1)=norm_r(1)*tangx(1)+norm_r(2)*tangx(2)+norm_r(3)*tangx(3)
	    norm_fix(2)=norm_r(1)*norm_e(1)+norm_r(2)*norm_e(2)+norm_r(3)*norm_e(3)
	    norm_wall(1)=cos(3.14159/2.0-thelta)
	    norm_wall(2)=sin(3.14159/2.0-thelta)

	    nn_2d=0
	    x_2d(:,:)=0.0
	    call solve_contact(x_fix,x_wall,norm_fix,norm_wall,x_2d,nn_2d)
	    do jcount=1,nn_2d
	       nn_contact_proc=nn_contact_proc+1
	       x_contact_proc(1,nn_contact_proc)=tangx(1)*x_2d(1,jcount)+tangy(1)*0.0+norm_e(1)*x_2d(2,jcount)
	       x_contact_proc(2,nn_contact_proc)=tangx(2)*x_2d(1,jcount)+tangy(2)*0.0+norm_e(2)*x_2d(2,jcount)
	       x_contact_proc(3,nn_contact_proc)=tangx(3)*x_2d(1,jcount)+tangy(2)*0.0+norm_e(3)*x_2d(2,jcount)
	       x_contact_proc(:,nn_contact_proc)=x_contact_proc(:,nn_contact_proc)+xr(:)
	    end do
	  end if ! end if of flag=0
      end if!end if of abs(anglt-static)

  end do  ! end of loop over loc_index




  local_nn(:)=0
  local_nn_temp(:)=0
  local_nn_temp(myid+1)=nn_contact_proc

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(local_nn_temp(1),local_nn(1),ncpus,mpi_integer,mpi_sum,mpi_comm_world,ierror)
  lower=0
  do icount=1,myid
     lower=lower+local_nn(icount)
  end do
  x_contact_temp(:,:)=0.0
  do icount=1,nn_contact_proc
     x_contact_temp(:,lower+icount)=x_contact_proc(:,icount)
  end do
  x_contact_proc(:,:)=0.0
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(x_contact_temp(1,1),x_contact_proc(1,1),nsd*maxmatrix,mpi_double_precision, &
		mpi_sum,mpi_comm_world,ierror)
  nn_contact_proc=0
  do icount=1,ncpus
     nn_contact_proc=nn_contact_proc+local_nn(icount)
  end do

  local_nn(:)=0
  local_nn_temp(:)=0
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
  call mpi_allreduce(ele_regen_temp(1),ele_regen(1),nn_con_ele,mpi_integer,mpi_sum,mpi_comm_world,ierror)
  ne_regen=0
  do icount=1,ncpus
     ne_regen=ne_regen+local_nn(icount)
  end do

  ne_regen_temp=ne_regen
  ele_regen_temp(:)=ele_regen(:)
  
  ne_regen=0
  do icount=1,ne_regen_temp
     flag=1
     do jcount=1,ne_regen
        if(ele_regen_temp(icount)==ele_regen(jcount)) flag=0
     end do
     if(flag==1) then
	ne_regen=ne_regen+1
	ele_regen(ne_regen)=ele_regen_temp(icount)
     end if
  end do

open(222,file='points_contact.dat',status='unknown')
if(myid==0) then
do i=1,nn_contact_proc
   write(222,*)x_contact_proc(:,i)
end do
end if



open(444,file='initialpoints.dat',status='unknown')
if(myid==0) then
do i=1,nn_inter
   write(444,*) x_inter(:,i)
end do
end if
close(444)

  do icount=1,nn_inter
     flag=0
!     do jcount=1,ne_regen
!        if(infdomain_inter(icount)==ele_regen(jcount)) flag=1
!     end do



!     do jcount=1,nn_con_ele
!	if(infdomain_inter(icount)==con_ele(jcount)) flag=1
!     end do


!     if(x_inter(3,icount).le.2.0*max_hg) flag=1
     temp=sqrt(x_inter(1,icount)**2+(x_inter(2,icount)-1.5)**2)
     if( (temp.gt.1.5*1.414-2.0*hg_sp).and.(x_inter(2,icount).lt.0.0) ) flag=1


     if(flag==0) then
        nn_contact_proc=nn_contact_proc+1
        x_contact_proc(:,nn_contact_proc)=x_inter(:,icount)
     end if
  end do


  nn_inter=nn_contact_proc
  x_inter(:,1:nn_inter)=x_contact_proc(:,1:nn_inter)

!open(222,file='points_contact.dat',status='unknown')
!if(myid==0) then
!do i=1,nn_contact_proc
!   write(222,*)x_contact_proc(:,i)
!end do
!end if





999 continue
end subroutine generate_contactline
