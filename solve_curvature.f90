

subroutine solve_curvature(x,x_inter,x_center,I_fluid,corr_Ip,I_fluid_center,curv_nn,hg,ien,&
			ne_local,ien_local,node_local,nn_local, &
		global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length,&
		sur_fluid,flag_domain)


  use fluid_variables, only:nsd,ne,nn,nen,ne_spbc
  use interface_variables
  use allocate_variables, only:center_domain,nn_center_domain,inter_ele,ne_inter
  use mpi_variables
  include 'mpif.h'

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,nn_center),I_fluid(nn)
  real(8) corr_Ip(maxmatrix),I_fluid_center(nn_center),curv_nn(nn),hg(ne)
  integer ien(nen,ne)
  real(8) sur_fluid(nsd,nn)

  integer ne_local,ien_local(ne_local)
  integer node_local(nn_local),nn_local
  integer nn_global_com,global_com(nn_global_com)
  integer nn_local_com,local_com(nn_local_com)
  integer ad_length,send_address(ad_length,2)



  real(8) II,dI(nsd),ddI(3*(nsd-1)),curv_a,norm_a(nsd)
  integer i,j,nit,ie,node,isd
  real(8) delta(nsd),xlocan(nsd),temp,err_p,dis

  integer nn_loc,base,top,loc_index

  integer flag_node(nn)

  integer nn_node_inter
  integer node_inter(4*ne_inter),fcurv_node_inter(4*ne_inter) !indicate the node with calculated curvature
  integer fcurv_node_inter_temp(4*ne_inter)

  real(8) curv_nn_temp(nn),dg(nn)
  integer id_curv(nn)


!=====================!
  real(8) w(nn),p(nn)
  integer flag_domain(ne)


!  curv_nn_old(:)=curv_nn(:)

  id_curv(:)=0
  curv_nn(:)=0.0
  curv_nn_temp(:)=0.0
  nn_node_inter=0
  flag_node(:)=0
  fcurv_node_inter(:)=0
  fcurv_node_inter_temp(:)=0
  do i=1,ne_inter
     ie=inter_ele(i)
     do j=1,nen
	node=ien(j,ie)
	flag_node(node)=1
     end do
  end do

  do i=1,nn
     if(flag_node(i)==1) then
	nn_node_inter=nn_node_inter+1
	node_inter(nn_node_inter)=i
     end if
  end do
  if(nn_node_inter.le.ncpus) then
    if(myid+1.le.nn_node_inter) then
      nn_loc=1
    else
      nn_loc=0
    end if
  else
    base=floor(real(nn_node_inter)/real(ncpus))
    top=nn_node_inter-base*ncpus
    if(myid+1.le.top) then
      nn_loc=base+1
    else
      nn_loc=base
    end if
  end if

  do loc_index=1,nn_loc
     i=node_inter(myid+1+(loc_index-1)*ncpus)
     nit=1
     err_p=999.0
     xlocan(:)=x(:,i)
!if(x(3,i).lt.0.001) goto 888
     delta(:)=0.0
     do while((nit.le.5).and.(err_p.gt.1.0e-6))
	temp=0.0
        do isd=1,nsd
	   temp=temp+delta(isd)**2
	end do
	temp=sqrt(temp)
	if(temp.gt.max_hg) delta(:)=delta(:)/temp*max_hg
	xlocan(:)=xlocan(:)+delta(:)
	if(nsd==2) then
	   call get_indicator_derivative_2D_1st(xlocan,x_inter,x_center,hg, &
		I_fluid_center,corr_Ip,II,dI,ddI,norm_a,curv_a)
	else
           call get_indicator_derivative_3D_1st(xlocan,x_inter,x_center,hg, &
                I_fluid_center,corr_Ip,II,dI,ddI,norm_a,curv_a)
	end if
	temp=0.0
	do isd=1,nsd
	   temp=temp+dI(isd)**2
	end do
	delta(:)=(0.5-II)*dI(:)/temp

	err_p=abs(II-0.5)
	nit=nit+1
     end do ! end of do while
     if(err_p.lt.1.0e-6) then
	if(nsd==2) then
	   call get_indicator_derivative_2D(xlocan,x_inter,x_center,hg,I_fluid_center, &
			corr_Ip,II,dI,ddI,norm_a,curv_a)
	else
           call get_indicator_derivative_3D(xlocan,x_inter,x_center,hg,I_fluid_center, &
                        corr_Ip,II,dI,ddI,norm_a,curv_a)
	end if
	dis=0.0
	do isd=1,nsd
	   dis=dis+(x(isd,i)-xlocan(isd))**2
	end do
	dis=sqrt(dis)
	curv_nn_temp(i)=-curv_a
	fcurv_node_inter_temp(myid+1+(loc_index-1)*ncpus)=1
      end if
888 continue
  end do  ! end loop over loc_index


  curv_nn(:)=0.0
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(curv_nn_temp(1),curv_nn(1),nn,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
  call mpi_allreduce(fcurv_node_inter_temp(1),fcurv_node_inter(1),4*ne_inter,mpi_integer,mpi_sum,mpi_comm_world,ierror)

!  do i=1,nn_node_inter
!     if(fcurv_node_inter(i)==1) then
!	j=node_inter(i)
!	curv_nn_old(j)=curv_nn(j)
!     end if
!  end do

!  curv_nn(:)=curv_nn_old(:)
!goto 123
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!   solve curvature using \nabla I \cdot \nabla k =0
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  flag_domain(:)=0
  id_curv(:)=0

!  call find_curv_domain(x,x_inter,ne_local,ien_local,flag_domain,id_curv,ien)
!  do i=1,nn_center_domain
!     ie=center_domain(i)
   do ie=1,ne
!   if(ie.le.ne) then
!     flag_domain(ie)=1
!
!
    if(flag_domain(ie)==1) then
     do j=1,nen
	node=ien(j,ie)
	id_curv(node)=1
     end do
    end if
  end do

  do i=1,nn_node_inter
     j=node_inter(i)
     if(fcurv_node_inter(i)==1)id_curv(j)=0 !set bc id
  end do


  p(:)=0.0
  w(:)=0.0

  call block_curv(flag_domain,x,sur_fluid,curv_nn,p,w,ien,ne_local,ien_local,I_fluid)
  call mpi_barrier(mpi_comm_world,ierror)

  call communicate_res_ad(p,1,nn,send_address,ad_length)
  call communicate_res_ad(w,1,nn,send_address,ad_length)

  call setid_pa(p,1,nn,id_curv,node_local,nn_local)
  dg(:)=0.0
do i=1,nn
   if(abs(w(i)).lt.1.0e-5) w(i)=1.0
end do

  call gmrescurv(x,w,p,dg,hg,ien,id_curv,&
                  ne_local,ien_local,node_local,nn_local, &
                   global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length,&
		  sur_fluid,flag_domain,I_fluid)
  call setid_pa(dg,1,nn,id_curv,node_local,nn_local)
  curv_nn(:)=curv_nn(:)+dg(:)
!p(:)=0.0
!  call block_curv(flag_domain,x,sur_fluid,curv_nn,p,w,ien,ne_local,ien_local,I_fluid)
!  call communicate_res_ad(p,1,nn,send_address,ad_length)
!temp=0.0
!  call setid_pa(p,1,nn,id_curv,node_local,nn_local)
!  call getnorm_pa(p,1,nn,node_local,nn_local,temp)
!if(myid==0)write(*,*)'temp=',sqrt(temp)

123 continue
end subroutine solve_curvature













