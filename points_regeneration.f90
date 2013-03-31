
! regenerate points for 2D and 3D. 

subroutine points_regeneration(x,x_inter,x_center,x_inter_regen,nn_inter_regen,I_fluid_center,corr_Ip,hg,ien,intflag,&
				regen_ele,ne_regen_ele)

  use fluid_variables, only:nsd,ne,nn,nen
  use interface_variables
  use mpi_variables
  include 'mpif.h'

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,nn_center)
  real(8) x_inter_regen(nsd,maxmatrix)
  integer nn_inter_regen
  real(8) I_fluid_center(nn_center),corr_Ip(maxmatrix),hg(ne)
  integer ien(nen,ne),intflag
  integer ne_regen_ele, regen_ele(ne_regen_ele)


  integer i,j,k,isd,jsd,ksd,inl,node,ie,icount,jcount,kcount
  real(8) x_fluid(nsd,nen) ! global coordinates of elements' nodes
  real(8) xlocan(nsd),xlocan_temp(nsd)
  
  integer, parameter :: nn_sub=5
  integer, parameter :: nn_ele=nn_sub**3
  real(8) sh(nen,nn_ele),x_loc_can(nsd,nn_ele)

  integer nit,flag_loc
  real(8) err_p,temp,delta(nsd)
  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) norm_p(nsd),curv_p,dcurv

  integer index_regen_nn(ncpus),index_regen_nn_temp(ncpus)
  integer nn_can
  real(8),dimension(:,:),allocatable ::x_can,x_can_temp

  integer nn_loc,base,top,loc_index,lower


  if((nsd==2).and.(nen==4)) then
          do i=1,nn_sub
             do j=1,nn_sub
		node=nn_sub*(i-1)+j
                x_loc_can(1,nn_sub*(i-1)+j)=2.0/real(nn_sub)*i-1.0-1.0/real(nn_sub)
                x_loc_can(2,nn_sub*(i-1)+j)=2.0/real(nn_sub)*j-1.0-1.0/real(nn_sub)

             sh(1,node)=0.25*(1-x_loc_can(1,node))*(1-x_loc_can(2,node))
             sh(2,node)=0.25*(1+x_loc_can(1,node))*(1-x_loc_can(2,node))
             sh(3,node)=0.25*(1+x_loc_can(1,node))*(1+x_loc_can(2,node))
             sh(4,node)=0.25*(1-x_loc_can(1,node))*(1+x_loc_can(2,node))

             end do
          end do
  end if
  if((nsd==3).and.(nen==8)) then
           do i=1,nn_sub
              do j=1,nn_sub
                 do k=1,nn_sub
		node=nn_sub*(nn_sub*(i-1)+j-1)+k
                    x_loc_can(1,nn_sub*(nn_sub*(i-1)+j-1)+k)=2.0/real(nn_sub)*i-1.0-1.0/real(nn_sub)
                    x_loc_can(2,nn_sub*(nn_sub*(i-1)+j-1)+k)=2.0/real(nn_sub)*j-1.0-1.0/real(nn_sub)
                    x_loc_can(3,nn_sub*(nn_sub*(i-1)+j-1)+k)=2.0/real(nn_sub)*k-1.0-1.0/real(nn_sub)
             sh(1,node)=0.125*(1-x_loc_can(1,node))*(1-x_loc_can(2,node))*(1-x_loc_can(3,node))
             sh(2,node)=0.125*(1+x_loc_can(1,node))*(1-x_loc_can(2,node))*(1-x_loc_can(3,node))
             sh(3,node)=0.125*(1+x_loc_can(1,node))*(1+x_loc_can(2,node))*(1-x_loc_can(3,node))
             sh(4,node)=0.125*(1-x_loc_can(1,node))*(1+x_loc_can(2,node))*(1-x_loc_can(3,node))
             sh(5,node)=0.125*(1-x_loc_can(1,node))*(1-x_loc_can(2,node))*(1+x_loc_can(3,node))
             sh(6,node)=0.125*(1+x_loc_can(1,node))*(1-x_loc_can(2,node))*(1+x_loc_can(3,node))
             sh(7,node)=0.125*(1+x_loc_can(1,node))*(1+x_loc_can(2,node))*(1+x_loc_can(3,node))
             sh(8,node)=0.125*(1-x_loc_can(1,node))*(1+x_loc_can(2,node))*(1+x_loc_can(3,node))

                 end do
              end do
           end do
  end if

  if(ne_regen_ele .le. ncpus) then
    if(myid+1 .le. ne_regen_ele) then
	nn_loc=1
    else
	nn_loc=0
    end if
  else
    base=floor(real(ne_regen_ele)/real(ncpus))
    top=ne_regen_ele-base*ncpus
    if(myid+1 .le. top) then
	nn_loc=base+1
    else
	nn_loc=base
    end if
  end if
  nn_inter_regen=0
  x_inter_regen(:,:)=0.0
  do loc_index=1,nn_loc
     i=myid+1+(loc_index-1)*ncpus

     ie=regen_ele(i)
     do inl=1,nen
	node=ien(inl,ie)
	x_fluid(1:nsd,inl)=x(1:nsd,node)
     end do
    
     do icount=1,nn_sub**nsd
        xlocan(:)=0.0
        do inl=1,nen
	   xlocan(1:nsd)=xlocan(1:nsd)+sh(inl,icount)*x_fluid(1:nsd,inl)
	end do

        if(nsd==2) then
           call get_indicator_derivative_2D_1st(xlocan,x_inter,x_center,hg,&
                                I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
        else if(nsd==3) then
           call get_indicator_derivative_3D_1st(xlocan,x_inter,x_center,hg,&
                                I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
        end if

	if(intflag==1) then
          if( (abs(II-0.5).lt.0.08) .and. (abs(II-0.5).gt.0.02)) then
	    nn_inter_regen=nn_inter_regen+1
	    x_inter_regen(1:nsd,nn_inter_regen)=xlocan(1:nsd)
	  end if
	elseif(intflag==2) then
	  if( (II.gt.0.5+0.02) .and. (II.lt.0.5+0.08)) then
	    nn_inter_regen=nn_inter_regen+1
	    x_inter_regen(1:nsd,nn_inter_regen)=xlocan(1:nsd)
	  end if
	elseif(intflag==3) then
	  if( (II.lt.0.5-0.02) .and. (II.gt.0.5-0.08)) then
	    nn_inter_regen=nn_inter_regen+1
	    x_inter_regen(1:nsd,nn_inter_regen)=xlocan(1:nsd)
	  end if
	end if
     end do  ! end of icount
  end do ! end of loc_index   
       
  index_regen_nn(:)=0
  index_regen_nn_temp(:)=0
  index_regen_nn_temp(myid+1)=nn_inter_regen

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(index_regen_nn_temp(1),index_regen_nn(1),ncpus,mpi_integer,mpi_sum,mpi_comm_world,ierror)

  nn_can=0
  do icount=1,ncpus
     nn_can=nn_can+index_regen_nn(icount)
  end do

  if(allocated(x_can)) then
     deallocate(x_can)
  end if
  allocate(x_can(nsd,nn_can))
  allocate(x_can_temp(nsd,nn_can))

  lower=0
  do icount=1,myid
     lower=lower+index_regen_nn(icount)
  end do
  
  x_can(:,:)=0.0
  x_can_temp(:,:)=0.0
  do icount=1,nn_inter_regen
     x_can_temp(1:nsd,lower+icount)=x_inter_regen(1:nsd,icount)
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(x_can_temp(1,1),x_can(1,1),nsd*nn_can,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

  deallocate(x_can_temp)
  call mpi_barrier(mpi_comm_world,ierror)
!goto 900
  nn_inter_regen=0
  x_inter_regen(:,:)=0.0

  if(nn_can .le. ncpus) then
    if(myid+1 .le. nn_can) then
        nn_loc=1
    else
        nn_loc=0
    end if
  else
    base=floor(real(nn_can)/real(ncpus))
    top=nn_can-base*ncpus
    if(myid+1 .le. top) then
        nn_loc=base+1
    else
        nn_loc=base
    end if
  end if

  do loc_index=1,nn_loc
     i=myid+1+(loc_index-1)*ncpus
     xlocan(1:nsd)=x_can(1:nsd,i)
     
     nit=1
     err_p=999
     delta(:)=0.0
     xlocan_temp(1:nsd)=xlocan(1:nsd)

     do while ((nit.le.7).and.(err_p.gt.1.0e-6))
	temp=0.0
	do isd=1,nsd
	   temp=temp+delta(isd)**2
	end do
	temp=sqrt(temp)
	if(temp.gt.max_hg) delta(:)=delta(:)/temp*max_hg
	xlocan(:)=xlocan(:)+delta(:)
	if(nsd==2) then
          call get_indicator_derivative_2D_1st(xlocan,x_inter,x_center,hg,&
              I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
        else
	  call get_indicator_derivative_3D_1st(xlocan,x_inter,x_center,hg,&
              I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
	end if
	temp=0.0
	do isd=1,nsd
	   temp=temp+dI(isd)**2
	end do
	delta(:)=(0.5-II)*dI(:)/temp
	err_p=abs(II-0.5)
	nit=nit+1
    end do

    flag_loc=0
    if(err_p.lt.1.0e-6) then
       if(nsd==2) then
              call get_curv_num_2D(xlocan,x_inter,x_center,hg,I_fluid_center,&
                        corr_Ip,dcurv,curv_p,norm_p)
       else if(nsd==3) then
              call get_curv_num_3D(xlocan,x_inter,x_center,hg,I_fluid_center,&
                        corr_Ip,dcurv,curv_p,norm_p)

       end if
       if((dcurv.lt.maxdcurv)) flag_loc=1
    end if

    if(flag_loc==1) then
!	do icount=1,nn_inter_regen
!	   temp=0.0
!	   do isd=1,nsd
!	      temp=temp+(x_inter_regen(isd,icount)-xlocan(isd))**2
!	   end do
!	   temp=sqrt(temp)
!	   if(temp.lt.max_hg/real(nn_sub+1)) goto 200
!	end do
	nn_inter_regen=nn_inter_regen+1
	x_inter_regen(:,nn_inter_regen)=xlocan(:)

200 continue
    end if

  end do ! end of loop loc_index


  index_regen_nn(:)=0
  index_regen_nn_temp(:)=0
  index_regen_nn_temp(myid+1)=nn_inter_regen
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(index_regen_nn_temp(1),index_regen_nn(1),ncpus,mpi_integer,mpi_sum,mpi_comm_world,ierror)

  nn_can=0
  do icount=1,ncpus
     nn_can=nn_can+index_regen_nn(icount)
  end do
  if(allocated(x_can)) then
     deallocate(x_can)
  end if
  allocate(x_can(nsd,nn_can))
  allocate(x_can_temp(nsd,nn_can))

  lower=0
  do icount=1,myid
     lower=lower+index_regen_nn(icount)
  end do
  x_can(:,:)=0.0
  x_can_temp(:,:)=0.0
  do icount=1,nn_inter_regen
     x_can_temp(1:nsd,lower+icount)=x_inter_regen(1:nsd,icount)
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(x_can_temp(1,1),x_can(1,1),nsd*nn_can,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

  deallocate(x_can_temp)

900 continue
  nn_inter_regen=nn_can
  x_inter_regen(:,1:nn_inter_regen)=x_can(:,1:nn_can)

end subroutine points_regeneration



























 
