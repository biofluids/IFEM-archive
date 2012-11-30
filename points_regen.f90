!========================================
!  Regenerate interfacial points
!========================================

subroutine points_regen(x,x_inter,x_center,x_inter_regen,nn_inter_regen,&
			I_fluid_center,corr_Ip,hg,ien,intflag)

  use fluid_variables, only:nsd,ne,nn,nen
  use interface_variables
  use allocate_variables, only:regen_ele_loc,ne_regen_ele_loc,&
		regen_ele,ne_regen_ele
  use mpi_variables
  include 'mpif.h'

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,ne)
  real(8) x_inter_regen(nsd,maxmatrix)
  integer nn_inter_regen
  real(8) I_fluid_center(ne),corr_Ip(maxmatrix)
  real(8) hg(ne)
  integer ien(nen,ne)
  integer intflag

  integer i,j,isd,inl,node,ie,icount,jcount
  real(8) x_fluid(nsd,nen) !global coordinates of elements' node
  integer nn_sub,nn_ele  !nn_sub*nn_sub=nn_ele
  real(8) x_loc_can(nsd,maxmatrix) !local coor of candidate points
  real(8) xlocan(nsd),xlocan_temp(nsd)

  real(8) sh(nen) !shape function
  real(8) x_glo(nsd),distance,hs,dx(nsd),Sp
  integer nn_cr,nn_local

  integer nit
  real(8) err_p,temp,delta(nsd)
  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) norm_p(nsd),curv_p,dcurv
  real(8) Ic_inter
  integer nn_inter_regen_loc
  integer nn_inter_ele_loc
  real(8) x_inter_regen_loc(nsd,maxmatrix) !regen points per processor
  real(8) x_inter_ele_loc(nsd,maxmatrix) !regen points per elemment

  integer flag_loc,d_flag

  integer index_regen_nn(ncpus),index_regen_nn_temp(ncpus)
  integer lower
  real(8) x_regen_temp(nsd,maxmatrix)

  nn_inter_regen=0
  x_inter_regen(:,:)=0.0
  nn_inter_regen_loc=0
  x_inter_regen_loc(:,:)=0.0
  nn_cr=10
  nn_local=0
  nn_sub=10
  Ic_inter=0.5
!******************************************************
  do jcount=1,ne_regen_ele_loc !begin element loop
     ie=regen_ele_loc(jcount)

     do inl=1,nen
	node=ien(inl,ie)
	x_fluid(1:nsd,inl)=x(1:nsd,node)
     end do
     nn_local=0
     nn_sub=10
!!!!!!**************************************************
     do while((nn_local.le.nn_cr).and.(nn_sub.le.10)) !begin regeneration loop
	nn_local=0
	x_inter_ele_loc(:,:)=0.0

        if((nsd==2) .and. (nen==4)) then  !2d4n
          do i=1,nn_sub
             do j=1,nn_sub
                x_loc_can(1,nn_sub*(i-1)+j)=2.0/nn_sub*i-1.0-1.0/nn_sub
                x_loc_can(2,nn_sub*(i-1)+j)=2.0/nn_sub*j-1.0-1.0/nn_sub
             end do
          end do
          nn_ele=nn_sub**2
        end if       !assign local coordinates for candidate points(2d4n)
!!!!!!!!!**********************************************
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
           Ic_inter=0.5

           call get_indicator_derivative_2D_1st(xlocan,x_inter,x_center,hg,&
                                I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)

       if(intflag==2) then

!           if(((II-Ic_inter).gt.0.).or.((Ic_inter-II).gt.0.02)) then
	    if( (abs(II-Ic_inter).gt.0.05) .or. (abs(II-Ic_inter).lt.0.02) .or. ((II-Ic_inter).gt.0)) then
                goto 200
           end if
       else if(intflag==3) then

            if( (abs(II-Ic_inter).gt.0.05) .or. (abs(II-Ic_inter).lt.0.02) .or. ((II-Ic_inter).lt.0)) then
                goto 200
           end if
       else if(intflag==1) then
!           if((II.gt.Ic_inter+0.25).or.(II.lt.Ic_inter-0.25).or.(abs(II-Ic_inter).lt.1.0e-3)) then
            if( (abs(II-Ic_inter).gt.0.05) .or. (abs(II-Ic_inter).lt.0.01)) then

!           if(II.ge.Ic_inter) then
		goto 200
	   end if

       end if

!=======================================================!
!++++++++++++point projection+++++++++++++++++++++++++++!
!=======================================================!
	   nit=1
	   err_p=999.0
	   delta(:)=0.0
	   Ic_inter=0.5
	   xlocan_temp(1:nsd)=xlocan(1:nsd)
	   do while((nit.le.5).and.(err_p.gt.1.0e-9))
	      xlocan(1:nsd)=xlocan(1:nsd)+delta(1:nsd)
	      call get_indicator_derivative_2D_1st(xlocan,x_inter,x_center,hg,&
				I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
	      if(II.gt.900) then
	        goto 200
	      end if
	      if(nsd==2) then
	        temp=dI(1)**2+dI(2)**2
	      delta(1)=(Ic_inter-II)*dI(1)/temp
	      delta(2)=(Ic_inter-II)*dI(2)/temp
	      end if
              err_p=abs(II-Ic_inter)
	      nit=nit+1
	   end do
!==========================================================!
	   flag_loc=0
	   if(err_p.lt.1.0e-9) then
	     distance=sqrt((xlocan(1)-xlocan_temp(1))**2+(xlocan(2)-xlocan_temp(2))**2)
!	     if(distance.lt.hg(ie)/nn_sub/1.5) then
	     if(distance.lt.1.0*max_hg) then
		call get_curv_num_2D(xlocan,x_inter,x_center,hg,I_fluid_center,&
			corr_Ip,dcurv,curv_p,norm_p)      
!                   dcurv=0.1
		  if((dcurv.lt.maxdcurv)) then
		    flag_loc=1
		  end if
	     end if
	   end if
	   if(flag_loc==1) then
             do icount=1,nn_local
                temp=(x_inter_ele_loc(1,icount)-xlocan(1))**2+(x_inter_ele_loc(2,icount)-xlocan(2))**2
                temp=sqrt(temp)
                if(temp.lt.hg(ie)/real(nn_cr)) then
                  goto 200
                end if
             end do
	     nn_local=nn_local+1
	     x_inter_ele_loc(1:nsd,nn_local)=xlocan(1:nsd)
	   end if

200 continue

	end do ! end of candidate points loop
	nn_sub=nn_sub+4     
   end do !end of do while regeneration loop
   do icount=1,nn_local
      nn_inter_regen_loc=nn_inter_regen_loc+1
      x_inter_regen_loc(1:nsd,nn_inter_regen_loc)=x_inter_ele_loc(1:nsd,icount)
   end do
300 continue
  end do ! end of element loop for each processor
  index_regen_nn(:)=0
  index_regen_nn_temp(:)=0
  index_regen_nn_temp(myid+1)=nn_inter_regen_loc
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(index_regen_nn_temp(1),index_regen_nn(1),ncpus,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
  call mpi_bcast(index_regen_nn(1),ncpus,mpi_integer,0,mpi_comm_world,ierror)
  nn_inter_regen=0
  call mpi_barrier(mpi_comm_world,ierror)
  
  do icount=1,ncpus
     nn_inter_regen=nn_inter_regen+index_regen_nn(icount)
  end do

  lower=0
  do icount=1,myid
     lower=lower+index_regen_nn(icount)
  end do
  x_regen_temp(:,:)=0.0
  do icount=1,nn_inter_regen_loc
     x_regen_temp(1:nsd,lower+icount)=x_inter_regen_loc(1:nsd,icount)
  end do
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(x_regen_temp(1,1),x_inter_regen(1,1),nsd*maxmatrix,mpi_double_precision, &
		mpi_sum,0,mpi_comm_world,ierror)
  call mpi_bcast(x_inter_regen(1,1),nsd*maxmatrix,mpi_double_precision,0,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

if(myid==0) then
  write(*,*)'nn_inter_regen=',nn_inter_regen
end if
111 format(f14.10,f14.10,f14.10)

end subroutine points_regen













