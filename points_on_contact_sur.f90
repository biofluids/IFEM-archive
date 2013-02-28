

subroutine points_on_contact_sur(x,x_inter,x_center,I_fluid,I_fluid_center,hg,corr_Ip,ien, &
			nn_con_ele,con_ele,norm_con_ele,rngface,nn_inter_regen,x_inter_regen)

  use fluid_variables, only:nsd,nn,ne,nen,neface,nnface,nrng,mapping,etype
  use interface_variables
!  use allocate_variables
  use mpi_variables
  include 'mpif.h'

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,nn_center),I_fluid(nn),I_fluid_center(nn_center),hg(ne)
  real(8) corr_Ip(maxmatrix)
  integer ien(nen,ne)

  integer nn_con_ele,con_ele(nn_con_ele)
  real(8) norm_con_ele(nsd,nn_con_ele)
  integer rngface(neface,ne)

  integer node(nnface)
  real(8) xloc(nsd,nnface),I_node(nnface)

  integer nn_regen_proc
  real(8) x_regen_proc(nsd,500)

  integer i,j,k,icount,jcount,kcount,isd,inl,irng,ie,ieface,inface

  integer nn_loc,base,top,loc_index
  integer lower, local_nn(ncpus),local_nn_temp(ncpus)
  integer flag_ca,flag,int_temp
  integer con_ele_index

  integer nn_inter_regen
  real(8) x_inter_regen(nsd,maxmatrix)
  integer index_regen_nn(ncpus),index_regen_nn_temp(ncpus)
  real(8) x_regen_temp(nsd,maxmatrix)

  if(myid==0)write(*,*)'begin finding points on contact surface for 3D'

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

  kcount=0
  nn_regen_proc=0
  x_regen_proc(:,:)=0.0
  do loc_index=1,nn_loc
     con_ele_index=myid+1+(loc_index-1)*ncpus
     ie=con_ele(myid+1+(loc_index-1)*ncpus)
     
     do ieface=1,neface
        irng=rngface(ieface,ie)
        if(irng==6)  then  ! right now slip bc is on face 6, needs to be correct later
	   do inface=1,nnface
		inl=mapping(ieface,inface,etype)
		node(inface)=ien(inl,ie)
		xloc(:,inface)=x(:,node(inface))
		I_node(inface)=I_fluid(node(inface))
	   end do
	   flag=1
	   int_temp=0
	   do inface=1,nnface
		if(I_node(inface).gt.0.55) int_temp=int_temp+1
	   end do
	   if(int_temp==4) flag=0
	   int_temp=0
           do inface=1,nnface
                if(I_node(inface).lt.0.45) int_temp=int_temp+1
           end do
           if(int_temp==4) flag=0
	   if(flag==1) then
	     call regen_points_element(xloc,x_inter,x_center,hg,I_fluid_center,corr_Ip, &
					norm_con_ele(:,con_ele_index),&
					nn_regen_proc,x_regen_proc)
	   end if


	end if ! end if of irng
      end do ! end of loop over ieface
  end do ! end of loop over nn_loc


  index_regen_nn(:)=0
  index_regen_nn_temp(:)=0
  index_regen_nn_temp(myid+1)=nn_regen_proc

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(index_regen_nn_temp(1),index_regen_nn(1),ncpus,mpi_integer,mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

  nn_inter_regen=0
  do icount=1,ncpus
     nn_inter_regen=nn_inter_regen+index_regen_nn(icount)
  end do

  lower=0
  do icount=1,myid
     lower=lower+index_regen_nn(icount)
  end do
  x_regen_temp(:,:)=0.0
  x_inter_regen(:,:)=0.0
  do icount=1,nn_regen_proc
     x_regen_temp(:,lower+icount)=x_regen_proc(:,icount)
  end do
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(x_regen_temp(1,1),x_inter_regen(1,1),nsd*maxmatrix,mpi_double_precision, &
			mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)


end subroutine points_on_contact_sur
