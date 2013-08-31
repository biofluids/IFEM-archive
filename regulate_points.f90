

subroutine regulate_points(x_inter_regen,x,nn_inter_regen,ien,hg,ne_intlocal,ien_intlocal)

  use fluid_variables, only:nsd,ne,nen,nn
  use interface_variables
  use allocate_variables, only:ne_inter,inter_ele
  use mpi_variables
  include 'mpif.h'

  real(8) x_inter_regen(nsd,maxmatrix),x(nsd,nn)
  integer nn_inter_regen,ien(nen,ne)
  integer ne_intlocal,ien_intlocal(ne_intlocal)
  real(8) hg(ne)

  integer infdomain_regen(maxmatrix)
  integer i,j,icount,jcount,inl,node,ie,isd,ncount,mcount

  integer nn_loc,base,top,loc_index

  real(8) x_ele(nsd,200),tol,x_regen_temp(nsd,maxmatrix)
  integer nn_ele,index_nn(ncpus),index_nn_temp(ncpus)
  integer nn_local,nn_local_new,lower
  real(8) x_local(nsd,nn_inter_regen)

  infdomain_regen(:)=0
  x_local(:,:)=0.0
  
  call search_inf_pa_inter(x_inter_regen,x,nn,nn_inter_regen,nsd,ne,nen,ien,infdomain_regen,&
			ne_intlocal,ien_intlocal)

  call find_inter(infdomain_regen,ien,nn_inter_regen)
  if(ne_inter.le.ncpus) then
    if(myid+1.le.ne_inter) then
      nn_loc=1
    else
      nn_loc=0
    end if
  else
    base=floor(real(ne_inter)/real(ncpus))
    top=ne_inter-base*ncpus
    if(myid+1.le.top) then
      nn_loc=base+1
    else
      nn_loc=base
    end if
  end if

  nn_local=0
  nn_local_new=0
  do loc_index=1,nn_loc
     i=myid+1+(loc_index-1)*ncpus
!     tol=hg(inter_ele(i))/8.0    
     tol=max_hg/5.0
     nn_ele=0
     do icount=1,nn_inter_regen
        if(inter_ele(i)==infdomain_regen(icount)) then
          nn_ele=nn_ele+1
if(nn_ele.gt.200)write(*,*)'exceed matrix limit in regulate_points.f90'
          x_ele(:,nn_ele)=x_inter_regen(:,icount)
        end if
     end do
     call points_deletion(x_ele,nn_ele,nsd,tol)
     nn_local=nn_local+nn_ele
     x_local(:,nn_local-nn_ele+1:nn_local)=x_ele(:,1:nn_ele)
  end do
  x_inter_regen(:,:)=0.0
  nn_inter_regen=0
  index_nn(:)=0
  index_nn_temp(:)=0
  index_nn_temp(myid+1)=nn_local

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(index_nn_temp(1),index_nn(1),ncpus,mpi_integer,mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)
  do icount=1,ncpus
     nn_inter_regen=nn_inter_regen+index_nn(icount)
  end do
  
  lower=0
  do icount=1,myid
     lower=lower+index_nn(icount)
  end do
  x_regen_temp(:,:)=0.0
  do icount=1,nn_local
     x_regen_temp(:,lower+icount)=x_local(:,icount)
  end do
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(x_regen_temp(1,1),x_inter_regen(1,1),nsd*maxmatrix,mpi_double_precision, &
			mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

end subroutine regulate_points



     
  
   
subroutine points_deletion(x_ele,nn_ele,nsd,tol)

  integer nn_ele,nn_ele_new,nsd
  real(8) x_ele(nsd,50),x_ele_new(nsd,50),tol

  integer i,j,icount,jcount,isd
  real(8) temp
  real(8) x1(nsd),x2(nsd)
  nn_ele_new=0
  do i=1,nn_ele
     x1(:)=x_ele(:,i)
     do j=1,nn_ele_new
        x2(:)=x_ele_new(:,j)
        temp=0.0
	do isd=1,nsd
	   temp=temp+(x1(isd)-x2(isd))**2
	end do
	temp=sqrt(temp)
	if(temp.le.tol) then
	   goto 234
	end if
     end do
     nn_ele_new=nn_ele_new+1
     x_ele_new(:,nn_ele_new)=x1(:)
234 continue
  end do
  nn_ele=nn_ele_new
  x_ele=x_ele_new
end subroutine points_deletion

     


     





