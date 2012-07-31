!=====================================
!find center & den calculation domain pa
!=====================================

subroutine find_domain_pa_V2(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local)

  use interface_variables, only:nn_inter,maxmatrix,hsp,max_hg
  use fluid_variables,only:nsd,ne,nn
!  use denmesh_variables, only:nn_den,ne_den,nen_den
  use allocate_variables, only:den_domain,center_domain,ne_den_domain,nn_center_domain
  use mpi_variables
  include 'mpif.h'

  real(8) x_center(nsd,ne),x(nsd,nn)
!  real(8) x_den(nsd,nn_den)
  real(8) x_inter(nsd,maxmatrix)
  integer ne_intlocal
  integer ien_intlocal(ne_intlocal)
  integer nn_domain_local
!  integer domain_local(ne_intlocal)
  integer domain_local(ne_intlocal)
  integer domain_nlocal(nn_local)
   real(8) hg(ne)
  integer nn_local,node_local(nn_local)
!  integer ne_local_den
!  integer ien_local_den(ne_local_den)
!  integer ien_den(nen_den,ne_den)

  integer i,j,icount,jcount,ie,isd,inl
  real(8) temp,support,den_center(nsd)
  integer flag

  integer index_local_temp(ncpus),index_local(ncpus)
  integer lower
  integer,dimension(:),allocatable :: domain_temp
!  support=4.0*maxval(hg(:))

!find center domain as the interpolation region
  support=20.0*hsp  !4*hg is enough for the regenration part
  nn_domain_local=0
  index_local_temp(:)=0
  index_local(:)=0
  do icount=1,ne_intlocal
     ie=ien_intlocal(icount)
     flag=0
     do i=1,nn_inter
        temp=0.0
        do isd=1,nsd
	   temp=temp+(x_center(isd,ie)-x_inter(isd,i))**2
	end do
	if(sqrt(temp).lt.support) then
	  flag=1
	  goto 100
	end if
     end do
100 continue
    if(flag==1) then
	nn_domain_local=nn_domain_local+1
	domain_local(nn_domain_local)=ie
    end if
  end do

  index_local_temp(myid+1)=nn_domain_local
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(index_local_temp(1),index_local(1),ncpus,mpi_integer,mpi_sum,&
		mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)
  nn_center_domain=0
  do icount=1,ncpus
     nn_center_domain=nn_center_domain+index_local(icount)
  end do
  if(allocated(center_domain)) then
    deallocate(center_domain)
  end if
  allocate(center_domain(nn_center_domain))

  if(allocated(domain_temp)) then
    deallocate(domain_temp)
  end if
  allocate(domain_temp(nn_center_domain))
  domain_temp(:)=0
  lower=0
  do icount=1,myid
     lower=lower+index_local(icount)
  end do
    
  do icount=1,index_local(myid+1)
     domain_temp(lower+icount)=domain_local(icount)
  end do
  center_domain(:)=0
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(domain_temp(1),center_domain(1),nn_center_domain,mpi_integer,mpi_sum,0,&
		mpi_comm_world,ierror)

  call mpi_bcast(center_domain(1),nn_center_domain,mpi_integer,0,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)
!**********************************************!

!**********************************************!
!find the narrow band of the dense mesh used for poisson equation
!  support=1.1*max_hg
  support=5.0*hsp
  nn_domain_local=0
  index_local_temp(:)=0
  index_local(:)=0
  do icount=1,ne_intlocal
     ie=ien_intlocal(icount)
!     den_center(1:nsd)=0.0
!     do isd=1,nsd
!	do inl=1,nen_den
!	   den_center(isd)=den_center(isd)+1.0/real(nen_den)*x_den(isd,ien_den(inl,ie))
!	end do
!     end do
     flag=0
     do i=1,nn_inter
        temp=0.0
        do isd=1,nsd
           temp=temp+(x_center(isd,ie)-x_inter(isd,i))**2
        end do
        if(sqrt(temp).lt.support) then
          flag=1
          goto 900
        end if
     end do
900 continue
    if(flag==1) then
        nn_domain_local=nn_domain_local+1
        domain_local(nn_domain_local)=ie
    end if
  end do

  index_local_temp(myid+1)=nn_domain_local
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(index_local_temp(1),index_local(1),ncpus,mpi_integer,mpi_sum,0,&
                mpi_comm_world,ierror)
  call mpi_bcast(index_local(1),ncpus,mpi_integer,0,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)
  ne_den_domain=0
  do icount=1,ncpus
     ne_den_domain=ne_den_domain+index_local(icount)
  end do
  if(allocated(den_domain)) then
    deallocate(den_domain)
  end if
  allocate(den_domain(ne_den_domain))

  if(allocated(domain_temp)) then
    deallocate(domain_temp)
  end if
  allocate(domain_temp(ne_den_domain))
!write(*,*)myid, 'ne_den_domain=',ne_den_domain 
 domain_temp(:)=0
  lower=0
  do icount=1,myid
     lower=lower+index_local(icount)
  end do

  do icount=1,index_local(myid+1)
     domain_temp(lower+icount)=domain_local(icount)
  end do
  den_domain(:)=0
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(domain_temp(1),den_domain(1),ne_den_domain,mpi_integer,mpi_sum,0,&
                mpi_comm_world,ierror)

  call mpi_bcast(den_domain(1),ne_den_domain,mpi_integer,0,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

!if(myid==0) then
!do icount=1,ne_den_domain
!  write(*,*)'ne_den_domain=',icount,'den_doamin=',den_domain(icount)
!end do  
!end if 
!stop

end subroutine find_domain_pa_V2
  
