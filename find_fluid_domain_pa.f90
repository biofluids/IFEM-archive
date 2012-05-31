!=====================================
!find center & den calculation domain pa
!=====================================

subroutine find_fluid_domain_pa(x,x_center,x_inter,ne_intlocal,ien_intlocal,hg,nn_local,node_local)

  use interface_variables, only:nn_inter,maxmatrix,hsp,max_hg
  use fluid_variables,only:nsd,ne,nn
!  use denmesh_variables, only:nn_den,ne_den,nen_den
  use allocate_variables, only:fluid_domain,nn_fluid_domain
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

!**********************************************!
  support=4.0*hsp  !4*hg is enough for the regenration part
  nn_domain_local=0
  index_local_temp(:)=0
  index_local(:)=0
  domain_nlocal(:)=0
  do icount=1,nn_local
     ie=node_local(icount)
     flag=0
     do i=1,nn_inter
        temp=0.0
        do isd=1,nsd
           temp=temp+(x(isd,ie)-x_inter(isd,i))**2
        end do
        if(sqrt(temp).lt.support) then
          flag=1
          goto 300
        end if
     end do
300 continue
    if(flag==1) then
        nn_domain_local=nn_domain_local+1
        domain_nlocal(nn_domain_local)=ie
    end if
  end do
  index_local_temp(myid+1)=nn_domain_local
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(index_local_temp(1),index_local(1),ncpus,mpi_integer,mpi_sum,&
                mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)
  nn_fluid_domain=0
  do icount=1,ncpus
     nn_fluid_domain=nn_fluid_domain+index_local(icount)
  end do
  if(allocated(fluid_domain)) then
    deallocate(fluid_domain)
  end if
  allocate(fluid_domain(nn_fluid_domain))

  if(allocated(domain_temp)) then
    deallocate(domain_temp)
  end if
  allocate(domain_temp(nn_fluid_domain))
  domain_temp(:)=0
  lower=0
  do icount=1,myid
     lower=lower+index_local(icount)
  end do

  do icount=1,index_local(myid+1)
     domain_temp(lower+icount)=domain_nlocal(icount)
  end do
  fluid_domain(:)=0
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(domain_temp(1),fluid_domain(1),nn_fluid_domain,mpi_integer,mpi_sum,&
                mpi_comm_world,ierror)

  call mpi_barrier(mpi_comm_world,ierror)

!**********************************************!
!find the narrow band of the dense mesh used for poisson equation

end subroutine find_fluid_domain_pa
  
