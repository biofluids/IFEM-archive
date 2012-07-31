

subroutine set_Ifluidcenter(I_fluid_center,I_fluid,ien,x_inter,x_center,hg,corr_Ip,flag)

  use mpi_variables
  use fluid_variables, only:ne,nen,nn,nsd
  use allocate_variables, only:den_domain, ne_den_domain
  use interface_variables
  include 'mpif.h'
  real(8) I_fluid_center(ne),I_fluid_temp(ne)
  real(8) I_fluid(nn)
  integer ien(nen,ne)
  real(8) x_inter(nsd,maxmatrix),x_center(nsd,ne),hg(ne),corr_Ip(maxmatrix)
  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) norm_a(nsd),curv_a
  real(8) coef
  integer ie,je,icount,inl,node,flag

  integer ne_loc,base,top,loc_index

  if(ne.le.ncpus) then
    if(myid+1.le.ne) then
      ne_loc=1
    else
      ne_loc=0
    end if
  else
     base=floor(real(ne)/real(ncpus))
     top=ne-base*ncpus
     if(myid+1.le.top) then
        ne_loc=base+1
     else
        ne_loc=base
     end if
  end if

  I_fluid_temp(:)=0.0
!  I_fluid_center(:)=0.0
if(flag==1)coef=1.01
if(flag==2) coef=1.01

!  do ie=1,ne
  do loc_index=1,ne_loc
     ie=myid+1+(loc_index-1)*ncpus
     I_fluid_temp(ie)=0.0


!     call get_indicator_derivative_2D_1st(x_center(:,ie),x_inter,x_center,hg,I_fluid_center,corr_Ip,&
!			II,dI,ddI,norm_a,curv_a)
!     I_fluid_temp(ie)=II
     do inl=1,nen
        I_fluid_temp(ie)=I_fluid_temp(ie)+1.0/real(nen)*I_fluid(ien(inl,ie))
     end do
     I_fluid_temp(ie)=(I_fluid_temp(ie)-0.5)*coef
     if(I_fluid_temp(ie).gt.0.499999) I_fluid_temp(ie)=0.5
     if(I_fluid_temp(ie).lt.-0.499999) I_fluid_temp(ie)=-0.5
     I_fluid_temp(ie)=I_fluid_temp(ie)+0.5
 



!     if(I_fluid_temp(ie).ge.0.99) then
!	I_fluid_temp(ie)=1.0
!     elseif(I_fluid_temp(ie).le.0.01) then
!	I_fluid_temp(ie)=0.0
!     end if
  end do
I_fluid_center(:)=0.0
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(I_fluid_temp(1),I_fluid_center(1),ne,mpi_double_precision, &
		mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

goto 999
  do icount=1,ne_den_domain
     ie=den_domain(icount)
!  do icount=1,ne_regen_ele
!     ie=regen_ele(icount)
     I_fluid_center(ie)=0.0
     do inl=1,nen
	I_fluid_center(ie)=I_fluid_center(ie)+1.0/real(nen)*I_fluid(ien(inl,ie))
     end do
  end do
999 continue
end subroutine set_Ifluidcenter
