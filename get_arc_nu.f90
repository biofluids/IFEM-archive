!get arclength for free surface type of cases

subroutine get_arc_nu(arc_inter,norm_inter,x_inter)

  use interface_variables
  use fluid_variables, only:nsd
  use mpi_variables
  include 'mpif.h'
  real(8) arc_inter(maxmatrix),arc_inter_temp(maxmatrix)
  real(8) norm_inter(nsd,maxmatrix)
  real(8) x_inter(nsd,maxmatrix)

  integer i,j,isd,jsd
  real(8) Bsum
  real(8) Sp,temp,dx(nsd)

  real(8) range_p, range_n
  real(8) dis,bound,support

!---used for mpi
  integer nn_inter_loc,base,top,loc_index

  if(nn_inter.le.ncpus) then
    if(myid+1.le.nn_inter) then
      nn_inter_loc=1
    else 
      nn_inter_loc=0
    end if 
  else
     base=floor(real(nn_inter)/real(ncpus))
     top=nn_inter-base*ncpus
     if(myid+1.le.top) then
        nn_inter_loc=base+1
     else 
        nn_inter_loc=base
     end if
  end if


  arc_inter(:)=0.0
  arc_inter_temp(:)=0.0

  total_length=0.0

!  do i=1,nn_inter
  do loc_index=1,nn_inter_loc
     i=myid+1+(loc_index-1)*ncpus
     range_p=0.0
     range_n=0.0
     Bsum=0.0
     do j=1,nn_inter
	dx(:)=x_inter(:,j)-x_inter(:,i)

	dis=0.0
	do isd=1,nsd
	   dis=dis+dx(isd)**2
	end do
	dis=sqrt(dis)
	call B_Spline_0order(dis,hsp,Sp)
	Bsum=Bsum+Sp

	temp=dx(1)*norm_inter(2,i)-dx(2)*norm_inter(1,i)
	if(temp.gt.0) then
	   if(dis.gt.range_p)then
	     range_p=dis
	   end if
	else
	   if(dis.gt.range_n)then
	     range_n=dis
	   endif
	end if
    end do
    bound=min(range_p,range_n)/hsp
    support=0.0
    if(bound.le.1) then
	support=1.0/120.0*(-1.0/6.0*(3.0-bound)**6+(2-bound)**6-15.0/6.0*(1-bound)**6)+0.5
    else if(bound.le.2) then
	support=1.0/120.0*(-1.0/6.0*(3.0-bound)**6+(2-bound)**6)+0.5
    else if(bound.le.3) then
	support=1.0/120.0*(-1.0/6.0*(3.0-bound)**6)+0.5
    else
	support=0.5
    end if

    support=support+0.5
    arc_inter_temp(i)=support*hsp/Bsum
!    total_length=total_length+arc_inter(i)
  end do

  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_reduce(arc_inter_temp(1),arc_inter(1),maxmatrix,mpi_double_precision, &
                mpi_sum,0,mpi_comm_world,ierror)
  call mpi_bcast(arc_inter(1),maxmatrix,mpi_double_precision,0,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)

  do i=1,nn_inter
     total_length=total_length+arc_inter(i)
  end do

 if(myid==0)write(*,*)'total_length=',total_length
end subroutine get_arc_nu


 
