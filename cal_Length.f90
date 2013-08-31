

subroutine cal_length(x_inter,dmass,x,I_fluid,ien,Length)

  use fluid_variables
  use interface_variables
  use mpi_variables
  use allocate_variables,only:ne_inter,inter_ele
  include 'mpif.h'

  real(8) x_inter(nsd,maxmatrix),dmass,x(nsd,nn)
  real(8) I_fluid(nn)
  integer ien(nen,ne)

  integer nn_loc,base,top,loc_index
  integer i,j,icount,jcount,inl,ie,node,isd,node1,node2,node3,node4
  real(8) II(4)
  integer nloc
  real(8) xloc(nsd,4),sh(4),x_edge(nsd,2)
  real(8) Length0, Length

  
  Length0=0.0
  Length=0.0
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

  do loc_index=1,nn_loc
!   do loc_index=1,1
     icount=myid+1+(loc_index-1)*ncpus
     ie=inter_ele(icount)
!     ie=icount
     node1=ien(1,ie)
     node2=ien(2,ie)
     node3=ien(3,ie)
     node4=ien(4,ie)
     II(1)=I_fluid(node1)
     II(2)=I_fluid(node2)
     II(3)=I_fluid(node3)
     II(4)=I_fluid(node4)
     nloc=0
     if((II(1)-0.5)*(II(2)-0.5).lt.0.0) then
       nloc=nloc+1
       xloc(1,nloc)=(II(1)+II(2)-1.0)/(II(1)-II(2))
       if(abs(xloc(1,nloc)).gt.1.0)then
          write(*,*)'calculate edge wrong'
	  stop
	end if
       xloc(2,nloc)=-1.0
     end if

     if((II(2)-0.5)*(II(3)-0.5).lt.0.0) then
       nloc=nloc+1
       xloc(2,nloc)=(II(2)+II(3)-1.0)/(II(2)-II(3))
       if(abs(xloc(2,nloc)).gt.1.0)then
         write(*,*)'calculate edge wrong'
         stop
       end if
       xloc(1,nloc)=1.0
     end if

     if((II(3)-0.5)*(II(4)-0.5).lt.0.0) then
       nloc=nloc+1
       xloc(1,nloc)=(II(4)+II(3)-1.0)/(II(4)-II(3))
       if(abs(xloc(1,nloc)).gt.1.0) then
         write(*,*)'calcualte edeg eworong '
	 stop
	end if
       xloc(2,nloc)=1.0
     end if

     if((II(4)-0.5)*(II(1)-0.5).lt.0.0) then
       nloc=nloc+1
       xloc(2,nloc)=(II(1)+II(4)-1.0)/(II(1)-II(4))
       if(abs(xloc(2,nloc)).gt.1.0) then
         write(*,*)'wrong'
	 stop
       end if
       xloc(1,nloc)=-1.0
     end if
!if(myid==0) then
!  write(*,*)'xloc1=',xloc(1:2,1)
!  write(*,*)'xloc2=',xloc(1:2,2)
!end if
!     if(nloc.ne.2) then
!       write(*,*)'exceedd points in edge'
!       write(*,*)'ie=',ie
!       write(*,*)'II=',II(1:4)
!       stop
!     end if
     if(nloc==2) then

     x_edge(:,:)=0.0
     do i=1,nloc
        sh(1)=0.25*(1-xloc(1,i))*(1-xloc(2,i))
	sh(2)=0.25*(1+xloc(1,i))*(1-xloc(2,i))
	sh(3)=0.25*(1+xloc(1,i))*(1+xloc(2,i))
	sh(4)=0.25*(1-xloc(1,i))*(1+xloc(2,i))
	do inl=1,nen
           x_edge(:,i)=x_edge(:,i)+sh(inl)*x(:,ien(inl,ie))
	end do
     end do
!     if(myid==0) then
!     write(*,*)'x1=',x(:,ien(1,ie))
!     write(*,*)'x2=',x(:,ien(2,ie))
!     write(*,*)'x3=',x(:,ien(3,ie))
!     write(*,*)'x4=',x(:,ien(4,ie))
!     write(*,*)'x_edge1=',x_edge(:,1)
!     write(*,*)'x_edge2=',x_edge(:,2)
!     end if
     Length0=Length0+sqrt((x_edge(1,1)-x_edge(1,2))**2+(x_edge(2,1)-x_edge(2,2))**2)
     end if
    if(nloc.ne.2) then
!       write(*,*)'ie=',ie,'I=',II(1:4)
     Length0=Length0+max_hg
    end if


  end do
!write(*,*)'myid=',myid,'Length0=',Length0
     call mpi_barrier(mpi_comm_world,ierror)
     call mpi_allreduce(Length0,Length,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
     call mpi_barrier(mpi_comm_world,ierror)

     if(myid==0)write(*,*)'lenght=',Length

end subroutine cal_Length


