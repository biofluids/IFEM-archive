

subroutine regulate_points(x_inter_regen,x,nn_inter_regen,ien,ne_intlocal,ien_intlocal)

  use fluid_variables, only:nsd,ne,nen,nn
  use interface_variables
  use allocate_variables, only:ne_inter
  use mpi_variables
  include 'mpif.h'

  real(8) x_inter_regen(nsd,maxmatrix),x(nsd,nn)
  integer nn_inter_regen, ien(nen,ne)
  integer ne_intlocal, ien_intlocal(ne_intlocal)

  integer infdomain_regen(maxmatrix)
  integer i,j,icount,jcount,inl,node,ie,isd,ncount,mcount

  integer index_ele(ne_inter*2),index_point(ne_inter*2,50),nn_point(50)

  integer nn_loc,base,top,loc_index

  integer nn_local,nn_local_temp,nn_ele
  infdomain_regen(:)=0
  index_ele(:)=0
  index_point(:,:)=0
  nn_point(:)=0
  call search_inf_pa_inter(x_inter_regen,x,nn,nn_inter_regen,nsd,ne,nen,ien,infdomain_regen, &
				ne_intlocal,ien_intlocal)
  ncount=0
  do i=1,nn_inter_regen
     do j=1,ncount
        if(infdomain_regen(i)==index_ele(j)) then
	  nn_point(j)=nn_point(j)+1
	  index_point(j,nn_point(j))=i
	  goto 100
	end if
     end do
     ncount=ncount+1
     index_ele(ncount)=infdomain_regen(i)
     nn_point(ncount)=nn_point(ncount)+1
     index_point(ncount,nn_point(ncount))=i
100 continue
  end do  

  if(ncount.le.ncpus) then
    if(myid+1.le.ncount) then
      nn_loc=1
    else
      nn_loc=0
    end if
  else
     base=floor(real(ncount)/real(ncpus))
     top=ncount-base*ncpus
     if(myid+1.le.top) then
        nn_loc=base+1
     else
        nn_loc=base
     end if
  end if

  nn_local_temp=0
  nn_local=0

  do loc_index=1,nn_loc
     i=myid+1+(loc_index-1)*ncpus
     nn_ele=0
     do icount=1,nn_point(i)
	node1=index_point(i,icount)
        do jcount=1,nn_ele
           temp=0.0
	   do isd=1,nsd
	      temp=temp+(x_inter_regen(isd,node1)-x_inter_loc(isd,jcount))**2
	   end do
	   temp=sqrt(temp)  







  
