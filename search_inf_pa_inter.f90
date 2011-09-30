subroutine search_inf_pa_inter(xyz_solid, xyz_fluid, nn_fluids,nn_solids, nsd, ne, nen, ien, infdomain, & 
			ne_intlocal, ien_intlocal )
  use mpi_variables ! call mpi variable module
  use interface_variables, only:maxmatrix
  include 'mpif.h'
integer nn_fluids
integer nn_solids
integer nsd
integer ne
integer nen
integer ien(nen,ne)
integer infdomain(maxmatrix)
real(8) xyz_solid(nsd,nn_solids)
real(8) xyz_fluid(nsd,nn_fluids)

integer finf
real(8) x(nsd)

integer inn
integer maxconn
integer nn_inter_loc,base,top,loc_index

! Additioanl variables to make search code parallel
integer ne_intlocal
integer ien_intlocal(ne_intlocal)
integer inf_tmp(nn_solids),inf_tmp2(nn_solids)
integer outflag
if (myid ==0) then
write(*,*) 'I am in search influence'
end if

  if(nn_solids.le.ncpus) then
    if(myid+1.le.nn_solids) then
      nn_inter_loc=1
    else
      nn_inter_loc=0
    end if
  else
     base=floor(real(nn_solids)/real(ncpus))
     top=nn_solids-base*ncpus
     if(myid+1.le.top) then
        nn_inter_loc=base+1
     else
        nn_inter_loc=base
     end if
  end if





maxconn=30
outflag=0
infdomain(:)=0   
inf_tmp(:)=0
inf_tmp2(:)=0
   do loc_index=1,nn_inter_loc
      inn=myid+1+(loc_index-1)*ncpus
      finf=0
      x(1:nsd)=xyz_solid(1:nsd,inn)
  !    call getinf_el_3d_pa(finf, x, xyz_fluid, nn_fluids, nsd, ne, nen, ien, maxconn, ne_intlocal, ien_intlocal)
call getinf_el_3d(finf, x, xyz_fluid, nn_fluids, nsd, ne, nen, ien, maxconn)
  !    if (finf .eq. 0) then
  !    write(*,*) 'errors! 1. solid is out of the fluid domain or 2. search code wrong'
  !   stop
  !    end if
      inf_tmp(inn)=finf
   end do

	      call mpi_barrier(mpi_comm_world,ierror)
              call mpi_reduce(inf_tmp(1),inf_tmp2(1),nn_solids,mpi_integer,mpi_sum,0,mpi_comm_world,ierror)
	      call mpi_bcast(inf_tmp2(1),nn_solids,mpi_integer,0,mpi_comm_world,ierror)
	      infdomain(1:nn_solids)=inf_tmp2(1:nn_solids)
	if (myid == 0) then
		outflag=minval(infdomain(1:nn_solids))
		if (outflag == 0) then
			 write(*,*) 'errors! 1. solid is out of the fluid domain or 2. search code wrong'
			 stop
		end if
	end if

      return
end
