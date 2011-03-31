subroutine search_inf_inter(xyz_solid, xyz_fluid, nn_fluids,nn_solids, nsd, ne, nen, ien, infdomain)
use interface_variables, only:maxmatrix
use mpi_variables
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
if(myid==0) then
write(*,*) 'I am in search influence inter'
endif
maxconn=30
   
   do inn=1,nn_solids
      infdomain(inn)=0
      finf=0
      x(1:nsd)=xyz_solid(1:nsd,inn)
  !    write(*,*) 'x= ',x(:)
      call getinf_el_3d(finf, x, xyz_fluid, nn_fluids, nsd, ne, nen, ien, maxconn)
      if (finf .eq. 0) then
      write(*,*) 'errors! 1. solid is out of the fluid domain or 2. search code wrong'
     stop
      end if
     ! write(*,*) 'finf=  ',finf
      infdomain(inn)=finf
   end do
      return
end
