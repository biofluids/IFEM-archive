subroutine rkpm_nodevolume(xna,nsd,nn,ien,ne,nen,ien_local,ne_local,&
	   dwjp_r,adist_r)
! calculate the node volume of fluid mesh for RKPM weight and support region radius
use fluid_variables, only : nquad,sq,wq
use mpi_variables
implicit none
include 'mpif.h'
real(8) xna(nsd,nn)
integer nsd
integer nn
integer ien(nen,ne)
integer nen
integer ne
integer ien_local(ne_local)
integer ne_local
!--------------------------------------------------
real(8) coef
integer ie_local
integer ie
integer inl
integer isd
real(8) xn(nsd,nen)
integer node
real(8) xmax(nsd)
integer iq
real(8) xr(nsd,nsd)
real(8) cf(nsd,nsd)
real(8) det
real(8) vol
integer nnum
!--------------------------------------------------
real(8) adist(nsd,nn)
real(8) dwjp(nn)
real(8) dwjp_r(nn)
real(8) adist_r(nsd,nn)

dwjp(:)=0.0d0
adist(:,:)=0.0d0
dwjp_r(:)=0.0d0
adist_r(:,:)=0.0d0


coef=0.8d0
if (myid == 0) write(*,*) '*** RKPM coefficient ***', coef
do ie_local=1,ne_local
	ie=ien_local(ie_local)
	do inl=1,nen
		do isd=1,nsd
			nnum=ien(inl,ie)
			xn(isd,inl)=xna(isd,nnum)
		end do
	end do

	do isd=1,nsd
		xmax(isd)=coef*(maxval(xn(isd,1:nen)) - minval(xn(isd,1:nen)))
        enddo

     do inl = 1,nen
        node = ien(inl,ie) 
        adist(1:nsd,node) = max(adist(1:nsd,node),xmax(1:nsd))
     enddo

 !...Calculate volume
        if (nsd == 2) then
          if (nen == 3) then
                        include "vol2d3n.fi"
          elseif (nen == 4) then
                        include "vol2d4n.fi"
          endif
        elseif (nsd == 3) then
      if (nen == 4) then
                        include "vol3d4n.fi"
      elseif (nen == 8) then
                        include "vol3d8n.fi"
          endif
    endif

     do inl = 1,nen
        nnum = ien(inl,ie)
        dwjp(nnum) = dwjp(nnum) + vol/nen
     enddo
  enddo

call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(dwjp(1),dwjp_r(1),nn,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call mpi_allreduce(adist(1,1),adist_r(1,1),nn*nsd,mpi_double_precision,mpi_max,mpi_comm_world,ierror)

!if (myid == 0) then
!open(unit=8406, file='dvolume_pa.txt', status='unknown')
!do isd=1,nn
!write(8406,*) 'dv',dwjp_r(isd),'radius',adist_r(1:nsd,isd)
!end do
!close(8406)
!end if

return

end

