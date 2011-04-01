subroutine solid_node_volume(x_solid,nsd,nn_solid,ien_solid,ne,nen,snvolume,adist)
! evaluate the nodal solid volume, to get snvolume which is equalivent to dwjp in RKPM initialize
! And the support radius of every solid node
use solid_variables, only: wq_solid, nquad_solid, xq_solid
implicit none

real(8) x_solid(nsd,nn_solid) ! solid coordinates
integer nsd ! # of dimension
integer nn_solid ! # of solid
integer ien_solid(ne,nen) ! connectivity of solid mesh
integer ne ! # of solid elements
integer nen ! # of nodes per solid element
real(8) snvolume(nn_solid) ! solid nodal volume
real(8) adist(nsd,nn_solid) ! support region of solid node
real(8) vol ! for *.fi
real(8) :: xr(nsd,nsd), cf(nsd,nsd)
real(8) det
integer iq
!------------------------
real(8) xn(nsd,nen)
integer ie
integer inl
integer isd
integer nnum
real(8) toxj(nsd,nsd)
real(8) toxji(nsd,nsd)
real(8) rs(nsd)
real(8) tot_vol
real(8) coef
real(8) xmax(nsd)

coef=1.0
snvolume(:)=0.0d0
tot_vol=0.0d0
do ie=1,ne
	do inl=1,nen
		do isd=1,nsd
			nnum=ien_solid(ie,inl)
			xn(isd,inl)=x_solid(isd,nnum)
		end do
	end do
	!...Calulate radius 
	do isd=1,nsd
                xmax(isd) = coef*(maxval(xn(isd,1:nen)) - minval(xn(isd,1:nen)))
        enddo
	do inl = 1,nen
        nnum = ien_solid(ie,inl) 
        adist(1:nsd,nnum) = max(adist(1:nsd,nnum),xmax(1:nsd))
   	enddo


	 !...Calculate volume
	vol=0.0
	do iq=1,nquad_solid
		rs(1:nsd)=xq_solid(1:nsd,iq)
		call r_element(rs)
		call r_jacob(xn,toxj,toxji,det)
		vol=vol+det*wq_solid(iq)
	end do
continue
	do inl=1,nen
		nnum=ien_solid(ie,inl)
		snvolume(nnum)=snvolume(nnum)+vol/nen
	end do
end do
	tot_vol=sum(snvolume)
!	write(*,*) adist(1:nsd,:)
return
end subroutine
