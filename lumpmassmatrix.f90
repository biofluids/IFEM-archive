! Jack Yang, Rensselaer Polytechnic Institute
! June 2014
! this subroutine evaluates the special lumped mass matrix
! see Hughes, p445
subroutine lumpmassmatrix(xloc, dloc, doloc, p, hk, ien, f_fluids,ne_local,ien_local,node_local,nn_local, &
                         fden,fvis,I_fluid,rngface)
  use global_constants
  use run_variables
  use fluid_variables
  use solid_variables, only: nn_solid
  use r_common, only: density_solid, vis_solid
  use lumpedmass
  implicit none

  integer ien(nen,ne)
  real* 8 xloc(nsd,nn)                      ! geometry(nodes' locations)
  real* 8 dloc(ndf,nn),doloc(ndf,nn)        ! unknowns(current/previous)
  real* 8 p(ndf,nn),hk(ne)                  ! p: residual; hk: element characteristic length

  real* 8 x(nsd,nen)
  real* 8 d(ndf,nen),d_old(ndf,nen)

  real* 8 eft0,det,effd,effm,effc
  real* 8 sh(0:nsd,nen),ph(0:nsd,nen)
  real* 8 xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real* 8 drt(ndf),drs(ndf),drtavg(ndf),drsavg(ndf)
  real* 8 dr(nsd,ndf),dravg(nsd,ndf)
  real* 8 u,v,w,pp,ug
  real* 8 tau(nsd,nsd)
  real* 8 hg,taum,tauc,vel,ree, taul
  real* 8 mu,nu,ro,g(nsd)
  real* 8 tempc(ndf),temp
  real* 8 dtinv,oma,ama
  integer inl, ie, isd, iq, node,jsd
  integer ieface,irng !,inface

  real* 8 f_fluids(nsd,nn)
  real* 8 fnode(nsd,nen),fq(nsd)
  integer nn_local
  integer node_local(nn_local)
!-------------------------------------------
  real(8) fden(nn)
  real(8) local_den(nen)
  real(8) fvis(nn)
  real(8) local_vis(nen)
!---------------------------------------------
  real(8) TC,ZC,RC,P0,dens0
  real(8) I_fluid(nn)
  real(8) kappa
  integer rngface(neface,ne)
!============================
! MPI varibalbes
  integer ne_local ! # of element on each processor
  integer ien_local(ne_local) ! subregion-->wholeregion element index
  integer ie_local ! loop parameter
  integer icount
!-----------------------------------------
! special lumped mass matrix variables
  real(8) alphalump
  real(8) vol, sumdiagcm, inddiagcm(nen)
  integer error_id1, error_id2, error_id3
!===================================================

if (allocated(mpqlumpmass)) then
    deallocate(mpqlumpmass)
endif
if (allocated(hattauij)) then
    deallocate(hattauij)
endif
if (allocated(flagdivision)) then
    deallocate(flagdivision)
endif
allocate(mpqlumpmass(nn)     ,stat=error_id1)
allocate(hattauij(nsd,nsd,nn),stat=error_id2)
allocate(flagdivision(nn)    ,stat=error_id3)

mpqlumpmass(1:nn)=0.0

do ie=1,ne     ! loop over elements
    do inl=1,nen
        node=ien(inl,ie)
        x(1:nsd,inl) = xloc(1:nsd,ien(inl,ie))
    enddo

    hg = hk(ie)
    vol=0.0
    sumdiagcm=0.0
    inddiagcm(1:nen)=0.0
    do iq=1,nquad  ! loop over the quadrature points in each element 
!...  calculate the shape function and the weight at quad point
        if (nsd==2) then
            if (nen.eq.3) then !calculate shape function at quad point
                include "sh2d3n.h"
            elseif (nen.eq.4) then
                include "sh2d4n.h"
            endif
        elseif (nsd==3) then
            if (nen.eq.4) then !calculate shape function at quad point
                include "sh3d4n.h"
            elseif (nen.eq.8) then
                include "sh3d8n.h"
            endif
        endif
!	write(*,*) 'shape functions=',sh(0,1),sh(0,2),sh(0,3),sh(0,4)
        eft0 = abs(det) * wq(iq) ! calculate the weight at each quad pt
        vol=vol+eft0
        ph(0:nsd,1:nen) = sh(0:nsd,1:nen)*eft0
        do inl=1,nen
            inddiagcm(inl)=inddiagcm(inl)+eft0*sh(0,inl)**2.0
            sumdiagcm=sumdiagcm+eft0*sh(0,inl)**2.0
        enddo
    enddo ! end of qudrature pts loop

    alphalump=vol/sumdiagcm

    do inl=1,nen
        node=ien(inl,ie)
        mpqlumpmass(node)=mpqlumpmass(node)+alphalump*inddiagcm(inl)
    enddo
enddo ! end of element loop

return
end subroutine lumpmassmatrix

