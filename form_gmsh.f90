module form_gmsh
  implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formd_gmsh(d,meshvel,ien_surf,bid)
!...applies boundary values to unknown vector d
  use fluid_variables, only: nn,ne,nsd,ndf,nrng,udf,vdf,wdf,static,bv,ne_surf,nen_surf
  implicit none

  real(8),intent(inout) :: d(ndf,nn)
  real(8),intent(in)    :: meshvel(nsd,nn)
  integer,intent(in)    :: ien_surf(1:nen_surf,1:ne_surf),bid(ne_surf)

  integer :: idf, irng, inn,ine,inen,isd

  do irng = 1,nrng
     do idf = 1,ndf
        if (bv(idf,irng) > -900.0d0) then
           do ine = 1,ne_surf
              if (bid(ine) == irng) then
                 do inen = 1,nen_surf
                    d(idf,ien_surf(inen,ine)) = bv(idf,irng)
                 enddo
              endif
           enddo
        endif
     enddo
  enddo

  !d(2,:) = 0.0d0

  if(static) then
     do inn = 1,nn
        do isd = 1,nsd 
           d(isd,inn) = 0.0
        enddo 
     enddo
  endif

  return
end subroutine formd_gmsh


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formdm_gmsh(disp,x,ien_surf,bid)
  use fluid_variables, only: nn,ne_surf,nen_surf,nrng,nsd,static,nsdpad,maxnsurf,bvd
  implicit none

  real(8),intent(inout) :: disp(nsd,nn)
  real(8),intent(in)    :: x(nsd,nn)
  integer,intent(in)    :: ien_surf(nen_surf,ne_surf)
  integer,intent(in)    :: bid(ne_surf)

  integer :: irng, inn,inen,ine,isd

  do irng = 1,nrng
     do isd = 1,nsd
        if (bvd(isd,irng) > -900.0d0) then
           do ine = 1,ne_surf
              if (bid(ine) == irng) then
                 do inen = 1,nen_surf
                    disp(isd,ien_surf(inen,ine)) = bvd(isd,irng)
                 enddo
              endif
           enddo
        endif
     enddo
  enddo

  if(static) then
     do inn = 1,nn
        do isd = 1,nsd 
           disp(isd,inn) = 0.0d0
        enddo 
     enddo
  endif

  return
end subroutine formdm_gmsh


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formid_gmsh(id,ien_surf,bid)
  use fluid_variables, only: nn,ndf,nrng,bv,nq,nsd,ndfpad,maxnsurf,static,ne_surf,nen_surf
  implicit none

  integer,intent(out) :: id(ndf,nn)
  integer,intent(in)  :: ien_surf(1:nen_surf,1:ne_surf),bid(ne_surf)

  integer :: idf,isd,irng,inn,ine,inen

  id(:,:) = 0
  do irng = 1,nrng
     do idf = 1,ndf
        if (bv(idf,irng) > -900.0d0) then
           do ine = 1,ne_surf
              if (bid(ine) == irng) then
                 do inen = 1,nen_surf
                    id(idf,ien_surf(inen,ine)) = 1
                 enddo
              endif
           enddo
        endif
     enddo
  enddo

  !id(2,:) = 1

  if(static) then
     do inn = 1,nn
        do isd = 1,nsd 
           id(isd,inn) = 1
        enddo 
     enddo
  endif

  nq = 0
  do inn=1,nn
     do idf=1,ndf
        if(id(idf,inn) == 0) then
           nq = nq + 1
           id(idf,inn) = nq
        else
           id(idf,inn) = 0
        endif
     enddo
  enddo

  if (nq == 0) then 
     nq = 1
     write(*,*) " all nodes are constrained !!! -> nq set to: ",nq
  endif

  return
end subroutine formid_gmsh


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formidm_gmsh(kid,ien_surf,bid)
  use fluid_variables, only: nn,nrng,nq,nsd,ndfpad,maxnsurf,static,ne_surf,nen_surf,bvd
  implicit none

  integer,intent(out) :: kid(nsd,nn)
  integer,intent(in)  :: ien_surf(1:nen_surf,1:ne_surf),bid(ne_surf)

  integer :: isd,irng,inn,ine,inen

  kid(:,:) = 0
  do irng = 1,nrng
     do isd = 1,nsd
        if (bvd(isd,irng) > -900.0d0) then
           do ine = 1,ne_surf
              if (bid(ine) == irng) then
                 do inen = 1,nen_surf
                    kid(isd,ien_surf(inen,ine)) = 1
                 enddo
              endif
           enddo
        endif
     enddo
  enddo

  if(static) then
     do inn = 1,nn
        do isd = 1,nsd 
           kid(isd,inn) = 1
        enddo 
     enddo
  endif

  nq = 0
  do inn=1,nn
     do isd=1,nsd
        if(kid(isd,inn) == 0) then
           nq = nq + 1
           kid(isd,inn) = nq
        else
           kid(isd,inn) = 0
        endif
     enddo
  enddo

  if (nq == 0) then
     nq = 1 !...if degree of freedom is 0, than nq = 1, so that gmres can calculate the residual  (division by zero error)
     write(*,*) " all mesh nodes are constrained !!! -> nq set to: ",nq
  endif

  return
end subroutine formidm_gmsh



end module form_gmsh


