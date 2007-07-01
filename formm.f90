module formm
  implicit none


contains

!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine formdm(x,ds,rngface,ien)
  use run_variables, only: tt
  use fluid_variables
  implicit none

  integer :: rngface(neface,ne),ien(nen,ne)
  real(8) :: x(nsd,nn)
  real(8) :: ds(nsd,nn)
  integer :: isd, inl, iec, irng, ieface, inface, inn

  real(8) :: eps1,eps2
  integer :: h(nrng,nn)


!  construct lien and mien mappings from ien array

  eps1 = -1000000.0d0
  eps2 = -10000.0d0

  h(1:nrng,1:nn) = 0

  !ds(1:nsd,1:nn) = eps1   !first guess for iteration
  ds(1:nsd,1:nn) = 0.0d0

  do ieface=1,neface
     do inface=1,nnface
        inl = mapping(ieface,inface,etype)
        do iec=1,ne
           irng = rngface(ieface,iec)
           if(irng.ne.0) h(irng,ien(inl,iec)) = 1  
        enddo
     enddo
  enddo

  do irng=1,nrng
     do inn=1,nn
        do isd=1,nsd
           if (h(irng,inn) == 1) then
              if (bcd(isd,irng) > 0) ds(isd,inn) = bvd(isd,irng)
           endif
        enddo
     enddo
  enddo

  do inn = 1,nn
     do isd = 1,nsd 
        if(ds(isd,inn) <= eps1) ds(isd,inn) = 0.0
     enddo 
  enddo

  return
end subroutine formdm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formidm(ids, rngface, ien)
  use fluid_variables
  implicit none

  integer :: ids(nsd,nn),rngface(neface,ne),ien(nen,ne)
  integer :: isd, inl, iec, irng, ieface, inface, inn

  real(8) :: epsr,epsl

  ids(1:nsd,1:nn) = 0

  epsr = 0.0001
  epsl = 0.000001

  do ieface=1,neface
     do iec=1,ne
        irng = rngface(ieface,iec)
        if(irng.ne.0) then
           do inface=1,nnface
              inl = mapping(ieface,inface,etype)
              do isd = 1,nsd
                 if(ids(isd,ien(inl,iec)) == 0)  ids(isd,ien(inl,iec)) = bcd(isd,irng)
              enddo
           enddo
        endif
     enddo
  enddo


  nq = 0
  do inn=1,nn
     do isd=1,nsd
        if(ids(isd,inn) == 0) then
           nq = nq + 1
           ids(isd,inn) = nq
        else
           ids(isd,inn) = 0
        endif
     enddo
  enddo


  return
end subroutine formidm


end module formm