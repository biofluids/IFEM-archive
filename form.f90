module form
  implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formd(ds,meshvel,rngface,ien)
  use global_constants
  use run_variables, only: tt
  use fluid_variables, only: nn,ne,ndf,nsd,nen,neface,nrng,nnface,mapping,bc,bv,etype,ic,static,udf,vdf,wdf

  implicit none

  integer  rngface(neface,ne),ien(nen,ne)
  !real(8) :: xn(nsd,nn),
  real(8) :: ds(ndf,nn),meshvel(nsd,nn)
  integer :: idf, inl, iec, irng, ieface, inface, inn

  real(8) :: eps1,eps2,tt_ramp

  real(8) :: hs(nrng,nn), h(nrng,nn)


  eps1 = -1000000.0 
  eps2 = -10000.0 

  h(:,:) = 0.0d0

  !ds(1:ndf,1:nn) = eps1
  

  do ieface=1,neface
     do inface=1,nnface
        inl = mapping(ieface,inface,etype)
        do iec=1,ne
           irng = rngface(ieface,iec)
           if (irng.ne.0) then
              h(irng,ien(inl,iec)) = h(irng,ien(inl,iec)) + 1.0
           endif
        enddo
     enddo
  enddo

  !call equal(h,hs,nrng*nn)
  hs(:,:) = h(:,:)

  tt_ramp = 0.4
  if ((tt >= 0).and.(tt <= tt_ramp)) then
       bv(udf,2) = 10.0d0  * tt/tt_ramp
       bv(udf,4) = 10.0d0  * tt/tt_ramp
       bv(udf,5) = 10.0d0  * tt/tt_ramp
       bv(udf,6) = 10.0d0  * tt/tt_ramp  !...inlet (6) velocity in direction #degree of freedom#
       bv(wdf,2) = 2.50d0  * tt/tt_ramp
       bv(wdf,4) = 2.50d0  * tt/tt_ramp
       !bv(wdf,5) = 1.0d0  * tt/tt_ramp
       !bv(wdf,6) = 1.0d0  * tt/tt_ramp
  elseif (tt > tt_ramp) then
       bv(udf,2) = 10.0d0
       bv(udf,4) = 10.0d0
       bv(udf,5) = 10.0d0
       bv(udf,6) = 10.0d0
       bv(wdf,2) = 2.50d0
       bv(wdf,4) = 2.50d0
       !bv(wdf,5) = 1.0d0
       !bv(wdf,6) = 1.0d0
  endif



  write(*,'("boundary 6 (inflow) x-velocity: ",f12.7," cm/s")')bv(udf,6)
  write(*,'("boundary 6 (inflow) z-velocity: ",f12.7," cm/s")')bv(wdf,6)


  do irng=1,nrng
     do inn=1,nn
        do idf=1,ndf
           if (hs(irng,inn).gt.1.0e-8) then
              if (bc(idf,irng) .gt. 0) then
                 if (idf <= 3) then
                    ds(idf,inn) = bv(idf,irng) !+ meshvel(idf,inn)
                 endif
              endif
           endif
        enddo
     enddo
  enddo

  do inn=1,nn
     do idf=1,ndf
        if(ds(idf,inn).lt.eps2) ds(idf,inn) = ic(idf)
     enddo
  enddo

  if(static) then
     do inn = 1,nn
        do idf = 1,nsd 
           ds(idf,inn) = 0.0
        enddo 
     enddo
  endif


  return
end subroutine formd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formid(ids, rngface, ien)
  use fluid_variables
  implicit none

  integer :: ids(ndf,nn),rngface(neface,ne),ien(nen,ne)
  integer :: idf, inl, iec, irng, ieface, inface, inn

  real(8) :: epsr,epsl

  real(8) :: ds(ndf,nn),d(ndf,nn)

  d(1:ndf,1:nn) = 0.0d0

  epsr = 0.0001
  epsl = 0.000001

  do ieface=1,neface
     do inface=1,nnface
        inl = mapping(ieface,inface,etype)
        do iec=1,ne
           irng = rngface(ieface,iec)
           if(irng.ne.0) then
              do idf = 1,ndf
                 if(d(idf,ien(inl,iec)).lt.epsr) d(idf,ien(inl,iec)) = bc(idf,irng)+epsl
              enddo
           endif
        enddo
     enddo
  enddo

  ds(:,:) = d(:,:)

  ids(1:ndf,1:nn) = ds(1:ndf,1:nn)
  

  if(static) then
     do inn = 1,nn
        do idf = 1,nsd 
           ids(idf,inn) = 1
        enddo 
     enddo
  endif

  nq = 0

  do inn=1,nn
     do idf=1,ndf
        if(ids(idf,inn).eq.0) then
           nq = nq + 1
           ids(idf,inn) = nq
        else
           ids(idf,inn) = 0
        endif
     enddo
  enddo

  return
end subroutine formid


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
           if(irng.ne.0) then
              h(irng,ien(inl,iec)) = 1  
           endif
        enddo
     enddo
  enddo

  do irng=1,nrng
     do inn=1,nn
        do isd=1,nsd
           if (h(irng,inn) == 1) then
              if (bcd(isd,irng) > 0) then
                 ds(isd,inn) = bvd(isd,irng)
              endif
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




end module form


