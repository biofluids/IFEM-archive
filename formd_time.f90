!=====================================================================
! Latest update: J. Yang, 05/20/2014
! This subroutine defines Dirichlet BC with time-varying velocity
!=====================================================================
subroutine formd_time(ds,rngface,ien)
  use global_constants
  use run_variables, only: tt ! physical time
  use fluid_variables, only: nn,ne,ndf,nsd,nen,neface,nrng,nnface,mapping,bc,bv,etype,ic,static,udf,vdf,wdf, maxconn
  implicit none
  integer  rngface(neface,ne),ien(nen,ne)
  real(8) :: ds(ndf,nn)
  integer :: idf, inl, iec, irng, ieface, inface, inn
!  real(8) :: eps1,eps2                    ! OBSOLETE
  real(8) :: hs(nrng+1,nn), h(nrng+1,nn)
!  eps1 = -1000000.0                       ! OBSOLETE
!  eps2 = -10000.0                         ! OBSOLETE
  h(:,:) = 0.0d0

    do ieface=1,neface
        do inface=1,nnface
            inl = mapping(ieface,inface,etype)
            do iec=1,ne
                irng = rngface(ieface,iec)
                if (irng.ne.0) h(irng,ien(inl,iec)) = h(irng,ien(inl,iec)) + 1.0
            enddo
        enddo
    enddo
    hs(:,:) = h(:,:)

!==========================================
! Jack 08/05/2013
! Velocity boundary condition with as a function of time
! bv(xdf,Nbc) is the velocity in xdf-direction on Boundary #Nbc
if (tt .lt. 1.0e-4) then
    bv(udf,4) = 6.0*exp(-0.5*((tt-0.5e-4)/0.15e-4)**2)
    !bv(udf,4) = 0.25*cos(3.1415927*tt/1.0e-4)
else
    bv(udf,4) = 0.0
endif
!==========================================

    do irng=1,nrng
        do inn=1,nn
            do idf=1,ndf
                if (hs(irng,inn) .gt. 1.0e-8) then
                    if (bc(idf,irng) .gt. 0) then
                        ds(idf,inn) = bv(idf,irng)
                    endif
                endif
            enddo
        enddo
    enddo

return
end subroutine formd_time

