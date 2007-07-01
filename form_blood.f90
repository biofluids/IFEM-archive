!this module is used to apply boundary conditions
module form
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formd(ds,rngface,ien,nodes_BC_fluid)
  use global_constants
  use run_variables, only: tt
  use fluid_variables, only: nn,ne,ndf,nsd,nen,neface,nrng,nnface,mapping,bc,bv,etype,ic,static,udf,vdf,wdf, maxconn
  use solid_variables, only: nn_solid

  implicit none

  integer  rngface(neface,ne),ien(nen,ne)
  real(8) :: ds(ndf,nn) !meshvel(nsd,nn)
  integer :: idf, inl, iec, irng, ieface, inface, inn
  real(8) :: eps1,eps2 !,tt_ramp
  real(8) :: hs(nrng+1,nn), h(nrng+1,nn)
  integer :: T0, k
  integer :: nodes_BC_fluid(1:nn_solid,1:maxconn)
  real(8) :: A10, A11, A12, A13, A14, A20, A21, A22, A23, A24, A30, A31, A32, A33, A34
  real(8) :: A40, A41, A42, A43, A44, Area

  eps1 = -1000000.0 
  eps2 = -10000.0 

  h(:,:) = 0.0d0

  ds(:,:) = eps1
  


  do ieface=1,neface
     do inface=1,nnface
        inl = mapping(ieface,inface,etype)
        do iec=1,ne
           irng = rngface(ieface,iec)
           if (irng.ne.0) h(irng,ien(inl,iec)) = h(irng,ien(inl,iec)) + 1.0
        enddo
     enddo
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mickael 02/28/2005
! BC on influence nodes
 ! do inn=1,nn_solid
!	do k=1,maxconn
!	  if (nodes_BC_fluid(inn,k)/=0) then
!	  	irng=nrng+1
!	    h(irng,nodes_BC_fluid(inn,k)) = h(irng,nodes_BC_fluid(inn,k)) + 1.0
!	  endif
!	enddo
 ! enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  hs = h
!!!!!!!!!!!!!
!Yili Gu


!bv(udf, 1) = 0.5*sin(pi*tt)
!bv(udf, 5) = 0.5*sin(pi*tt)

!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!
! Blood flow
Area=19.634954d-2 !( pi*r^2, r=2.5 mm)
        
  ! Parameters

  ! 0<t<0.2
  A10=4d-4
  A11=-42.5d-4
  A12=-1241.66667d-4
  A13=37000d-4
  A14=-143333.3333d-4

  ! 0.2<t<0.45
  A20=708.55357d-4
  A21=-9049.52381d-4
  A22=42905.83335d-4
  A23=-88833.33337d-4
  A24=67666.66669d-4

  ! 0.45<0.65
  A30=-507.74964d-4
  A31=2416.33065d-4
  A32=-3004.99258d-4
  A33=-533.34237d-4
  A34=2000.00411d-4

  ! 0.65<0.9
  A40=-5.4204d-4
  A41=13.14305d-4
  A42=5.13981d-4
  A43=3.85125d-4
  A44=-14.29531d-4


  ! 11 CARDIAC CYCLES
  if (tt<=0.9) then
        T0=0
  elseif ((tt>0.9).and.(tt<=1.8)) then
        T0=1
  elseif ((tt>1.8).and.(tt<=2.7)) then
        T0=2
  elseif ((tt>2.7).and.(tt<=3.6)) then
        T0=3
  elseif ((tt>3.6).and.(tt<=4.5)) then
        T0=4
  elseif ((tt>4.5).and.(tt<=5.4)) then
        T0=5
  elseif ((tt>5.4).and.(tt<=6.3)) then
        T0=6
  elseif ((tt>6.3).and.(tt<=7.2)) then
        T0=7
  elseif ((tt>7.2).and.(tt<=8.1)) then
        T0=8
  elseif ((tt>8.1).and.(tt<=9.0)) then
        T0=9
  elseif ((tt>9.0).and.(tt<=9.9)) then
        T0=10
  elseif ((tt>9.9).and.(tt<=10.8)) then
        T0=11
  elseif ((tt>10.8).and.(tt<=11.7)) then
        T0=12
  endif

  if ((tt >= (T0*0.9+0.00)).and.(tt <= (T0*0.9+0.2))) then
                bv(udf,1)=(A10+A11*(tt-T0*0.9)+A12*((tt-T0*0.9)**2)+ &
                A13*((tt-T0*0.9)**3)+A14*((tt-T0*0.9)**4))/Area
  elseif ((tt > (T0*0.9+0.2)).and.(tt <= (T0*0.9+0.45))) then
                bv(udf,1)=(A20+A21*(tt-T0*0.9)+A22*((tt-T0*0.9)**2)+ &
                A23*((tt-T0*0.9)**3)+A24*((tt-T0*0.9)**4))/Area
  elseif ((tt > (T0*0.9+0.45)).and.(tt <= (T0*0.9+0.65))) then
                bv(udf,1)=(A30+A31*(tt-T0*0.9)+A32*((tt-T0*0.9)**2)+ &
               A33*((tt-T0*0.9)**3)+A34*((tt-T0*0.9)**4))/Area
  elseif ((tt > (T0*0.9+0.65)).and.(tt <= (T0*0.9+0.9))) then
                bv(udf,1)=(A40+A41*(tt-T0*0.9)+A42*((tt-T0*0.9)**2)+ &
                A43*((tt-T0*0.9)**3)+A44*((tt-T0*0.9)**4))/Area
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
 
do irng=1,nrng
     do inn=1,nn
        do idf=1,ndf
           if (hs(irng,inn).gt.1.0e-8) then
              if (bc(idf,irng) .gt. 0) then
                ds(idf,inn) = bv(idf,irng)
			  endif
           endif
        enddo
     enddo
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mickael 02/28/2005
! BC on influence nodes
 ! do inn=1,nn_solid
!	do k=1,maxconn
!	  if (nodes_BC_fluid(inn,k)/=0) then
!	  	irng=nrng+1
!		do idf=1,ndf
!			if (hs(irng,nodes_BC_fluid(inn,k)).gt.1.0e-8) then
 !             IF (BC(IDF,IRng) .gt. 0) then
  !              ds(idf,nodes_BC_fluid(inn,k)) = bv(idf,irng)
!			  endif
 !           endif
!		enddo
!	  endif
!	enddo
 ! enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  do inn=1,nn
     do idf=1,ndf
        if(ds(idf,inn).lt.eps2) then
			ds(idf,inn) = ic(idf)
		endif
     enddo
  enddo




  if(static) ds(:,:)=0.0

  return
end subroutine formd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formid(ids, rngface, ien, nodes_BC_fluid)
  use fluid_variables
  use solid_variables, only: nn_solid

  implicit none

  integer :: ids(ndf,nn),rngface(neface,ne),ien(nen,ne)
  integer :: idf, inl, iec, irng, ieface, inface, inn,k
  integer :: nodes_BC_fluid(1:nn_solid,1:maxconn)

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
                 if(d(idf,ien(inl,iec)).lt.epsr) then
					d(idf,ien(inl,iec)) = bc(idf,irng)+epsl
					
				  endif
              enddo
           endif
        enddo
     enddo
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mickael 02/28/2005
! BC on influence nodes
!  do inn=1,nn_solid
!	do k=1,maxconn
!	  if (nodes_BC_fluid(inn,k)/=0) then
!	  	irng=nrng+1
!		do idf=1,ndf
!			if(d(idf,nodes_BC_fluid(inn,k)).lt.epsr) then
!				d(idf,nodes_BC_fluid(inn,k)) = bc(idf,irng)+epsl
 !   		endif
!		enddo
!	  endif
!	enddo
 ! enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ds = d
  ids = ds



  if(static) ids(:,:) = 1

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


end module form
