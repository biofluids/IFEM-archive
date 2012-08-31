!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!calculate total mass!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  cal_mass(xloc,I_fluid,ien,nn,ne,nen,mass)
  use global_constants
  use run_variables
  use fluid_variables,only:nsd,iquad
  use mpi_variables
  include 'mpif.h'
!  use centermesh_variables
!  use denmesh_variables, only:nn_den,ne_den,nen_den
!  implicit none

  integer nn,ne,nen
  integer ien(nen,ne)
  real(8) xloc(nsd,nn)
  real(8) I_fluid(nn),mass

  integer ne_loc,base,top,loc_index

  real(8) x(nsd,nen)
  real(8) eft0,det,effd,effm,effc
  real(8) sh(0:nsd,nen),ph(0:nsd,nen)
  real(8) xr(nsd,nsd),cf(nsd,nsd),sx(nsd,nsd)

  real(8) hg,temp
  real(8) M1,M1_temp
  
  integer inl,ie,isd,node,iq,ine

  integer,parameter :: ndfpad=5,nsdpad=3,nenpad=8,nquadpad=8
  integer nquad
  real(8) sq(0:nsdpad,nenpad,nquadpad)
  real(8) xq(nsdpad,nquadpad),wq(nquadpad)

  if(nsd==3) then
!     write(*,*)'no 3d in block laplace'
!     stop
          if (nen.eq.4) then
                call quad3d4n(iquad, nquad, xq, wq, nsdpad, nquadpad)
          else if (nen.eq.8) then
                call quad3d8n(iquad, nquad, xq, wq, nsdpad, nquadpad)
          end if
      do iq=1,nquad
                if(nen.eq.4) then
                  sq(0,1,iq) = xq(1,iq)
                  sq(0,2,iq) = xq(2,iq)
                  sq(0,3,iq) = xq(3,iq)
                  sq(0,4,iq) = 1 - xq(1,iq) - xq(2,iq) - xq(3,iq)
               else
                  sq(0,1,iq) = (1 - xq(1,iq))   &
                           * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
                  sq(0,2,iq) = (1 + xq(1,iq))   &
                           * (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
                  sq(0,3,iq) = (1 + xq(1,iq))   &
                           * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
                  sq(0,4,iq) = (1 - xq(1,iq))   &
                           * (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
                  sq(0,5,iq) = (1 - xq(1,iq))   & 
                           * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8 
                  sq(0,6,iq) = (1 + xq(1,iq))   & 
                           * (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8 
                  sq(0,7,iq) = (1 + xq(1,iq))   & 
                           * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8 
                  sq(0,8,iq) = (1 - xq(1,iq))   & 
                           * (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8 
                  sq(1,1,iq) = - (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
                  sq(1,2,iq) = + (1 - xq(2,iq)) * (1 - xq(3,iq)) / 8
                  sq(1,3,iq) = + (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
                  sq(1,4,iq) = - (1 + xq(2,iq)) * (1 - xq(3,iq)) / 8
                  sq(1,5,iq) = - (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
                  sq(1,6,iq) = + (1 - xq(2,iq)) * (1 + xq(3,iq)) / 8
                  sq(1,7,iq) = + (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
                  sq(1,8,iq) = - (1 + xq(2,iq)) * (1 + xq(3,iq)) / 8
                  sq(2,1,iq) = - (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
                  sq(2,2,iq) = - (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
                  sq(2,3,iq) = + (1 + xq(1,iq)) * (1 - xq(3,iq)) / 8
                  sq(2,4,iq) = + (1 - xq(1,iq)) * (1 - xq(3,iq)) / 8
                  sq(2,5,iq) = - (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
                  sq(2,6,iq) = - (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
                  sq(2,7,iq) = + (1 + xq(1,iq)) * (1 + xq(3,iq)) / 8
                  sq(2,8,iq) = + (1 - xq(1,iq)) * (1 + xq(3,iq)) / 8
                  sq(3,1,iq) = - (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
                  sq(3,2,iq) = - (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
                  sq(3,3,iq) = - (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
                  sq(3,4,iq) = - (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
                  sq(3,5,iq) = + (1 - xq(1,iq)) * (1 - xq(2,iq)) / 8
                  sq(3,6,iq) = + (1 + xq(1,iq)) * (1 - xq(2,iq)) / 8
                  sq(3,7,iq) = + (1 + xq(1,iq)) * (1 + xq(2,iq)) / 8
                  sq(3,8,iq) = + (1 - xq(1,iq)) * (1 + xq(2,iq)) / 8
                 endif
          enddo



  else if(nsd==2) then
!=====================================================================
  if (nen==3) then
       call quad2d3n(iquad, nquad, xq, wq, nsdpad, nquadpad)
  else if (nen==4) then
       call quad2d4n(iquad, nquad, xq, wq, nsdpad, nquadpad)
  end if
  do iq=1,nquad
       if(nen==3) then
                  sq(0,1,iq) = xq(1,iq)
                  sq(0,2,iq) = xq(2,iq)
                  sq(0,3,iq) = 1 - xq(1,iq) - xq(2,iq)
        elseif (nen==4) then
                  sq(0,1,iq) = (1 - xq(1,iq)) * (1 - xq(2,iq)) / 4
                  sq(0,2,iq) = (1 + xq(1,iq)) * (1 - xq(2,iq)) / 4
                  sq(0,3,iq) = (1 + xq(1,iq)) * (1 + xq(2,iq)) / 4
                  sq(0,4,iq) = (1 - xq(1,iq)) * (1 + xq(2,iq)) / 4
                
                  sq(1,1,iq) = - (1 - xq(2,iq)) / 4
                  sq(1,2,iq) = + (1 - xq(2,iq)) / 4
                  sq(1,3,iq) = + (1 + xq(2,iq)) / 4
                  sq(1,4,iq) = - (1 + xq(2,iq)) / 4
        
                  sq(2,1,iq) = - (1 - xq(1,iq)) / 4
                  sq(2,2,iq) = - (1 + xq(1,iq)) / 4
                  sq(2,3,iq) = + (1 + xq(1,iq)) / 4
                  sq(2,4,iq) = + (1 - xq(1,iq)) / 4

        endif
  enddo
!==========================================================
  end if

  if(ne.le.ncpus) then
    if(myid+1.le.ne) then
      ne_loc=1
    else
      ne_loc=0
    end if
  else
    base=floor(real(ne)/real(ncpus))
    top=ne-base*ncpus
    if(myid+1.lt.top) then
      ne_loc=base+1
    else
      ne_loc=base
    end if
  end if

  M1_temp=0.0
  M1=0.0
  do loc_index=1,ne_loc
     ie=myid+1+(loc_index-1)*ncpus

     do inl=1,nen
        x(1:nsd,inl)=xloc(1:nsd,ien(inl,ie))
     end do



     do iq=1,nquad  ! loop over the quadrature points in each element
!!!!!!!!!!!!!!!!!!!!calculate shape function!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (nsd==2) then
	   if(nen.eq.3) then
		include 'sh2d3n.h'
	   elseif(nen.eq.4) then
		include 'sh2d4n.h'
	   end if
	elseif(nsd==3) then
	   if(nen.eq.4) then
		include 'sh3d4n.h'
	   elseif(nen.eq.8) then
		include 'sh3d8n.h'
	   end if
	end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     eft0 = abs(det)*wq(iq)   ! calculate weight at each quad point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ph(0,1:nen)=sh(0,1:nen)*eft0
     temp=0
     do inl=1,nen
        node=ien(inl,ie)
        temp=temp+sh(0,inl)*I_fluid(node)
     end do
     do inl=1,nen
       node=ien(inl,ie)
       M1_temp=M1_temp+ph(0,inl)*temp
     end do


     end do ! end of quad loop

  end do ! end of ele loop
  call mpi_barrier(mpi_comm_world,ierror)
  call mpi_allreduce(M1_temp,M1,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
  call mpi_barrier(mpi_comm_world,ierror)
mass=M1
if(myid==0)write(*,*)'mass=',mass
!write(*,*)'w_inter=',w_inter(:)
end subroutine cal_mass










































