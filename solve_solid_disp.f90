subroutine solve_solid_disp(x,x_curr,kid,ien,node_sbc,xpre1,solid_prevel,solid_preacc,ien_sbc,&
			solid_stress,solid_bcvel,mtype)
! subroutine to solve the solid displacement at every time step
! x  -- initial solid configuration
! x_curr --- current solid configuration
! disp --- solid displacement 
! kid --- 1st type boundary condition nodes id
! ien --- solid connective matrix
! xpre1 --- last 1 time step solid position
! xpre2 --- last 2 time step solid position
! solid_acc --- solid acceleration at n+1
! solid_preacc --- solid acceleration at n
! solid_vel --- solid velocity at n+1
! solid_prevel --- solid velocity at n
! solid_stress --- solid traction B.C.
! solid_bcvel --- solid velocity B.C.

use mpi_variables, only: myid
use solid_variables
use run_variables, only: dt
include 'mpif.h'
!implicit none

real(8) x(nsd_solid,nn_solid)
real(8) x_curr(nsd_solid,nn_solid)
real(8) disp(nsd_solid,nn_solid)
integer kid(nsd_solid,nn_solid)
integer ien(ne_solid,nen_solid)
integer node_sbc(nn_sbc)
real(8) xpre1(nsd_solid,nn_solid)
real(8) solid_prevel(nsd_solid,nn_solid)
real(8) solid_preacc(nsd_solid,nn_solid)
real(8) solid_acc(nsd_solid,nn_solid)
real(8) solid_vel(nsd_solid,nn_solid)
real(8) solid_bcvel(nsd_solid,nn_solid)
integer mtype(ne_solid)
!real(8) xpre2(nsd_solid,nn_solid)
!real(8) disp1(nsd_solid,nn_solid)
!real(8) disp2(nsd_solid,nn_solid)
integer ien_sbc(ne_sbc,nen_solid+2)
real(8) solid_stress(nsd_solid*2,nn_solid)
!-------------------------------------------------
real(8) sq_solid(0:3,8,8)
integer iq
!real(8) :: xq_solid(nsdpad_solid,nquadpad_solid),wq(nquadpad_solid)
!-------------------------------------------------
integer inner
integer outer
parameter (inner = 100) ! solid equation inner 50 should be sufficient
parameter (outer = 10)  ! solid equation outer 5 should be sufficient
!-------------------------------------------------
real(8) dg(nsd_solid,nn_solid) ! disp correction
real(8) w(nsd_solid,nn_solid)  ! pre-conditioner
real(8) p(nsd_solid,nn_solid)  ! residual vector
real(8) res
real(8) del
integer i
real(8) time
!----------------------------
        real(8) alpha
        real(8) beta
        real(8) gama
!----------------------------
! define the numerical parameters
alpha = -0.05 ! -1/3 < alpha < 0 and alpha == 0 is Newmark method
gama = (1.0 - 2.0 * alpha) * 0.5
beta = ( (1.0 - alpha)**2 ) * 0.25
!beta = (1.0 - alpha**2 ) * 0.25
!----------------------------------

! Get sq for solid mesh to calculate sh in block_solid.f90
if (nsd_solid == 3) then
     do iq=1,nquad_solid
                if(nen_solid.eq.4) then
                  sq_solid(0,1,iq) = xq_solid(1,iq)
                  sq_solid(0,2,iq) = xq_solid(2,iq)
                  sq_solid(0,3,iq) = xq_solid(3,iq)
                  sq_solid(0,4,iq) = 1 - xq_solid(1,iq) - xq_solid(2,iq) - xq_solid(3,iq)
        else
                  sq_solid(0,1,iq) = (1 - xq_solid(1,iq))   &
                           * (1 - xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(0,2,iq) = (1 + xq_solid(1,iq))   &
                           * (1 - xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(0,3,iq) = (1 + xq_solid(1,iq))   &
                           * (1 + xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(0,4,iq) = (1 - xq_solid(1,iq))   &
                           * (1 + xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(0,5,iq) = (1 - xq_solid(1,iq))   &
                           * (1 - xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(0,6,iq) = (1 + xq_solid(1,iq))   &
                           * (1 - xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(0,7,iq) = (1 + xq_solid(1,iq))   &
                           * (1 + xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(0,8,iq) = (1 - xq_solid(1,iq))   &
                           * (1 + xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(1,1,iq) = - (1 - xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(1,2,iq) = + (1 - xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(1,3,iq) = + (1 + xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(1,4,iq) = - (1 + xq_solid(2,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(1,5,iq) = - (1 - xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(1,6,iq) = + (1 - xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(1,7,iq) = + (1 + xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(1,8,iq) = - (1 + xq_solid(2,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(2,1,iq) = - (1 - xq_solid(1,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(2,2,iq) = - (1 + xq_solid(1,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(2,3,iq) = + (1 + xq_solid(1,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(2,4,iq) = + (1 - xq_solid(1,iq)) * (1 - xq_solid(3,iq)) / 8
                  sq_solid(2,5,iq) = - (1 - xq_solid(1,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(2,6,iq) = - (1 + xq_solid(1,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(2,7,iq) = + (1 + xq_solid(1,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(2,8,iq) = + (1 - xq_solid(1,iq)) * (1 + xq_solid(3,iq)) / 8
                  sq_solid(3,1,iq) = - (1 - xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 8
                  sq_solid(3,2,iq) = - (1 + xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 8
                  sq_solid(3,3,iq) = - (1 + xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 8
                  sq_solid(3,4,iq) = - (1 - xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 8
                  sq_solid(3,5,iq) = + (1 - xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 8
                  sq_solid(3,6,iq) = + (1 + xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 8
                  sq_solid(3,7,iq) = + (1 + xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 8
                  sq_solid(3,8,iq) = + (1 - xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 8
        endif
          enddo
elseif (nsd_solid == 2) then
                  
      do iq=1,nquad_solid
                if(nen_solid==3) then
                  sq_solid(0,1,iq) = xq_solid(1,iq)
                  sq_solid(0,2,iq) = xq_solid(2,iq)
                  sq_solid(0,3,iq) = 1 - xq_solid(1,iq) - xq_solid(2,iq) 
        elseif (nen_solid==4) then 
                  sq_solid(0,1,iq) = (1 - xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 4
                  sq_solid(0,2,iq) = (1 + xq_solid(1,iq)) * (1 - xq_solid(2,iq)) / 4
                  sq_solid(0,3,iq) = (1 + xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 4
                  sq_solid(0,4,iq) = (1 - xq_solid(1,iq)) * (1 + xq_solid(2,iq)) / 4
                           
                  sq_solid(1,1,iq) = - (1 - xq_solid(2,iq)) / 4
                  sq_solid(1,2,iq) = + (1 - xq_solid(2,iq)) / 4
                  sq_solid(1,3,iq) = + (1 + xq_solid(2,iq)) / 4
                  sq_solid(1,4,iq) = - (1 + xq_solid(2,iq)) / 4
                  
                  sq_solid(2,1,iq) = - (1 - xq_solid(1,iq)) / 4
                  sq_solid(2,2,iq) = - (1 + xq_solid(1,iq)) / 4
                  sq_solid(2,3,iq) = + (1 + xq_solid(1,iq)) / 4
                  sq_solid(2,4,iq) = + (1 - xq_solid(1,iq)) / 4
                  
        endif     
          enddo  
end if
!------------------------------------------------------------------
solid_acc(:,:) = 0.0d0
solid_vel(:,:) = 0.0d0
solid_acc(:,:) = solid_preacc(:,:)

!do i=1,nn_sbc     
!        solid_acc(:,node_sbc(i))=(solid_bcvel(:,node_sbc(i)) - solid_prevel(:,node_sbc(i))) / dt
!end do

w(:,:)=0.0d0
p(:,:)=0.0d0
dg(:,:)=0.0d0 

call block_solid(x,solid_acc,w,p,ien,nsd_solid,nen_solid,ne_solid,&
		nn_solid,nquad_solid,wq_solid,sq_solid,xpre1,&
		solid_prevel,solid_preacc,ien_sbc,ne_sbc,solid_stress,mtype)
!do i=1,nn_sbc     
!        p(1,node_sbc(i))=p(1,node_sbc(i)) + x(1,node_sbc(i))*0.01
!end do
call setsolid_id(p,kid,nsd_solid)
call getnorm(p,p,nsd_solid*nn_solid,res)
res=sqrt(res)
if (myid == 0) write(*,*) '===Initial error for solid displacement===', res

!if ( res .gt. 1e-6) then  ! solid disp need to be solved

!-----------------------
! Take the invese of w as the preconditioner
! Help a lot!!!! However, have not figured out why ...
!w(:,:)=1.0d0/w(:,:)
!----------------------
!time=mpi_wtime()
call gmres_solid(x,w,p,dg,ien,kid,nsd_solid,nn_solid,ne_solid,nen_solid,inner,outer,&
		nquad_solid,wq_solid,sq_solid,xpre1,&
		solid_prevel,solid_preacc,solid_stress,ne_sbc,ien_sbc,mtype)
!time=mpi_wtime()-time
!if (myid==0) write(*,*) '*****COSTING TIME********', time

call getnorm(dg,dg,nsd_solid*nn_solid,del)
 del = sqrt(del)
if (myid == 0) write(*,*) '===solid displacement correction norm===', del

!do i=1,nn_sbc
!	write(*,*) 'dg', dg(:,node_sbc(i)), ' at node',node_sbc(i), disp(:,node_sbc(i))
!        dg(:,node_sbc(i))=disp(:,node_sbc(i))
!end do
solid_acc(:,:) = solid_acc(:,:) + dg(:,:)
x_curr(:,:) = xpre1(:,:) + dt*solid_prevel(:,:) +&
		 (dt**2)*0.5*( (1.0-2.0*beta)*solid_preacc(:,:) + 2.0*beta*solid_acc(:,:) )
solid_vel(:,:) = solid_prevel(:,:) + dt*( (1-gama)*solid_preacc(:,:) + gama*solid_acc(:,:) )

!do i=1,nn_sbc
!        solid_vel(:,node_sbc(i))=solid_bcvel(:,node_sbc(i))
!end do


solid_prevel(:,:) = solid_vel(:,:)
solid_preacc(:,:) = solid_acc(:,:)


return
end 
