

subroutine regen_points_element(xn,x_inter,x_center,hg,I_fluid_center,corr_Ip,&
				norm_e,nn_regen_proc,x_regen_proc,ie)

  use fluid_variables,only:nsd,nn,ne,nen
  use interface_variables
  use mpi_variables

  real(8) xn(nsd,4),x_inter(nsd,maxmatrix),x_center(nsd,nn_center),hg(ne)
  real(8) I_fluid_center(nn_center),corr_Ip(maxmatrix)
  real(8) norm_e(nsd)
  integer nn_regen_proc
  real(8) x_regen_proc(nsd,500)


  real(8) norm(nsd),tangx(nsd),tangy(nsd)

  integer i,j,k,icount,jcount,kcount,isd,ie
  real(8) temp,modmod

  real(8) xloc(nsd,4)

  integer, parameter :: nn_sub=3
  real(8) x_loc_can(nsd,nn_sub**2),sh(4),xlocan(nsd)
  real(8) x_tran(nsd),vec_temp(nsd)

  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) norm_p(nsd),curv_p,dcurv
  real(8) err_p,delta(nsd)
  integer nit
  real(8) distance_sign
  real(8) xlocan_ini(nsd)  

!goto 100
  x_loc_can(:,:)=0.0

  tangx(:)=xn(:,2)-xn(:,1)
  modmod=sqrt(tangx(1)**2+tangx(2)**2+tangx(3)**2)
  tangx(:)=tangx(:)/modmod

  tangy(:)=xn(:,3)-xn(:,2)
  modmod=sqrt(tangy(1)**2+tangy(2)**2+tangy(3)**2)
  tangy(:)=tangy(:)/modmod

  norm(1)=tangx(2)*tangy(3)-tangx(3)*tangy(2)
  norm(2)=tangx(3)*tangy(1)-tangx(1)*tangy(3)
  norm(3)=tangx(1)*tangy(2)-tangx(2)*tangy(1)

  modmod=sqrt(norm(1)**2+norm(2)**2+norm(3)**2)
  norm(:)=norm(:)/modmod

  temp=norm(1)*norm_e(1)+norm(2)*norm_e(2)+norm(3)*norm_e(3)
  if(temp.gt.0.0) distance_sign=1.0
  if(temp.le.0.0) distance_sign=-1.0

  tangy(1)=norm(2)*tangx(3)-norm(3)*tangx(2)
  tangy(2)=norm(3)*tangx(1)-norm(1)*tangx(3)
  tangy(3)=norm(1)*tangx(2)-norm(2)*tangx(1)

  modmod=sqrt(tangy(1)**2+tangy(2)**2+tangy(3)**2)
  tangy(:)=tangy(:)/modmod

  xloc(1,:)=tangx(1)*(xn(1,:)-xn(1,1))+tangx(2)*(xn(2,:)-xn(2,1))+tangx(3)*(xn(3,:)-xn(3,1))
  xloc(2,:)=tangy(1)*(xn(1,:)-xn(1,1))+tangy(2)*(xn(2,:)-xn(2,1))+tangy(3)*(xn(3,:)-xn(3,1))
  xloc(3,:)=norm(1) *(xn(1,:)-xn(1,1))+norm(2) *(xn(2,:)-xn(2,1))+norm(3) *(xn(3,:)-xn(3,1)) 

  do i=1,nn_sub
     do j=1,nn_sub
        x_loc_can(1,nn_sub*(i-1)+j)=2.0/nn_sub*i-1.0-1.0/nn_sub
        x_loc_can(2,nn_sub*(i-1)+j)=2.0/nn_sub*j-1.0-1.0/nn_sub
     end do
  end do
  do i=1,nn_sub**2
     sh(1)=0.25*(1-x_loc_can(1,i))*(1-x_loc_can(2,i))
     sh(2)=0.25*(1+x_loc_can(1,i))*(1-x_loc_can(2,i))
     sh(3)=0.25*(1+x_loc_can(1,i))*(1+x_loc_can(2,i))
     sh(4)=0.25*(1-x_loc_can(1,i))*(1+x_loc_can(2,i))

     xlocan(:)=0.0
     do j=1,4
	xlocan(1:2)=xlocan(1:2)+sh(j)*xloc(1:2,j)
     end do
     xlocan(3)=distance_sign*hg_sp/20.0  ! slightly lift the point

!     xlocan(:)=xlocan(:)+xn(:,1)
     x_tran(1)=tangx(1)*xlocan(1)+tangy(1)*xlocan(2)+norm(1)*xlocan(3)
     x_tran(2)=tangx(2)*xlocan(1)+tangy(2)*xlocan(2)+norm(2)*xlocan(3)
     x_tran(3)=tangx(3)*xlocan(1)+tangy(3)*xlocan(2)+norm(3)*xlocan(3)
     xlocan(:)=x_tran(:)+xn(:,1)
xlocan_ini(:)=xlocan(:)
!nn_regen_proc=nn_regen_proc+1
!x_regen_proc(:,nn_regen_proc)=xlocan(:)
!+++++++++++++start projection++++++++++++++!
     nit=1
     err_p=999.0
     delta(:)=0.0
     do while((nit.lt.5).and.(err_p.gt.1.0e-6))

              temp=sqrt(delta(1)**2+delta(2)**2+delta(3)**2)
              if(temp.gt.hg_sp) delta(:)=delta(:)/temp*hg_sp

	xlocan(1:nsd)=xlocan(1:nsd)+delta(1:nsd)
        call get_indicator_derivative_3D_1st(xlocan,x_inter,x_center,hg,&
                    I_fluid_center,corr_Ip,II,dI,ddI,norm_p,curv_p)
!write(*,'(6f6.2)')norm_e(:),norm_p(:)

        vec_temp(1)=norm_e(2)*norm_p(3)-norm_e(3)*norm_p(2)
        vec_temp(2)=norm_e(3)*norm_p(1)-norm_e(1)*norm_p(3)
        vec_temp(3)=norm_e(1)*norm_p(2)-norm_e(2)*norm_p(1)

        norm_p(1)=vec_temp(2)*norm_e(3)-vec_temp(3)*norm_e(2)
        norm_p(2)=vec_temp(3)*norm_e(1)-vec_temp(1)*norm_e(3)
        norm_p(3)=vec_temp(1)*norm_e(2)-vec_temp(2)*norm_e(1)

        modmod=sqrt(norm_p(1)**2+norm_p(2)**2+norm_p(3)**2)
        norm_p(:)=norm_p(:)/modmod
        temp=dI(1)*norm_p(1)+dI(2)*norm_p(2)+dI(3)*norm_p(3)
if(abs(temp).lt.1.0e-6)write(*,*)'something wrong with points regen in ele'
        delta(:)=(0.5-II)*norm_p(:)/temp
        nit=nit+1
        err_p=abs(II-0.5)
     end do
!+++++++++++++end of projection+++++++++++++++!
     if(err_p.lt.1.0e-6) then
                call get_curv_num_3D(xlocan,x_inter,x_center,hg,I_fluid_center,&
                        corr_Ip,dcurv,curv_p,norm_p)
      if(dcurv.lt.20.0) then
       nn_regen_proc=nn_regen_proc+1
       x_regen_proc(:,nn_regen_proc)=xlocan(:)

      end if
     end if
  end do !! end do of loop over nn_sub**2
100 continue
end subroutine regen_points_element
