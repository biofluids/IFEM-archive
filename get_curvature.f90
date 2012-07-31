subroutine get_curvature(x,x_inter,x_center,I_fluid,corr_Ip,I_fluid_center,sur_fluid,hg,curv_nn)

  use fluid_variables, only:nsd,ne,nn,nen,den_liq
  use interface_variables
  use mpi_variables
  include 'mpif.h'

  real(8) x(nsd,nn),x_inter(nsd,maxmatrix),x_center(nsd,ne)
  real(8) I_fluid(nn),corr_Ip(maxmatrix),I_fluid_center(ne)
  real(8) sur_fluid(nsd,nn),sur_fluid_temp(nsd,nn)

  real(8) den_p,den_f
  real(8) eps,hg(ne)

  real(8) II,dI(nsd),ddI(3*(nsd-1))
  real(8) curv_a,norm_a(nsd)
  real(8) curv_nn(nn),curv_nn_temp(nn)
  integer flag
  integer i,j,nit
  real(8) delta(nsd),x_ini(nsd),dis,temp,normx(1:nsd),xlocan(nsd)
  real(8) coef
  integer nn_loc,base,top,loc_index

  if(nn.le.ncpus) then
    if(myid+1.le.nn) then
      nn_loc=1
    else
      nn_loc=0
    end if
  else
    base=floor(real(nn)/real(ncpus))
    top=nn-base*ncpus
    if(myid+1.le.top) then
      nn_loc=base+1
    else
      nn_loc=base
    end if
  end if

  sur_fluid_temp(:,:)=0.0
  den_p=0.5*(den_inter+den_liq)
  eps=0.005
  curv_nn_temp(:)=0.0

  do loc_index=1,nn_loc
     i=myid+1+(loc_index-1)*ncpus
     if( (I_fluid(i).lt.(1.0-eps)) .and. (I_fluid(i).gt.eps)) then
      if(nsd==2) then
       call get_indicator_derivative_2D_1st(x(:,i),x_inter,x_center,hg,I_fluid_center,corr_Ip, &
                                        II,dI,ddI,norm_a,curv_a)
       normx(1:nsd)=dI(1:nsd)
!       temp=dI(1)**2+dI(2)**2
       nit=1
       delta(:)=0.0
       x_ini(:)=x(:,i)
       xlocan(:)=x(:,i)
       err_p=999.0
!if(i==5623)write(*,*)'x5355=',x(:,i)
       do while((nit.le.15).and.(err_p.gt.1.0e-7))
temp=sqrt(delta(1)**2+delta(2)**2)
coef=1.0
if(temp.gt.max_hg) coef=max_hg/temp

          xlocan(1:nsd)=xlocan(1:nsd)+delta(1:nsd)*coef
          call get_indicator_derivative_2D_1st(xlocan,x_inter,x_center,hg,&
                           I_fluid_center,corr_Ip,II,dI,ddI,norm_a,curv_a)
!          temp=normx(1)*dI(1)+normx(2)*dI(2)
!          delta(1)=(0.5-II)*normx(1)/temp
!          delta(2)=(0.5-II)*normx(2)/temp
           temp=dI(1)**2+dI(2)**2
           delta(1)=(0.5-II)*dI(1)/temp
           delta(2)=(0.5-II)*dI(2)/temp
           
          err_p=abs(II-0.5)
!if(i==5623)write(*,*)'mmmmmmmmmmmmmmmmmmmmm'
!if(i==5623)write(*,*)'nit=',nit
!if(i==5623)write(*,*)'dI=',dI(1:nsd)
!if(i==5623)write(*,*)'temp=',temp
!if(i==5623)write(*,*)'err_p=',err_p
!if(i==5623)write(*,*)'delta=',delta(:)
!if(i==5623)write(*,*)'coef=',max_hg,sqrt(delta(1)**2+delta(2)**2),max_hg/sqrt(delta(1)**2+delta(2)**2)
!if(i==5623)write(*,*)'xproject=',xlocan(:)+delta(:)*max_hg/sqrt(delta(1)**2+delta(2)**2)
!if(i==5623)write(*,*)'mmmmmmmmmmmmmmmmmmmmmmm'
          nit=nit+1
       end do
       if(err_p.lt.1.0e-6) then
           call get_indicator_derivative_2D(xlocan,x_inter,x_center,hg,&
                           I_fluid_center,corr_Ip,II,dI,ddI,norm_a,curv_a)
           dis=sqrt((x_ini(1)-xlocan(1))**2+(x_ini(2)-xlocan(2))**2)
           if(I_fluid(i).lt.0.5) then
             curv_nn_temp(i)=1.0/(1.0/(-curv_a)+dis)
           else
             curv_nn_temp(i)=1.0/(1.0/(-curv_a)-dis)
           end if
       else
           write(*,*)'i=',i,' not project',err_p,normx(1:nsd)
       end if





      elseif(nsd==3) then
              call get_indicator_derivative_3D(x(:,i),x_inter,x_center,hg,I_fluid_center,corr_Ip, &
                                        II,dI,ddI,norm_a,curv_a)
      end if
!       den_f=den_liq+(den_inter-den_liq)*I_fluid(i)
!       sur_fluid_temp(1:nsd,i)=dI(1:nsd)!*sur_tension*(-curv_a)*den_f/den_p
!  sur_fluid_temp(1:nsd,i)=dI(1:nsd)*sur_tension*den_f/den_p
!       curv_nn_temp(i)=-curv_a
     end if
  end do
  call mpi_barrier(mpi_comm_world,ierror)
!if(flag==1) then
!  call mpi_allreduce(sur_fluid_temp(1,1),sur_fluid(1,1),nsd*nn,mpi_double_precision, &
!                mpi_sum,mpi_comm_world,ierror)
!else
  call mpi_allreduce(curv_nn_temp(1),curv_nn(1),nn,mpi_double_precision, &
                mpi_sum,mpi_comm_world,ierror)
!end if


  call mpi_barrier(mpi_comm_world,ierror)

end subroutine get_curvature

