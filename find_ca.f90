!================
!find the capillary number
!================

subroutine find_ca(x_inter,x,vel_fluid,vol_nn,ca,norm_p,thelta,anglet,flag_ca)

  use fluid_variables, only:nsd,ne,nn,nen
  use interface_variables
  use mpi_variables
  include 'mpif.h'

  real(8) x_inter(nsd)
  real(8) x(nsd,nn)
  real(8) vel_fluid(nsd,nn)
  real(8) vol_nn(nn)
  real(8) ca,V_cl(nsd),norm_p(nsd),thelta
  real(8) cap
  integer i,j,icount,jcount
  real(8) dx(nsd),Sp,hs,temp

  real(8) M(nsd+1,nsd+1),B(nsd+1),P(nsd+1)
  real(8) vec(nsd+1)
  integer IP(nsd+1),INFO
  real(8) anglet
  integer flag_ca

  M(:,:)=0.0
  B(:)=0.0
  P(:)=0.0
  P(1)=1.0
  vec(1)=1.0
  V_cl(:)=0.0
  flag_ca=1
  do j=1,nn
     dx(:)=abs(x_inter(:)-x(:,j))
     call B_Spline1(dx,hsp,nsd,Sp,INFO)
     if(INFO==1) then
     vec(2:nsd+1)=x_inter(:)-x(:,j)
     do icount=1,nsd+1
           do jcount=1,nsd+1
                 M(icount,jcount)=M(icount,jcount)+vec(icount)*vec(jcount)*Sp/(hsp**nsd)*vol_nn(j)
           end do
     end do
     end if
  end do
  call DGESV(nsd+1,1,M,nsd+1,IP,P,nsd+1,INFO)
  B(:)=P(:)
  do j=1,nn
     dx(:)=abs(x_inter(:)-x(:,j))
     call B_Spline1(dx,hsp,nsd,Sp,INFO)
     if(INFO==1) then
     vec(2:nsd+1)=x_inter(:)-x(:,j)
     temp=0.0
     do icount=1,nsd+1
        temp=temp+vec(icount)*B(icount)
     end do
     V_cl(:)=V_cl(:)+vel_fluid(:,j)*temp*Sp/(hsp**nsd)*vol_nn(j)
     end if
  end do

!  temp=0.0
!  do icount=1,nsd
!     temp=temp+V_cl(icount)**2
!  end do
!  temp=sqrt(temp)
!  temp=abs(V_cl(1))
!  ca=vis_inter*temp/sur_tension
!if(myid==0)write(*,*)'ca=',ca,'V_cl=',temp
!if(myid==0)write(*,*)'x-inter=',x_inter(:)
  temp=0.0
  do icount=1,nsd
     temp=temp+V_cl(icount)*norm_p(icount)
  end do

  ca=vis_inter*abs(temp)/sur_tension
   if((temp.ge.0.0).and.(anglet.lt.(static_angle+ad_re_angle))) then
     flag_ca=0
!     write(*,*)'advancing but angle less then theta ad',myid
!     thelta=3.14159/180.0*static_angle
     goto 999
   end if
   if((temp.le.0.0).and.(anglet.gt.(static_angle-ad_re_angle))) then
     flag_ca=0
!     write(*,*)'receding but angle greater than theta_re',myid
!    thelta=3.14159/180.0*static_angle 
    goto 999
   end if

   if(temp.ge.0.0) then
    cap=Hoff_ad+ca
   else
    cap=Hoff_re+ca
   end if
   thelta=acos(1-2*tanh(5.16*(cap/(1+1.31*(cap**0.99)))**0.706))


   if(temp.ge.0.0) then
     thelta=thelta
   else
     thelta=3.14159/180.0*(static_angle-ad_re_angle)*2.0-thelta
   end if

if(abs(thelta/3.14159*180.0-anglet).gt.10.0) then
  thelta=(thelta/3.14159*180.0-anglet)/abs(thelta/3.14159*180.0-anglet)*10.0+anglet
  thelta=thelta/180.0*3.14159
end if


write(*,*)'ca=',ca,'thelta=',thelta/3.14159*180.0,'adv?=',temp,myid
999 continue
end subroutine find_ca









