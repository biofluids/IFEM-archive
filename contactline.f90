!used to construct contact line

subroutine contactline(x_inter,x,norm_inter)

  use fluid_variables, only:nsd,nn
  use interface_variables

  real(8) x_inter(nsd,maxmatrix),x(nsd,nn),norm_inter(nsd,maxmatrix)

  real(8) x_inter_temp(nsd,maxmatrix),x_con(nsd,maxmatrix),norm_con(nsd,maxmatrix)
  integer i,j,icount,jcount,nn_inter_temp,nn_con
  real(8) x_fix(nsd,2),norm_fix(nsd,2),norm_wall(nsd,2),x_wall(nsd,2)
  real(8) thelta,pi
  real(8) y_limit
  real(8) yl,yr
  integer nn_con_regen
  real(8) x_con_regen(nsd,maxmatrix)

  nn_con_regen=0
  pi=3.1415926
  y_limit=-0.5+2.5*max_hg
!  thelta=pi/2.0!+5*pi/180.0
!  norm_wall(1,1)=-sin(thelta)
!  norm_wall(2,1)=-cos(thelta)
!  norm_wall(1,2)=sin(thelta)
!  norm_wall(2,2)=-cos(thelta)
  nn_inter_temp=0
  nn_con=0
!  x_wall(2,:)=-0.5
  do i=1,nn_inter
     if(abs(x_inter(2,i)).lt.abs(y_limit)) then
       nn_inter_temp=nn_inter_temp+1
       x_inter_temp(:,nn_inter_temp)=x_inter(:,i)
     else
       nn_con=nn_con+1
       x_con(:,nn_con)=x_inter(:,i)
       norm_con(:,nn_con)=norm_inter(:,i)
     end if
  end do
if(minval(x_inter(2,1:nn_inter)).lt.-0.5+2.0*max_hg) then
  y_limit=-0.5+2.5*max_hg
  thelta=pi/2.0+45*pi/180.0
  norm_wall(1,1)=-sin(thelta)
  norm_wall(2,1)=-cos(thelta)
  norm_wall(1,2)=sin(thelta)
  norm_wall(2,2)=-cos(thelta)
  x_wall(2,:)=-0.5

  x_fix(:,1)=x_con(:,1)
  x_fix(:,2)=x_con(:,1)
  norm_fix(:,1)=norm_con(:,1)
  norm_fix(:,2)=norm_con(:,1)
  yr=999.0
  yl=999.0
  do i=1,nn_con
     if(x_con(1,i).lt.0.0) then
       if(abs(x_con(2,i)-y_limit).lt.yl) then
         x_fix(:,1)=x_con(:,i)
         norm_fix(:,1)=norm_con(:,i)
	 yl=abs(x_con(2,i)-y_limit)
       end if
     end if

     if(x_con(1,i).gt.0.0) then
       if(abs(x_con(2,i)-y_limit).lt.yr)then
         x_fix(:,2)=x_con(:,i)
         norm_fix(:,2)=norm_con(:,i)
	 yr=abs(x_con(2,i)-y_limit)
       end if
     end if
  end do
  
  call solve_contact(x_fix(:,1),x_wall(:,1),norm_fix(:,1),norm_wall(:,1),x_con_regen,nn_con_regen)
  call solve_contact(x_fix(:,2),x_wall(:,2),norm_fix(:,2),norm_wall(:,2),x_con_regen,nn_con_regen)
end if



if(maxval(x_inter(2,1:nn_inter)).gt.0.5-2.0*max_hg) then
  y_limit=0.5-2.5*max_hg
  thelta=pi/2.0+45*pi/180.0
  norm_wall(1,1)=-sin(thelta)
  norm_wall(2,1)=cos(thelta)
  norm_wall(1,2)=sin(thelta)
  norm_wall(2,2)=cos(thelta)
  x_wall(2,:)=0.5

  x_fix(:,1)=x_con(:,1)
  x_fix(:,2)=x_con(:,1)
  norm_fix(:,1)=norm_con(:,1)
  norm_fix(:,2)=norm_con(:,1)
  yr=999.0
  yl=999.0
  do i=1,nn_con
     if(x_con(1,i).lt.0.0) then
        if(abs(x_con(2,i)-y_limit).lt.yl) then
          x_fix(:,1)=x_con(:,i)
          norm_fix(:,1)=norm_con(:,i)
          yl=abs(x_con(2,i)-y_limit)
        end if
     end if

     if(x_con(1,i).gt.0.0) then
         if(abs(x_con(2,i)-y_limit).lt.yr)then
            x_fix(:,2)=x_con(:,i)
            norm_fix(:,2)=norm_con(:,i)
            yr=abs(x_con(2,i)-y_limit)
          end if
      end if
    end do
    call solve_contact(x_fix(:,1),x_wall(:,1),norm_fix(:,1),norm_wall(:,1),x_con_regen,nn_con_regen)
    call solve_contact(x_fix(:,2),x_wall(:,2),norm_fix(:,2),norm_wall(:,2),x_con_regen,nn_con_regen)
end if  





  nn_inter=nn_inter_temp
  do i=1,nn_inter_temp
     x_inter(:,i)=x_inter_temp(:,i)
  end do
  do i=1,nn_con_regen
     nn_inter=nn_inter+1
     x_inter(:,nn_inter)=x_con_regen(:,i)
  end do

end subroutine contactline
