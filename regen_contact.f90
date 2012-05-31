

subroutine regen_contact(xc,x,x_inter,x_center,hg,I_fluid_center,corr_Ip,nn_org_regen,x_org_regen,gridy,norm_e,tang,&
			nn_con_ele,con_ele,ne_regen,ele_regen,ien,norm_p,nn_bottom,x_bottom)

  use fluid_variables,only:nn,nsd,ne,nen
  use interface_variables, only:maxmatrix,static_angle,ad_re_angle,max_hg
  use mpi_variables
  real(8) xlocan(nsd),x_inter(nsd,maxmatrix),x_center(nsd,ne),hg(ne),xc(nsd),x(nsd,nn)
  real(8) I_fluid_center(ne),corr_Ip(maxmatrix),xlocan_temp(nsd)
  integer nn_org_regen,nn_bottom
  real(8) x_org_regen(nsd,60),gridy,norm_p(nsd),x_bottom(nsd,60)
  integer ien(nen,ne)
  integer nn_con_ele,con_ele(nn_con_ele),ne_regen,ele_regen(100)

  integer i,j,icount,jcount,inl,nit,max_n
  real(8) err_p,delta(nsd)
  real(8) II,dI(nsd),ddI(3*(nsd-1)),norm_e(nsd),curv_p,norm_g(nsd),tang(nsd)
 
  real(8) interval,temp(nsd),distance
  real(8) vector_static(nsd),pi,thelta
  integer finf,maxconn,flag
  real(8) x_fix(nsd)

write(*,*)'myid in regen_contact=',myid
  maxconn=30
  xlocan(:)=xc(:)
  pi=3.1415926
  max_n=40
  interval=max_hg/real(max_n)

  x_org_regen(:,:)=0.0
  nn_org_regen=0
  nn_bottom=0
  icount=0
  temp(:)=xc(:)
  x_fix(:)=xc(:)
  flag=0
  do i=1,50
     nit=1
     err_p=999.0
     delta(:)=0.0
     if(flag==0) then
        nn_bottom=nn_bottom+1
        x_bottom(:,nn_bottom)=temp(:)
     end if

     xlocan(:)=temp(:)+interval*icount*norm_e(:)
     xlocan_temp(:)=xlocan(:)

     do while((nit.lt.5).and.(err_p.gt.1.0e-8))
	xlocan(:)=xlocan(:)+delta(:)
        call get_indicator_derivative_2D_1st(xlocan,x_inter,x_center,hg,&
             I_fluid_center,corr_Ip,II,dI,ddI,norm_g,curv_p)
!	delta=(0.5-II)/dI(1)
        delta(1)=(0.5-II)*tang(1)/(tang(1)*dI(1)+tang(2)*dI(2))
        delta(2)=(0.5-II)*tang(2)/(tang(1)*dI(1)+tang(2)*dI(2))
	err_p=abs(0.5-II)
	nit=nit+1
	if(abs(tang(1)*dI(1)+tang(2)*dI(2)).lt.1.0e-5) then
	   err_p=999.0
	   nit=999
	end if
     end do
     if(err_p.gt.1.0e-6) then
        icount=icount+1
        goto 200
      end if
       distance=(xlocan_temp(1)-xlocan(1))**2+(xlocan_temp(2)-xlocan(2))**2
       distance=sqrt(distance)
       if(distance.gt.1.0*max_hg) then
         icount=icount+1
          goto 200
        end if
     finf=0
     call getinf_el_3d_con(finf,xlocan,x,nn,nsd,ne,nen,ien,maxconn,nn_con_ele,con_ele)
     if( abs(((xlocan(1)-x_fix(1))*norm_e(1)+(xlocan(2)-x_fix(2))*norm_e(2))).gt.gridy-interval) then
        if(finf==0)goto 999

     flag=flag+1
     end if
     if(flag==1) then
       xc(:)=xlocan(:)
       norm_p(:)=norm_g(:)
     end if
213 continue

     icount=1
     temp(:)=xlocan(:)
    if(flag.gt.0) then
       nn_org_regen=nn_org_regen+1
       x_org_regen(:,nn_org_regen)=xlocan(:)
     end if
if(finf.ne.0) then
     if(ne_regen==0) then
	ne_regen=1
        ele_regen(ne_regen)=finf
     else if(finf.ne.ele_regen(ne_regen)) then
        ne_regen=ne_regen+1
        ele_regen(ne_regen)=finf
     end if
end if
200 continue
  end do
write(*,*)'finishi regen myid=',myid
999 continue
end subroutine regen_contact



