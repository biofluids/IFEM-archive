subroutine outnormal_hex(x,nsd,p,outnorm,area)
! calculate the outward normal for the surface  on the solid interface 
! and the area of the surface for tet element cases only 
real(8) x(nsd,4)
real(8) p(nsd)
integer nsd
real(8) outnorm(nsd)
real(8) area
!---------------------------
real(8) ab(nsd)
real(8) ac(nsd)
real(8) ad(nsd)

real(8) ap(nsd)
real(8) dot_pro
real(8) temp(nsd)
real(8) temp_area

ab(:)=x(:,2)-x(:,1)
ac(:)=x(:,3)-x(:,1)
ad(:)=x(:,4)-x(:,1)

ap(:)=p(:)-x(:,1)

outnorm(1)=ab(2)*ac(3)-ab(3)*ac(2)
outnorm(2)=-ab(1)*ac(3)+ab(3)*ac(1)
outnorm(3)=ab(1)*ac(2)-ab(2)*ac(1)

area=outnorm(1)**2+outnorm(2)**2+outnorm(3)**2
area=0.5*sqrt(area) ! area of abc

dot_pro=ap(1)*outnorm(1)+ap(2)*outnorm(2)+ap(3)*outnorm(3)

if (dot_pro > 0) then
outnorm(:)=-outnorm(:)
end if

outnorm(:)=outnorm(:)/area*0.5

temp(1)=ac(2)*ad(3)-ac(3)*ad(2)
temp(2)=-ac(1)*ad(3)+ac(3)*ad(1)
temp(3)=ac(1)*ad(2)-ac(2)*ad(1)

temp_area=temp(1)**2+temp(2)**2+temp(3)**2
temp_area=0.5*sqrt(temp_area) ! area of acd

area=area+temp_area

continue

end

