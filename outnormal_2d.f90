subroutine outnormal_2d(x_solid,edge_el,normal,nsd_solid,nen_solid,len)
! calculate the outward normal for the edge on the solid interface
integer nsd_solid ! dimension of solid
integer nen_solid ! number of nodes in solid element
real(8) x_solid(nen_solid,nsd_solid) ! coordinates of solid node in the element
integer edge_el(nen_solid) ! flag matrix that deside if the node is on the interface
real(8) normal(nsd_solid) ! output for outward normal
real(8) len ! length of the edge on the boundary
real(8) a(nsd_solid)
real(8) b(nsd_solid)
real(8) c(nsd_solid)
real(8) d(nsd_solid) ! check notes for definition 
real(8) tmp(2,nsd_solid)
integer i
integer j
real(8) f(2,2)
real(8) invf(2,2)
real(8) p(2)
real(8) det
real(8) l
j=1
do i=1,nen_solid
	if (edge_el(i) .eq. 1) then ! nodes on the edge
		tmp(j,1:nsd_solid)=x_solid(i,1:nsd_solid)
		j=j+1
	end if

	if (edge_el(i) .eq. 0) then ! nodes not on the edge
		c(1:nsd_solid)=x_solid(i,1:nsd_solid)
	end if
end do
a(1:nsd_solid)=tmp(1,1:nsd_solid)
b(1:nsd_solid)=tmp(2,1:nsd_solid)


f(1,1)=b(2)-a(2)
f(1,2)=-(b(1)-a(1))
f(2,1)=b(1)-a(1)
f(2,2)=b(2)-a(2)

p(1)=a(1)*(b(2)-a(2))-a(2)*(b(1)-a(1))
p(2)=c(1)*(b(1)-a(1))+c(2)*(b(2)-a(2))

det=f(1,1)*f(2,2)-f(1,2)*f(2,1)

invf(1,1)=f(2,2)/det
invf(2,2)=f(1,1)/det
invf(1,2)=-f(1,2)/det
invf(2,1)=-f(2,1)/det

d(:)=0.0d0

do i=1,nsd_solid
	do j=1,nsd_solid
		d(i)=d(i)+invf(i,j)*p(j)
	end do
end do

normal(1:nsd_solid)=d(1:nsd_solid)-c(1:nsd_solid)
l=normal(1)**2+normal(2)**2
l=sqrt(l)
normal(1:nsd_solid)=normal(1:nsd_solid)/l

len=(a(1)-b(1))**2+(a(2)-b(2))**2
len=sqrt(len)

!write(*,*) 'normal inward', normal(:)
return
end subroutine outnormal_2d
