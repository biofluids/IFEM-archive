subroutine normal2d_gen(x_el,nnode1,nnode2,normal,nsd,nen,len)
! calculate the outward normal for the edge on the solid interface
integer nsd ! dimension of solid
integer nen ! number of nodes in solid element
real(8) x_el(nsd,nen) ! coordinates of solid node in the element
integer nnode1,nnode2
real(8) normal(nsd) ! output for outward normal
real(8) len ! length of the edge on the boundary
real(8) a(nsd)
real(8) b(nsd)
real(8) c(nsd)
real(8) d(nsd) ! check notes for definition 
integer i
integer j
real(8) f(2,2)
real(8) invf(2,2)
real(8) p(2)
real(8) det
real(8) l

do i=1,nen
    if ((i .ne. nnode1) .and. (i .ne. nnode2)) then ! nodes not on the edge
        c(1:nsd)=x_el(1:nsd,i)
    endif
enddo

a(1:nsd)=x_el(1:nsd,nnode1)
b(1:nsd)=x_el(1:nsd,nnode2)


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

do i=1,nsd
    do j=1,nsd
        d(i)=d(i)+invf(i,j)*p(j)
    enddo
enddo

normal(1:nsd)=d(1:nsd)-c(1:nsd)
l=normal(1)**2+normal(2)**2
l=sqrt(l)
normal(1:nsd)=normal(1:nsd)/l

len=(a(1)-b(1))**2+(a(2)-b(2))**2
len=sqrt(len)

!write(*,*) 'normal inward', normal(:)
return
end subroutine normal2d_gen
