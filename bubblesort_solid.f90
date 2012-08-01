subroutine bubblesort_solid(a,n,m,row)
! preform bubble sort for send_address matrix 
! matirx a to be sorted
! n --- lengith of a
! m ---width of a
! row --- which row to be sort
use mpi_variables

implicit none

integer a(n,m)
integer n
integer m
integer row
integer b


integer i
integer j
integer count
integer tmp(1,m)


do i=n,1,-1
	do j=1,i-1
	if (a(j,row) .lt. a(j+1,row)) then
		tmp(1,:) = a(j,:)
		a(j,:) = a(j+1,:)
		a(j+1,:) = tmp(1,:)
	end if
	end do
end do

b=a(1,row)
countrow_solid=1

do i=1,n
	if (a(i,row) .eq. b) then
		count=count+1
	else
	   countrow_solid=countrow_solid+1
	   b=a(i,row)
	end if
end do

if (myid == 0) write(*,*) 'countow', countrow

allocate(sub_address_solid(countrow_solid,2))
! sub_address_solid(proc id, number of nodes share by this proc)

b=a(1,row)
count=0
countrow_solid=1
do i=1,n
        if (a(i,row) .eq. b) then
	count=count+1
	else
	sub_address_solid(countrow_solid,1)=b
	sub_address_solid(countrow_solid,2)=count
	countrow_solid=countrow_solid+1
	b=a(i,row)
	count=1
	end if
end do

sub_address_solid(countrow_solid,1)=b
sub_address_solid(countrow_solid,2)=count

!if (myid == 0) then
!do i=1,countrow
!write(*,*) 'sub_address', sub_address(i,:)
!end do
!end if
!continue

end
