subroutine gaussj(a,n,npp,b,m,mp)
  implicit none

  integer :: n,npp,m,mp
  integer,parameter :: nmax = 10000
  real(8) ::  a(npp,npp),b(npp,mp)
  integer,dimension(nmax) :: ipiv,indxr,indxc
  integer :: i,j,k,l,ll
  integer :: irow,icol
  real(8) :: big,dum,pivinv

  do j=1,n
     ipiv(j)=0
  enddo

  do i=1,n
     big=0.0d0
     do j=1,n
        if(ipiv(j).ne.1)then
           do k=1,n
              if (ipiv(k).eq. 0) then
                 if (dabs(a(j,k)) .ge. big)then
                    big=dabs(a(j,k))
                    irow=j
                    icol=k
                 endif
              else if (ipiv(k).gt.1) then
                 write(*,*) 'singular matrix'
                 stop
              endif
           enddo
        endif
     enddo

     ipiv(icol)=ipiv(icol)+1
     
     if (irow.ne.icol) then
        do l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        enddo
        do l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        enddo
     endif

     indxr(i)=irow
     indxc(i)=icol
     if (a(icol,icol) .eq. 0.0d0) then
        write(*,*) 'singular matrix.'
        stop
     endif

     pivinv=1.0d0/a(icol,icol)
     a(icol,icol)=1.0d0

     do l=1,n
        a(icol,l)=a(icol,l)*pivinv
     enddo
     do l=1,m
        b(icol,l)=b(icol,l)*pivinv
     enddo

     do ll=1,n
        if(ll.ne.icol)then
           dum=a(ll,icol)
           a(ll,icol)=0.0d0
           do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
           enddo
           do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
           enddo
        endif
     enddo

  enddo

  do l=n,1,-1
     if(indxr(l).ne.indxc(l))then
        do k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
        enddo
     endif
  enddo
  return
  end

