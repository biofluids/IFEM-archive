      subroutine slnpd(n,ifail,a,b,m)
c-------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
c.... solve system of equations by gauss elimination process and
c.... return result in vector  b
c
      dimension a(m,m),b(m)
c
      n1=n-1

      do 100 k=1,n1
        k1=k+1
        cc=a(k,k)
        if(abs(cc)) 1,1,3

1       continue

        do 7 j=k1,n

          if(abs(a(j,k))) 7,7,5

5         continue

          do 6 l=k,n
            cc=a(k,l)
            a(k,l)=a(j,l)
            a(j,l)=cc
6         continue

          cc=b(k)
          b(k)=b(j)
          b(j)=cc
          cc=a(k,k)

          go to 3

7       continue

c       write(*,2) k
2       format(/,'* * * singularity in row',i5,' * * * ')

        ifail=1

        go to 300

3       continue

        do 4 j=k1,n
          a(k,j)=a(k,j)/cc
4       continue

        b(k)=b(k)/cc

        do 10 i=k1,n
          cc=a(i,k)
          do 9 j=k1,n
            a(i,j)=a(i,j)-cc*a(k,j)
9         continue
          b(i)=b(i)-cc*b(k)
10      continue

100   continue

      if(abs(a(n,n))) 8,8,101

    8 continue

c     write(*,2) k

      ifail=1

      go to 300

  101 continue

      b(n)=b(n)/a(n,n)
      do 200 l=1,n1
        k=n-l
        k1=k+1
        do 200 j=k1,n
          b(k)=b(k)-a(k,j)*b(j)
200   continue

      ifail=0

300   continue

      return
      end

