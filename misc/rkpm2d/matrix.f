      subroutine scalemul_matx(source_matx,scale,target_matx,m,n)
      
      implicit none
      integer m,n
      real*8  source_matx(m,n),target_matx(m,n)
      real*8  scale
      
      integer i,j
c      
      do i=1,m
         do j=1,n
            target_matx(i,j)=source_matx(i,j)*scale
         enddo
      enddo
c      
      end	!ends scalemul_matx
c      
c---------------------------------------------------------      
c      
      subroutine add_matx(a,b,c,m,n)
      
      implicit none
      integer m,n
      real*8 a(m,n),b(m,n),c(m,n)
      
      integer i,j
      
      do i=1,m
         do j=1,n
            c(i,j)=a(i,j)+b(i,j)
         enddo
      enddo
      
      end	!ends add_matx
c      
c---------------------------------------------------------      
      subroutine minus_matx(a,b,c,m,n)
      
      implicit none
      integer m,n
      real*8 a(m,n),b(m,n),c(m,n)
      
      integer i,j
      
      do i=1,m
         do j=1,n
            c(i,j)=a(i,j)-b(i,j)
         enddo
      enddo
      
      end	!ends minus_matx
      
c---------------------------------------------------------      
     
      subroutine invers_matx33(a,b,det)
c*** inverse a 3x3 matx

      implicit none
      
      integer n
      parameter (n=3)
      real*8 a(n,n),b(n,n)
      
      real*8 det,invdet
      real*8 epis_det
      real*8 t1,t2,t3
      data epis_det/1.e-8/
      
      t1=a(1,1)*( a(2,2)*a(3,3)-a(2,3)*a(3,2) )
      t2=-a(1,2)*( a(2,1)*a(3,3)-a(3,1)*a(2,3) )
      t3=a(1,3)*( a(2,1)*a(3,2)-a(3,1)*a(2,2) )
      det=t1+t2+t3
c      write(*,*) 'det=',det
      
      if (abs(det).lt.epis_det) then
         write(*,*) 'singular matx in invers_matx33'
         write(*,*) 't1=',t1,' t2=',t2,' t3=',t3
         stop
      endif
      
      invdet=1./det
      b(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))*invdet
      b(1,2)=-(a(1,2)*a(3,3)-a(3,2)*a(1,3))*invdet
      b(1,3)=(a(1,2)*a(2,3)-a(2,2)*a(1,3))*invdet
      
      b(2,1)=-(a(2,1)*a(3,3)-a(3,1)*a(2,3))*invdet
      b(2,2)=(a(1,1)*a(3,3)-a(3,1)*a(1,3))*invdet
      b(2,3)=-(a(1,1)*a(2,3)-a(2,1)*a(1,3))*invdet
      
      b(3,1)=(a(2,1)*a(3,2)-a(3,1)*a(2,2))*invdet
      b(3,2)=-(a(1,1)*a(3,2)-a(3,1)*a(1,2))*invdet
      b(3,3)=(a(1,1)*a(2,2)-a(2,1)*a(1,2))*invdet
      
      end	!ends  invers_matx
           
c---------------------------------------------------------     
      subroutine mul_matx(a_matx,b_matx,c_matx,m,n,L)
      
      implicit none
      integer m,n,L
      real*8 a_matx(m,n),b_matx(n,L),c_matx(m,L)
      
      integer i,j,k
      
      do i=1,m
         do j=1,L
            c_matx(i,j) = 0.0
            
            do k=1,n
               c_matx(i,j) = c_matx(i,j)
     &         + a_matx(i,k)*b_matx(k,j)
            enddo
         enddo
      enddo
      
      end	!ends mul_matx
      
c--------------------------------------------------------- 
 
      subroutine trans_matx(a,b,m,n)
      
      implicit none
      integer m,n
      real*8 a(m,n),b(n,m)
      
      integer i,j
      
      
      do i=1,m
         do j=1,n
            b(j,i)=a(i,j)
         enddo
      enddo
      
      end	! ends trans_matx
      
