c     
c     calculation cauchy stress
c     
      subroutine r_scauchy(det,todet,xto,lx,ly,ne)
      implicit real*8 (a-h,o-z) 
      include 'r_common'
      dimension xto(2,2),ssb(3,3),ss(3,3),tt(3,3)
c
c      ssb(1,1)=btos(1)
c      ssb(2,2)=btos(2)
c      ssb(1,3)=0.0d0
c      ssb(3,1)=0.0d0
c      ssb(2,3)=0.0d0
c      ssb(3,2)=0.0d0
c      ssb(1,2)=btos(3)
c      ssb(2,1)=ss(1,2)
c      ssb(3,3)=btos(4)
c
      ss(1,1)=tos(1)
      ss(2,2)=tos(2)
      ss(1,3)=0.0d0
      ss(3,1)=0.0d0
      ss(2,3)=0.0d0
      ss(3,2)=0.0d0
      ss(1,2)=tos(3)
      ss(2,1)=ss(1,2)
      ss(3,3)=tos(4)
c
      do 20 i=1,4
         cstr(i,ne,lx,ly)=0.0d0
   20 continue
      do 12 i=1,2
         do 13 j=1,2
            tt(i,j)=xto(i,j)
   13    continue
   12 continue
      tt(1,3)=0.0d0
      tt(2,3)=0.0d0
      tt(3,2)=0.0d0
      tt(3,1)=0.0d0
      tt(3,3)=1.0d0
c
c      b1=0.0d0
c      b2=0.0d0
c      b3=0.0d0
c      b4=0.0d0
c
      do 10 m=1,3
         do 11 n=1,3
            cstr(1,ne,lx,ly)=cstr(1,ne,lx,ly)+
     $           todet/det*tt(1,m)*ss(m,n)*tt(1,n)
            cstr(2,ne,lx,ly)=cstr(2,ne,lx,ly)+
     $           todet/det*tt(2,m)*ss(m,n)*tt(2,n)
            cstr(3,ne,lx,ly)=cstr(3,ne,lx,ly)+
     $           todet/det*tt(1,m)*ss(m,n)*tt(2,n)
            cstr(4,ne,lx,ly)=cstr(4,ne,lx,ly)+
     $           todet/det*tt(3,m)*ss(m,n)*tt(3,n)
c
c            b1=b1+todet/det*tt(1,m)*ssb(m,n)*tt(1,n)
c            b2=b2+todet/det*tt(2,m)*ssb(m,n)*tt(2,n)
c            b3=b3+todet/det*tt(1,m)*ssb(m,n)*tt(2,n)
c            b4=b4+todet/det*tt(3,m)*ssb(m,n)*tt(3,n)
 11      continue
 10   continue
c
      return
      end




