      program makecube
      parameter (nsd=2,nen=4,neface=4)
      integer i,m,ie,n,ien(nen,10000),ien2(nen,10000),rng(neface,10000)
      integer numelex,numeley,nn,l
      integer top,bot,left,right,node
      character*10 letter
      real x(nsd,10000),alpha,x3
      real lengthx,lengthy,dx,dy
      open(unit=10, file='mxyz.in')
      open(unit=11, file='mien.in')
      open(unit=15, file='mrng.in')
!      open(unit=31, file='tout.dat')
      lengthx=10
      lengthy=2
      numelex=100
      numeley=20
      dx=lengthx/numelex
      dy=lengthy/numeley
      nn=(numelex+1)*(numeley+1)
      ne=numelex*numeley
      node=0
         do j=1,numeley+1
            do i=1,numelex+1
               node=node+1
c               write(*,*) node
               x(1,node)=(i-1)*dx
               x(2,node)=(j-1)*dy
c               write(*,*) node,x(1,node),x(2,node)
            enddo
         enddo
      
 20   format(f7.3,f7.3)
 26   format(f7.3,f7.3)
      do i=1,nn
         write(10,20) x(1,i),x(2,i)
!         write(31,26) x(1,i),x(2,i)
      enddo
      l=0
      do j=1,numeley
      	do k=1,numelex
               l=l+1
               ien(1,l)=l+(j-1) 
               ien(4,l)=ien(1,l)+numelex+1
               ien(3,l)=ien(4,l)+1
               ien(2,l)=ien(1,l)+1
          enddo
       enddo
    

      
      do i=1,ne
         do j=1,neface
            rng(j,i)=0
         enddo
      enddo

       l=0
    
         do i=1,numelex
            l=l+1
            bot=l
            top=ne+1-bot
            rng(1,bot)=1
            rng(3,top)=3
         enddo
     
       l=0
     
         do i=1,numeley
            l=l+1
            right=numelex*l
            left=right-(numelex-1)
            rng(2,right)=2
            rng(4,left)=4
         enddo




     
      do i=1,ne
         write(11,21) ien(1,i),ien(2,i),ien(3,i),ien(4,i)
!         write(31,27) ien(1,i),ien(2,i),ien(3,i),ien(4,i)
      enddo
      do i=1,ne
         write(15,22) rng(1,i),rng(2,i),rng(3,i),rng(4,i)
      enddo
	write(*,*) 'nn=',nn
	write(*,*) 'ne=',ne
      
 21   format(i7,i7,i7,i7)
 22   format(i4,i4,i4,i4)
 27   format(i7,i7,i7,i7)
      end

