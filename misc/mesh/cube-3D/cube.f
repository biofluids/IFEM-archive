      program makecube
      parameter (nsd=3,nen=8,neface=6)
      integer i,m,ie,n,ien(nen,30000),ien2(nen,30000),rng(neface,30000)
      integer numelex,numeley,numelez,nn,l
      integer top,bot,left,right,front,back,node
      character*10 letter
      real*8 x(nsd,30000),alpha,x3
      real*8 lengthx,lengthy,lengthz,length,dx,dy,dz
	real*8 xshift,yshift,zshift
      open(unit=10, file='mxyz.in')
      open(unit=11, file='mien.in')
      open(unit=15, file='mrng.in')
      open(unit=31, file='tout.dat')
	open(unit=32, file='mesh.info')

!  for ifem fluids only
	length=1.0d0
      lengthx=5.0d0*length
	lengthy=0.25d0*length
      lengthz=2.0d0*length

      lengthx=5.0d0*length
	lengthy=1.0d0*length
      lengthz=1.0d0*length

	xshift=0.0d0
	yshift=0.0d0
	zshift=0.0d0
!	yshift=-lengthy/2.0d0
!	zshift=6.0d0

	dx=0.125d0
!	dx=0.125d0
	dy=dx
	dz=dx
! for rubber 3-D testing case only
c	length=1.0d0
c	lengthx=2.0d0*length
c	lengthy=2.0d0*length
c	lengthz=2.0d0*length
c	xshift=0.0d0
c	yshift=0.0d0
c	zshift=0.0d0
c	dx=1.0d0
c	dy=dx
c	dz=dx
!--------------------------------------
	numelex=lengthx/dx
	numeley=lengthy/dy
      numelez=lengthz/dz	

	nn=(numelex+1)*(numeley+1)*(numelez+1)
      ne=numelex*numeley*numelez
      node=0

	do k=1,numelez+1  !the order CANNOT be switched
	   do j=1,numeley+1
            do i=1,numelex+1
               node=node+1
               x(1,node)=(i-1)*dx+xshift
               x(2,node)=(j-1)*dy+yshift
               x(3,node)=(k-1)*dz+zshift
            enddo
         enddo
      enddo
 900  format(i5, 4f10.5)
 901 	format(3f10.5)
      do i=1,node
         write(10,901) x(1,i),x(2,i),x(3,i)
         write(31,900) i,x(1,i),x(2,i),x(3,i),0.0
      enddo
      l=0
      do i=1,numelez
         do j=1,numeley
            do k=1,numelex
               l=l+1
               ien(1,l)=l+(j-1)+(i-1)*(numelex+numeley+1) 
               ien(2,l)=ien(1,l)+1
               ien(3,l)=ien(2,l)+(numelex+1)*(numeley+1)
               ien(4,l)=ien(3,l)-1
               ien(5,l)=ien(1,l)+numelex+1
               ien(6,l)=ien(2,l)+numelex+1
               ien(7,l)=ien(3,l)+numelex+1
               ien(8,l)=ien(4,l)+numelex+1
            enddo
         enddo
      enddo

      
      do i=1,ne
         do j=1,neface
            rng(j,i)=0
         enddo
      enddo
      do j=1,numelez
         do i=1,numelex
            bot=i+(j-1)*(numelex*numeley)
            top=numelex*(numeley-1)+bot
            rng(1,bot)=1  !front if y and z are 
            rng(6,top)=6
         enddo
      enddo
      l=0
      do j=1,numeley
         do i=1,numelex
            l=l+1
            front=l
            back=ne+1-front
            rng(2,front)=2
            rng(4,back)=4
         enddo
      enddo
      l=0
      do j=1,numelez
         do i=1,numeley
            l=l+1
            right=numelex*l
            left=right-(numelex-1)
            rng(3,right)=3
            rng(5,left)=5
         enddo
      enddo
     
      do i=1,ne
!	write to "mien.dat"
         write(11,923) ien(1,i),ien(2,i),ien(3,i),ien(4,i),
     +              ien(5,i),ien(6,i),ien(7,i),ien(8,i)
!	write to "tout.dat"
         write(31,921) i, ien(1,i),ien(2,i),ien(3,i),ien(4,i),
     +              ien(5,i),ien(6,i),ien(7,i),ien(8,i)
      enddo
      do i=1,ne
	!	write to 'mrng.dat'
         write(15,922) rng(1,i),rng(2,i),rng(3,i),rng(4,i),
     +              rng(5,i),rng(6,i)
      enddo
      
 921  format(10i7)
 922  format(6i4)
 923	format(8i7)

	write(32,*) 'nn=',node
	write(32,*) 'ne=',ne
	write(32,*) 'rng=',6
	write(32,*) 'lengthx=',lengthx
	write(32,*) 'lengthy=',lengthy
	write(32,*) 'lengthz=',lengthz
	write(32,*) 'numelex=',numelex
	write(32,*) 'numeley=',numeley
	write(32,*) 'numelez=',numelez
	write(32,*) 'dx=',dx
	write(32,*) 'dy=',dy
	write(32,*) 'dz=',dz

	write(32,*) 'rng=1: bottom -> z=',0.0d0+zshift
	write(32,*) 'rng=2: front  -> y=',0.0d0+yshift
	write(32,*) 'rng=3: right  -> x=',lengthx+xshift
	write(32,*) 'rng=4: back   -> y=',lengthy+yshift
	write(32,*) 'rng=5: left   -> x=',0.0d0+xshift
	write(32,*) 'rng=6: top    -> z=',lengthz+zshift

	close(10)
	close(11)
	close(15)
	close(31)
	close(32)

      end

