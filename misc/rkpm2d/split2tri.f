c
	subroutine split2tri(nelem,lnods,xm)
c
c*** this subroutine subdivides each 4-node 
c    quadrilateral element into
c    two triangular elements
c
c
        implicit none
        include 'parameter.h'
c        
        integer mnode,ndime
        parameter (mnode=4,ndime=2)
        
        integer nelem
        real*8 xm(2,maxNumnp)
        integer lnods(maxNode,maxElem)
        
           !** local vars
	real*8 corde(2,mnode)
        integer ltemp(mnode)
        real*8 diag1,diag2,difer,epis_diag
        integer kount,notal,index
        integer ielem,inode,idime
        
        data epis_diag/1.e-6/

	write(*,*) '*** Get into sub SPLIT ***'
	call checklim(nelem*2,maxElem,'nelem')

	kount=0
	do ielem=1,nelem
   	   notal=nelem+ielem
c	   matno(notal)=matno(ielem)
        
	   do inode=1,mnode
	      lnods(inode,notal)=lnods(inode,ielem)
	   enddo
        enddo

	do ielem=1,nelem
	   notal=nelem+ielem
	   do inode=1,mnode
	      index=lnods(inode,notal)
	      ltemp(inode)=index
 
	      do idime=1,ndime
   	         corde(idime,inode)=xm(idime,index)
              enddo
           enddo

	   diag1=sqrt( (corde(1,1)-corde(1,3))**2
     &                +(corde(2,1)-corde(2,3))**2 )
	   diag2=sqrt( (corde(1,2)-corde(1,4))**2
     &                +(corde(2,2)-corde(2,4))**2 )
c
c**** divide across the shorter diagonal
c
	   difer=diag1-diag2
	   if( (difer/max(diag1,diag2)).gt.-epis_diag) then
	      kount=kount+1
	      lnods(1,kount)=ltemp(1)
	      lnods(2,kount)=ltemp(2)
	      lnods(3,kount)=ltemp(4)
c	      matno(kount)=matno(notal)
	      kount=kount+1
	      lnods(1,kount)=ltemp(2)
	      lnods(2,kount)=ltemp(3)
	      lnods(3,kount)=ltemp(4)
c	      matno(kount)=matno(notal)
           else
	      kount=kount+1
	      lnods(1,kount)=ltemp(1)
	      lnods(2,kount)=ltemp(2)
	      lnods(3,kount)=ltemp(3)
c	      matno(kount)=matno(notal)
	      kount=kount+1
	      lnods(1,kount)=ltemp(1)
	      lnods(2,kount)=ltemp(3)
	      lnods(3,kount)=ltemp(4)
c	      matno(kount)=matno(notal)
           endif
        enddo
c        
	nelem=2*nelem
c
        return        
	end	! endp split2tri
c
