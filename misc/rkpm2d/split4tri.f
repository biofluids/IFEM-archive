
	subroutine split4tri(nelem,numnp,lnods,xm)
c*** this subroutine subdivides each 4-node quadrilateral element into
c    four triangular elements
c

c** input parameters:
c   nelem
c   numnp
c   lnods
c   xm
c** output parameters:
c   nelem : the new one
c   numnp : the new one
c   lnods : the new one
c   xm    : the new one

        implicit none
        include 'parameter.h'
        
        integer mnode,ndim
        parameter (mnode=4,ndim=2)
        
        integer nelem,numnp
        real*8 xm(2,maxNumnp)
        integer lnods(maxNode,maxElem)
        
           !** local vars
        integer ltemp(5)
        integer newNodeInd
        integer kount,notal,index
        integer ielem,inode,idim
        
	write(*,*) '*** Get into sub SPLIT4 ***'
	call checklim(nelem*4,maxElem,'nelem')

	kount=0
	do ielem=1,nelem
   	   notal=nelem*3+ielem
c	   matno(notal)=matno(ielem)
        
	   do inode=1,mnode
	      lnods(inode,notal)=lnods(inode,ielem)
	   enddo
        enddo

        newNodeInd=numnp
	do ielem=1,nelem
	   notal=nelem*3+ielem
	   do inode=1,mnode
	      index=lnods(inode,notal)
	      ltemp(inode)=index
            enddo

           newNodeInd = newNodeInd + 1
      	   call checklim(newNodeInd,maxNumnp,'Sub split4:newNodeInd')
           
           ltemp(5) = newNodeInd
           
           do idim=1,ndim
              xm(idim,newNodeInd) = 0.0
              do inode=1,mnode
                 xm(idim,newNodeInd)=xm(idim,newNodeInd)
     &                              +xm(idim,ltemp(inode))
	      enddo
              xm(idim,newNodeInd)=xm(idim,newNodeInd)/real(mnode)
           enddo

	   kount=kount+1
	   lnods(1,kount)=ltemp(1)
	   lnods(2,kount)=ltemp(2)
	   lnods(3,kount)=ltemp(5)
c	   matno(kount)=matno(notal)

	   kount=kount+1
	   lnods(1,kount)=ltemp(2)
	   lnods(2,kount)=ltemp(3)
	   lnods(3,kount)=ltemp(5)
c	   matno(kount)=matno(notal)

	   kount=kount+1
	   lnods(1,kount)=ltemp(3)
	   lnods(2,kount)=ltemp(4)
	   lnods(3,kount)=ltemp(5)
c	   matno(kount)=matno(notal)

	   kount=kount+1
	   lnods(1,kount)=ltemp(4)
	   lnods(2,kount)=ltemp(1)
	   lnods(3,kount)=ltemp(5)
c	   matno(kount)=matno(notal)

        enddo
        
	nelem=4*nelem
        numnp=newNodeInd
        
	end	! endp split4tri
