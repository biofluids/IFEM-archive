c     
c     green-lagrangian strain and derivative
c     
      subroutine r_sstrain(toc,xto,lx,ly,lz,ne)
      implicit real*8 (a-h,o-z)
      include 'r_common'
      dimension xto(3,3),toc(3,3)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dge(i,j,k)      -- i -- strain, j   -- dir,   k -- node
c    ddge(i,j,m,k,n)  -- i -- strain, j,m -- dir, k,n -- node
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,3
         ge(i,ne,lx,ly,lz)=0.5d0*(toc(i,i)-1.0d0)
	enddo

      ge(4,ne,lx,ly,lz)=toc(2,3)
      ge(5,ne,lx,ly,lz)=toc(1,3)
	ge(6,ne,lx,ly,lz)=toc(1,2)
c  for 3-D only
	do i=1,6
	  do k=1,nis
	    do j=1,3
            if (i .eq. 4) then
                dge(i,j,k)=xto(j,2)*bd(3,k)+xto(j,3)*bd(2,k)
	      elseif (i .eq. 5) then
                dge(i,j,k)=xto(j,1)*bd(3,k)+xto(j,3)*bd(1,k)
		  elseif (i .eq. 6) then
                dge(i,j,k)=xto(j,1)*bd(2,k)+xto(j,2)*bd(1,k)
            elseif (i. le. 3) then
                dge(i,j,k)=xto(j,i)*bd(i,k)
            endif

		  do m=1,3
			do n=1,nis
			  if (m.eq.j) then
                  if (i .eq. 4) then
	              ddge(i,j,m,k,n)=bd(2,n)*bd(3,k)+bd(3,n)*bd(2,k)
                  elseif (i .eq.5) then
				  ddge(i,j,m,k,n)=bd(1,n)*bd(3,k)+bd(3,n)*bd(1,k)
                  elseif (i .eq.6) then
				  ddge(i,j,m,k,n)=bd(1,n)*bd(2,k)+bd(2,n)*bd(1,k)
				else
				  ddge(i,j,m,k,n)=bd(i,n)*bd(i,k)
                  endif
			  else
				ddge(i,j,m,k,n)=0.0d0
			  endif
			enddo
		  enddo
		enddo
	  enddo
	enddo

c  for 2-D only(original version)
c      do i=1,3
c         do k=1,nis
c            do j=1,2
c               if (i .eq. 3) then
c                  dge(i,j,k)=xto(j,1)*bd(2,k)+xto(j,2)*bd(1,k)
c               else 
c                  dge(i,j,k)=xto(j,i)*bd(i,k)
c               endif
cccccccccccccccccccc
c               do m=1,2
c                  do n=1,nis
c                     if (m .eq. j) then
c                        if (i .eq. 3) then
c                           ddge(i,j,m,k,n)=bd(1,n)*bd(2,k)+
c     $		bd(2,n)*bd(1,k)
c                        else
c                           ddge(i,j,m,k,n)=bd(i,n)*bd(i,k)
c                        endif
c                     else
c                        ddge(i,j,m,k,n)=0.0d0
c                     endif
c			    enddo
c			 enddo
c		   enddo
c		enddo
c	enddo
      return
      end


