c     
c     green-lagrangian strain and derivative
c     
      subroutine r_sstrain(toc,xto,lx,ly,ne)
      implicit real*8 (a-h,o-z)
      include 'r_common'
      dimension xto(2,2),toc(3,3)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     dge(i,j,k)      -- i -- strain, j   -- dir,   k -- node
c     dge(i,j,m,k,n)  -- i -- strain, j,m -- dir, k,n -- node
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 10 i=1,2
         ge(i,ne,lx,ly)=0.5d0*(toc(i,i)-1.0d0)
   10 continue
      ge(3,ne,lx,ly)=toc(1,2)
c     ge(4,ne,lx,ly)=0.0d0
      do 11 i=1,3
         do 12 k=1,nis
            do 13 j=1,2
               if (i .eq. 3) then
                  dge(i,j,k)=xto(j,1)*bd(2,k)+xto(j,2)*bd(1,k)
               else 
                  dge(i,j,k)=xto(j,i)*bd(i,k)
               endif
cccccccccccccccccccc
               do 14 m=1,2
                  do 15 n=1,nis
                     if (m .eq. j) then
                        if (i .eq. 3) then
                           ddge(i,j,m,k,n)=bd(1,n)*bd(2,k)+
     $		bd(2,n)*bd(1,k)
                        else
                           ddge(i,j,m,k,n)=bd(i,n)*bd(i,k)
                        endif
                     else
                        ddge(i,j,m,k,n)=0.0d0
                     endif
   15             continue
   14          continue
   13       continue
   12    continue
   11 continue
      return
      end



