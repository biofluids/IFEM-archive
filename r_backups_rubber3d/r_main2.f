      subroutine r_main2(nn,klok)
      implicit real*8 (a-h,o-z)
      include 'r_common'

      include 'main_common'            
  
c     output to tecplot 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     smooth element --> nodal average
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,nnd
         do k=1,6
            tstress(k,i)=0.0d0
            tstrain(k,i)=0.0d0
	   enddo
         ntem=0
	   pave(i)=0
         do j=1,numel
            do m=1,nis
               if (nea(j,m) .eq. i) then
                  ntem=ntem+1
	            pave(i)=pave(i)+pre(1,j)
                  do n=1,6
                     tstress(n,i)=tstress(n,i)+cstr(n,j,1,1,1)
                     tstrain(n,i)=tstrain(n,i)+ge(n,j,1,1,1)
				enddo
                  go to 541
               endif
		  enddo
541	   enddo
         do k=1,6
            tstress(k,i)=tstress(k,i)/ntem
            tstrain(k,i)=tstrain(k,i)/ntem
	   enddo
	    pave(i)=pave(i)/ntem
	enddo

c     
	if (klok.eq.0) then
	      write(9001,850) klok,nnda,numela 
 850  format('ZONE   T= "',i5,'", N=',i5,', E=',i5,
     $     ', F=FEPOINT, ET=BRICK')
	else
	      write(9001,853) klok,nnda,numela 
 853  format('ZONE   T= "',i5,'", N=',i5,', E=',i5,
     $     ', F=FEPOINT, ET=BRICK, D=(FECONNECT)')
	endif
	     
      do i=1,nnd
         write(9001,851) (coor(i,1)+dis(1,i)),
     $        (coor(i,2)+dis(2,i)),
     $		(coor(i,3)+dis(3,i)),
     $        dis(1,i),
     $		dis(2,i),
     $		dis(3,i),
     $        tstress(1,i),
     $        tstress(2,i),
     $        tstress(3,i),
     $        tstress(4,i),
     $        tstress(5,i),
     $        tstress(6,i),
     $        pave(i),
     $        tstrain(1,i), tstrain(2,i), tstrain(3,i),
     $        tstrain(4,i), tstrain(5,i), tstrain(6,i)
 851     format(19(e12.4,', '))
	enddo
c     
	if (klok.eq.0) then
      do i=1,numela
         write(9001,852) (nea(i,j),j=1,nis)
 852     format(8i8)
	enddo
	endif

      return
      end