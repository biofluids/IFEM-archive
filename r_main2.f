      subroutine r_main2(nn,klok)
      implicit real*8 (a-h,o-z)
      include 'r_common'

      include 'main_common'            
  
c     output to tecplot 
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     smooth element --> nodal average
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=(nslay-1)*nnda+1,nslay*nnda
         do k=1,4
            tstress(k,i)=0.0d0
            tstrain(k,i)=0.0d0
	   enddo
         ntem=0
         pave(i)=0.0d0
         do j=(nslay-1)*numela+1,nslay*numela
            do m=1,nis
               if (nea(j,m) .eq. i) then
                  ntem=ntem+1
c
c     4/1
c     
                  pave(i)=pave(i)+pre(1,j)
                  do n=1,4
                     tstress(n,i)=tstress(n,i)+cstr(n,j,1,1)
                     tstrain(n,i)=tstrain(n,i)+ge(n,j,1,1)
				enddo
                  go to 541
               endif
		  enddo
 541    enddo
         do k=1,4
            tstress(k,i)=tstress(k,i)/ntem
            tstrain(k,i)=tstrain(k,i)/ntem
	   enddo
         pave(i)=pave(i)/ntem
	enddo
c     
      write(9001,850) klok,nnda,numela 
c 850  format('ZONE N=',i5,',',' E=',i5,',',
 850  format('ZONE   T= "',i5,'", N=',i5,', E=',i5,
     $     ', F=FEPOINT, ET=QUADRILATERAL')
c     
      do i=(nslay-1)*nnda+1,nslay*nnda
         write(9001,851) (coor(i,1)+dis(1,i))*unit_length,
     $        (coor(i,2)+dis(2,i))*unit_length,
     $        dis(1,i)*unit_length,dis(2,i)*unit_length,
     $        tstress(1,i)*unit_pressure,
     $        tstress(2,i)*unit_pressure,
     $        tstress(3,i)*unit_pressure,
     $        tstress(4,i)*unit_pressure,
     $        pave(i),
     $        tstrain(1,i),
     $        tstrain(2,i),
     $        tstrain(3,i),
     $        tstrain(4,i)
 851     format(13(e14.6,', '))
	enddo
c     

      do i=1,numela
         write(9001,852) (nea(i,j),j=1,nis)
 852     format(i4,' ',i4,' ',i4,' ',i4)

	enddo
      return
      end