c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c modified from r_main2.f file to 
c generate ensight stress and strain file
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zfem_ensStr(klok, 
     $     dlptlocal_number,
     $     dlptlocal_head, 
     $     dlptlocal_tail)
      implicit real*8 (a-h,o-z)
      include 'r_common'
      include 'main_common'

      integer dlptlocal_number
      integer dlptlocal_head
      integer dlptlocal_tail

c      real*8       cstr(4,nel,3,3)
c      real*8       ge(4,nel,3,3) 
      
      dimension strain(8,max_point)

      character*5  fileroot
      character*17 stress_file1
      character*17 stress_file2
      character*17 stress_file3
      character*17 stress_file4
	character*17 stress_file5
      character*17 stress_file6

      character*17 strain_file1
      character*17 strain_file2
      character*17 strain_file3
      character*17 strain_file4
      character*17 strain_file5
      character*17 strain_file6

      parameter (ifileunit = 15)

      integer part
      data part/1/

c%%%%%%%%%%%%%%%%%%%%%%
c modified from r_main2.f, interpolation using gauss poins 
c%%%%%%%%%%%%%%%%%%%%%%%%
      do i=1,nnda
         do k=1,6
            tstress(k,i)=0.0d0
            tstrain(k,i)=0.0d0
	   enddo
         ntem=0
         do j=1,numela
            do m=1,nis
               if (nea(j,m) .eq. i) then
                  ntem=ntem+1
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
	enddo

c%%%%%%%%%%%%%%%%%%%%%%%
c ensight output
c%%%%%%%%%%%%%%%%%%%%%%


 200  format(a5   )
 201  format(a4,i1)
 202  format(a3,i2)
 203  format(a2,i3)
 204  format(a1,i4)
 205  format(   i5)

      if (klok .eq. 0) then
         write(fileroot, 200) '00000'
      elseif (klok .lt. 10) then
         write(fileroot, 201) '0000',klok
      elseif (klok .lt. 100) then
         write(fileroot, 202) '000' ,klok
      elseif (klok .lt. 1000) then
         write(fileroot, 203) '00'  ,klok
      elseif (klok .lt. 10000) then
         write(fileroot, 204) '0'   ,klok
      elseif (klok .lt. 100000) then    
         write(fileroot, 205)   ''  ,klok
      else
         write(0,*) 'klok .ge. 100000: modify subroutine createfileroot'
         call exit(1)
      endif

      write(stress_file1,'(A12, A5)')  'fem.stressXX', fileroot
      write(stress_file2,'(A12, A5)')  'fem.stressYY', fileroot
      write(stress_file3,'(A12, A5)')  'fem.stressZZ', fileroot
      write(stress_file4,'(A12, A5)')  'fem.stressYZ', fileroot
      write(stress_file5,'(A12, A5)')  'fem.stressXZ', fileroot
      write(stress_file6,'(A12, A5)')  'fem.stressXY', fileroot

      write(strain_file1,'(A12, A5)')  'fem.strainXX', fileroot
      write(strain_file2,'(A12, A5)')  'fem.strainYY', fileroot
      write(strain_file3,'(A12, A5)')  'fem.strainZZ', fileroot
      write(strain_file4,'(A12, A5)')  'fem.strainYZ', fileroot
      write(strain_file5,'(A12, A5)')  'fem.strainXZ', fileroot
      write(strain_file6,'(A12, A5)')  'fem.strainXY', fileroot

      write(*,*) 'file name is ', stress_file1
      write(*,*) 'file name is ', stress_file2
      write(*,*) 'file name is ', stress_file3
      write(*,*) 'file name is ', stress_file4
      write(*,*) 'file name is ', stress_file5
      write(*,*) 'file name is ', stress_file6

      write(*,*) 'file name is ', strain_file1
      write(*,*) 'file name is ', strain_file2
      write(*,*) 'file name is ', strain_file3
      write(*,*) 'file name is ', strain_file4
      write(*,*) 'file name is ', strain_file5
      write(*,*) 'file name is ', strain_file6

109	format('part',i2)
101	format (6e12.5)
c stress xx
      open(ifileunit, file=stress_file1, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress xx'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstress(1,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c stress yy
      open(ifileunit, file=stress_file2, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress yy'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstress(2,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c stress zz
      open(ifileunit, file=stress_file3, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress zz'
 	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstress(3,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c stress yz
      open(ifileunit, file=stress_file4, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress yz'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstress(4,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c stress xz
      open(ifileunit, file=stress_file5, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress yz'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstress(5,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c stress xy
      open(ifileunit, file=stress_file4, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress xy'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstress(6,i)*unit_pressure,i=1,nnda)
      close(ifileunit)


c strain xx
      open(ifileunit, file=strain_file1, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain xx'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstrain(1,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c strain yy
      open(ifileunit, file=strain_file2, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain yy'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstrain(2,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c strain zz
      open(ifileunit, file=strain_file3, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain zz'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstrain(3,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c strain yz
      open(ifileunit, file=strain_file4, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain yz'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstrain(4,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c strain xz
      open(ifileunit, file=strain_file5, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain xz'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstrain(5,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

c strain xy
      open(ifileunit, file=strain_file6, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain xy'
	write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,101) (tstrain(6,i)*unit_pressure,i=1,nnda)
      close(ifileunit)

      return
      end




