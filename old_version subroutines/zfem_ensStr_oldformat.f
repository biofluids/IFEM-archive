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

      character*17 strain_file1
      character*17 strain_file2
      character*17 strain_file3
      character*17 strain_file4

      parameter (ifileunit = 15)

      integer part
      data part/1/

c%%%%%%%%%%%%%%%%%%%%%%
c modified from r_main2.f, interpolation using gauss poins 
c%%%%%%%%%%%%%%%%%%%%%%%%
      do 540 i=1,nnda
         do 542 k=1,4
            tstress(k,i)=0.0d0
            tstrain(k,i)=0.0d0
 542     continue
         ntem=0
         do 541 j=1,numela
            do 543 m=1,nis
               if (nea(j,m) .eq. i) then
                  ntem=ntem+1
                  do 545 n=1,4
                     tstress(n,i)=tstress(n,i)+cstr(n,j,1,1)
                     tstrain(n,i)=tstrain(n,i)+ge(n,j,1,1)
 545              continue
                  go to 541
               endif
 543        continue
 541     continue
         do 546 k=1,4
            tstress(k,i)=tstress(k,i)/ntem
            tstrain(k,i)=tstrain(k,i)/ntem
 546     continue
 540  continue

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
      write(stress_file4,'(A12, A5)')  'fem.stressXZ', fileroot

      write(strain_file1,'(A12, A5)')  'fem.strainXX', fileroot
      write(strain_file2,'(A12, A5)')  'fem.strainYY', fileroot
      write(strain_file3,'(A12, A5)')  'fem.strainZZ', fileroot
      write(strain_file4,'(A12, A5)')  'fem.strainXZ', fileroot

      write(*,*) 'file name is ', stress_file1
      write(*,*) 'file name is ', stress_file2
      write(*,*) 'file name is ', stress_file3
      write(*,*) 'file name is ', stress_file4
      write(*,*) 'file name is ', strain_file1
      write(*,*) 'file name is ', strain_file2
      write(*,*) 'file name is ', strain_file3
      write(*,*) 'file name is ', strain_file4

c stress xx
      open(ifileunit, file=stress_file1, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress xx'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'coordinates'
c      do 2000 i = 1, dlptlocal_tail 
      do 2000 i = 1, nnda 
         write(ifileunit,'(E12.5)') tstress(1,i)*unit_pressure
 2000  continue
      close(ifileunit)

c stress yy
      open(ifileunit, file=stress_file2, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress yy'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'coordinates'
c      do 2001 i = 1, dlptlocal_tail 
      do 2001 i = 1, nnda
         write(ifileunit,'(E12.5)') tstress(4,i)*unit_pressure
 2001 continue
      close(ifileunit)

c stress zz
      open(ifileunit, file=stress_file3, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress zz'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'coordinates'
c      do 2003 i = 1, dlptlocal_tail 
      do 2003 i = 1, nnda
         write(ifileunit,'(E12.5)') tstress(2,i)*unit_pressure
 2003 continue
      close(ifileunit)

c stress xz
      open(ifileunit, file=stress_file4, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: stress xz'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'coordinates'
c      do 2004 i = 1, dlptlocal_tail 
      do 2004 i = 1, nnda 
         write(ifileunit,'(E12.5)') tstress(3,i)*unit_pressure
 2004 continue
      close(ifileunit)


c strain xx
      open(ifileunit, file=strain_file1, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain xx'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'coordinates'
c      do 2100 i = 1, dlptlocal_tail 
      do 2100 i = 1, nnda
         write(ifileunit,'(E12.5)') tstrain(1,i)*unit_pressure
 2100  continue
      close(ifileunit)

c strain yy
      open(ifileunit, file=strain_file2, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain yy'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'coordinates'
c      do 2101 i = 1, dlptlocal_tail 
      do 2101 i = 1, nnda
         write(ifileunit,'(E12.5)') tstrain(4,i)*unit_pressure
 2101 continue
      close(ifileunit)

c strain zz
      open(ifileunit, file=strain_file3, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain zz'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'coordinates'
c      do 2103 i = 1, dlptlocal_tail 
      do 2103 i = 1, nnda
         write(ifileunit,'(E12.5)') tstrain(2,i)*unit_pressure
 2103 continue
      close(ifileunit)

c strain xz
      open(ifileunit, file=strain_file4, form='formatted')
      write(ifileunit, '(A)') 'after interpolation: strain xz'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'coordinates'
c      do 2104 i = 1, dlptlocal_tail 
      do 2104 i = 1, nnda
         write(ifileunit,'(E12.5)') tstrain(3,i)*unit_pressure
 2104 continue
      close(ifileunit)

      return
      end




