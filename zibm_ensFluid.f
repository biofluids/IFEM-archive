c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c modified form io11.f file to generate 
c ensight fluid field file
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zibm_ensFluid(d,klok)

      implicit real*8 (a-h,o-z)
	include 'global.h'

	real* 8 d(ndf,nn)

      character*5  fileroot
      character*12 name_file1
      character*12 name_file2

      parameter (ifileunit = 15)

c%%%%%%%%%%%%%%%%%%%%%%%
      integer part
      data part/2/

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

      write(name_file1,'(A7, A5)')  'ibm.vel', fileroot
      write(name_file2,'(A7, A5)')  'ibm.pre', fileroot

      write(*,*) 'file name is ', name_file1
      write(*,*) 'file name is ', name_file2
 
c     Write velocity output in ens_movie.vel*
109	format('part',i2)
110	format(6e12.5)
      open(ifileunit, file=name_file1, form='formatted')
      write(ifileunit, '(A)') 'fluid field: velocity'
      write(ifileunit, 109) part
	write(ifileunit, *) 'block'	
	write(ifileunit,110) (d(1,i),d(2,i),d(3,i),i=1,nn)

      close(ifileunit)

c     Write pressure output in ens_movie.pre*
      open(ifileunit, file=name_file2, form='formatted')
      write(ifileunit, '(A)') 'fluid field: pressure'
      write(ifileunit, 109) part
	write(ifileunit, *) 'block'
	write(ifileunit,110) (d(4,i),i=1,nn)
      close(ifileunit)

      return
      end




