c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c modified form io11.f file to generate 
c ensight fluid field file
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zfem_ensFluid(d,klok)

      implicit real*8 (a-h,o-z)
      include 'global.h'

	real* 8 d(ndf,nn)

      character*5  fileroot
      character*12 name_file1
      character*12 name_file2
      character*12 name_file3

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

      write(name_file1,'(A7, A5)')  'fem.vel', fileroot
      write(name_file2,'(A7, A5)')  'fem.pre', fileroot
      write(name_file3,'(A7, A5)')  'fem.vor', fileroot

      write(*,*) 'file name is ', name_file1
      write(*,*) 'file name is ', name_file2
      write(*,*) 'file name is ', name_file3
 
      nskip = 2

c     Write velocity output in ens_movie.vel*
      open(ifileunit, file=name_file1, form='formatted')
      write(ifileunit, '(A)') 'fluid field: velocity'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'block'
      
      do i=1,nn
	  write(ifileunit,'(E12.5)') d(1,i)
	enddo
	
	do i=1,nn
	  write(ifileunit,'(E12.5)') d(2,i)
	enddo
	
	do i=1,nn
	  write(ifileunit,'(E12.5)') d(3,i)
	enddo

      close(ifileunit)


c     Write pressure output in ens_movie.pre*
      open(ifileunit, file=name_file2, form='formatted')
      write(ifileunit, '(A)') 'fluid field: pressure'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'block'
      
	do i=1,nn
	  write(ifileunit,'(E12.5)') d(4,i)
	enddo

      close(ifileunit)

      j=0
c++++++++
c     calculate vorticity
c++++++++

c         do 400 k = mnk, mxk
c            do 403 i= dmnlat1+3, dmxlat1
c               deltau = (u(i,j,k+1) - u(i,j,k-1))*unit_velocity
c               deltaz = 2.0d0*unit_length
c               deltaw = (w(i+1,j,k) - w(i-1,j,k))*unit_velocity
c               deltax = 2.0d0*unit_length
c               vort(i,j,k) = deltau/deltaz - deltaw/deltax
c 403        continue
c 400     continue

c     Write vorticity output in ens_movie.vor*
      open(ifileunit, file=name_file3, form='formatted')

	write(ifileunit, '(A)') 'NOT BEEN CALCULATED NOW, CHECK ZFEM_ENSFL
     +UID.F'
      write(ifileunit, '(A)') 'fluid field: vorticity'
      write(ifileunit, '(A)') 'part'
      write(ifileunit, '(I10)') part
      write(ifileunit, '(A)') 'block'

      
c       do 2200 k = mnk,mxk, nskip
c          do 2201 j =  0,dmxlat2
c             do 2202 i= 3+nbou,dmxlat1-nbou, nskip
c               write(ifileunit,'(E12.5)') 
c     $               vort(i,0,k)
c 2202        continue
c 2201     continue
c 2200  continue
      close(ifileunit)

      return
      end
