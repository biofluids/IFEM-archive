!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Module: meshgen_gmsh.f90
!
!  Axel Gerstenberger, NWU, October 2003
!
!  provides subroutines to read mesh information in the native Gmsh format
!  for more information on the output format see the Gmsh documentation
!
!  - read node, element and surface element number
!  - read mesh from Gmsh input format
!  - provides volume and surface mesh information

module meshgen_gmsh
  implicit none
  save

  integer,parameter :: gmsh_input = 1

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_gmsh_node_elem_numbers(nn,ne,ne_surf,filename_length,filename)
  implicit none

  integer,intent(out) :: nn,ne,ne_surf
  integer,intent(in)  :: filename_length
  character(len=filename_length) :: filename

  integer :: file,i,idummy,ntest,ne_temp
  character(len=7) :: char7
  real(8) :: rdummy

  file=21
  open(file, FILE=filename, STATUS="old",readonly)

  read(file,*) char7
  if (char7 == "$NOD   ") then
     read(file,*) nn
     do i = 1,nn
         read(file,*) rdummy,rdummy,rdummy,rdummy
     enddo
  else
     write(*,*) "missing node entries"
     stop
  endif
  
  read(file,*) char7
  if (char7 /= "$ENDNOD") then
     write(*,*) "wrong node number"
     stop
  endif

  read(file,*) char7
  if (char7 == "$ELM   ") then
     ne = 0
     read(file,*) ne_temp
     do i = 1,ne_temp
        read(file,*) idummy,ntest
        if (ntest == 4) then       !...tetraeder
           ne = ne + 1
        endif
        if (ntest == 2) then       !...triangle
           ne_surf = ne_surf + 1
        endif
     enddo
     !write(*,*) ne
     !write(*,*) ne_surf
  else
     write(*,*) "missing element entries"
     stop
  endif

  close(file)

end subroutine read_gmsh_node_elem_numbers




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_gmsh(xyz,ien,ien_surf,bid,nn,ne,ne_surf,nen,nen_surf,filename_length,filename)
  use fluid_variables, only: nsd
  implicit none

  integer,intent(in)  :: nn,ne,ne_surf,nen,nen_surf
  real(8),intent(out) :: xyz(nsd,nn)
  integer,intent(out) :: ien(nen,ne)
  integer,intent(out) :: ien_surf(nen_surf,ne_surf)
  integer,intent(out) :: bid(ne_surf)
  integer,intent(in)  :: filename_length
  character(len=filename_length) :: filename

  integer :: file,i,in,ne_temp,elem_count,elem_surf_count,idummy,ntest,nen_test,bid_test,iien,itest,ien_test(1:9)
  integer :: connect(nn)  !...connects place in vector with node number -> sort nodes from 1 to nn
  character(len=7) :: char7

  file=21
  open(file, FILE=filename, STATUS="old",readonly)
  
  if (nen == 4) then   !...tetraeder

     read(file,*) char7  !...$NOD
     read(file,*) char7  
     do i = 1,nn
        read(file,*) connect(i),xyz(1:nsd,i)
        !write(*,*) "node: ",i," sucessful read"
     enddo

     read(file,*) char7  !...ENDNOD
     write(*,*) " number of nodes read           : ",nn
     read(file,*) char7  !...ELM
     read(file,*) ne_temp

     elem_count = 0
     elem_surf_count = 0
     bid(1:ne_surf) = 0
     do i = 1,ne_temp     
        read(file,*) itest,ntest,bid_test,idummy,nen_test,ien_test(1:nen_test)
        !write(*,*) "element: ",itest

       !...order node numbers from 1 to nn
        do iien = 1,nen_test
           loop1: do in = 1,nn
              if ( ien_test(iien) == connect(in) ) then
                 ien_test(iien) = in
                 exit loop1
              endif
           enddo loop1
        enddo

        !...read volume and surface elements
        if (ntest == 4) then     !...tetraeder

           if (nen_test /= nen) then
              write(*,*) "nodes per elements mismatch in <read_gmsh>"
              stop
           endif

           elem_count = elem_count + 1
           ien(1:nen_test,elem_count) = ien_test(1:nen_test)
           
        elseif(ntest == 2) then  !...triangle

           if (nen_test /= nen_surf) then
              write(*,*) "nodes per surface elements mismatch in <read_gmsh>"
              stop
           endif

        elem_surf_count = elem_surf_count + 1
        ien_surf(1:nen_surf,elem_surf_count) = ien_test(1:nen_test)

        bid(elem_surf_count) = bid_test

        endif
     enddo

     if (ne /= elem_count) then
         write(*,*) "wrong element numbers <read_gmsh>"
         stop
     endif
     if (ne_surf /= elem_surf_count) then
         write(*,*) "wrong surface element numbers <read_gmsh>"
         stop
     endif

     write(*,*) " number of elements read        : ",ne
     write(*,*) " number of surface elements read: ",ne_surf


  elseif (nen == 8) then
     write(*,*) "not implemented"
     stop
  endif

  close(file)

end subroutine read_gmsh

end module meshgen_gmsh
