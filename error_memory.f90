module error_memory
  implicit none
  save

  integer :: error_id

contains

subroutine alloc_error(var_name,file_name,error_id)
  
  character(len=50) :: var_name,file_name
  integer :: error_id


  if (error_id .ne. 0) then
     write(*,*) "Error while trying to allocate memory!!!"
     write(*,*) " File      :",file_name
	 write(*,*) " Variable  :",var_name
	 stop
  end if

  return
end subroutine alloc_error




subroutine dealloc_error(var_name,file_name,error_id)
  
  character(len=50) :: var_name,file_name
  integer :: error_id

  if (error_id .ne. 0) then
     write(*,*) "Error while trying to deallocate memory!!!"
     write(*,*) " File      :",file_name
	 write(*,*) " Variable  :",var_name
	 stop
  end if

  return
end subroutine dealloc_error

end module error_memory