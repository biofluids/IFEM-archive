subroutine test(test_1,n_test)
  integer n_test
  integer,allocatable :: test_1(:)
  integer i
  integer error_id

  allocate(test_1(n_test),stat=error_id)
  do i=1,n_test
     test_1(i)=i
  end do

end subroutine
