subroutine test_inf(infdomain_sub)
use interface_variables

integer infdomain_sub(maxmatrix)
write(*,*)'inf_test=',infdomain_sub(1:nn_inter)

stop
end subroutine test_inf
