subroutine setidpml_pa(res,ndf,nn,nn_PML,id,node_local,nn_local)
! set essenstial BC parallelly 
    integer ndf
    integer nn, nn_PML
    real(8) res(ndf,nn+nn_PML)
    integer id(ndf,nn)
    integer nn_local
    integer node_local(nn_local)

    integer node
    integer i
    integer j

    do i=1,nn_local
        node=node_local(i)
        do j=1,ndf
            if (id(j,node) .eq. 0) res(j,node)=0.0
        enddo
    enddo

return
end subroutine setidpml_pa