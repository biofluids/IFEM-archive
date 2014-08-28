subroutine equalpml_pa(x,y,ndf,nn,node_local,nn_local)
    use pml_variables
! let y=x for the nodes on its own proc
! y has to be clear to zero first

    real(8) x(ndf,nn+nn_PML)
    real(8) y(ndf,nn+nn_PML)
    integer nn
    integer ndf
    integer node_local(nn_local)
    integer nn_local

    integer i
    integer j
    integer node

    do i=1,nn_local
        node=node_local(i)
        do j=1,ndf
            y(j,node)=x(j,node)
        enddo
        if (seqcPML(node) > 0) then
            do j=1,ndf
                y(j,nn+seqcPML(node))=x(j,nn+seqcPML(node))
            enddo
        endif
    enddo

end subroutine equalpml_pa
