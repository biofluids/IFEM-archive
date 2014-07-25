if (myid==0) write(*,*) "nn_PML=", nn_PML
if (nn_PML > 0) then
    deallocate(p)
    deallocate(w)
    deallocate(dg)
    allocate(p(ndf,nn+nn_PML),stat=error_id)
    allocate(w(ndf,nn+nn_PML),stat=error_id)
    allocate(dg(ndf,nn+nn_PML),stat=error_id)
endif