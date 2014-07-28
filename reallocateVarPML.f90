if (myid==0) write(*,*) "nn_PML=", nn_PML
if (nn_PML > 0) then
    if (myid==0) write(*,*) "reallocating p, w, dg"
    deallocate(p)
    deallocate(w)
    deallocate(dg)
    allocate(p(ndf,nn+nn_PML),stat=error_id)
    allocate(w(ndf,nn+nn_PML),stat=error_id)
    allocate(dg(ndf,nn+nn_PML),stat=error_id)
else
    if (myid==0) write(*,*) "no reallocation of p, w, dg necessary"
endif