subroutine gmres(x,d,dold,w,bg,dg,hg,ien,fext,id,&
                 ne_local,ien_local,node_local,nn_local,&
                 global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length,&
                 fden,fvis,I_fluid,rngface)
    use fluid_variables, only: nsd,nn,ne,nen,ndf,inner,outer,neface
    use solid_variables, only: nn_solid
    use mpi_variables
    use pml_variables
    implicit none
    include 'mpif.h'
	real* 8 x(nsd,nn)
    integer id(ndf,nn)
	real* 8 d(ndf,nn), dold(ndf,nn),hg(ne),fext(ndf,nn)
    integer ien(nen,ne)
	real* 8 bg(ndf*(nn+nn_PML)), dg(ndf*(nn+nn_PML)), w(ndf*(nn+nn_PML))
	real* 8 Hm(inner+1,inner) !Henssenberg matrix
!=========================================================
! reduce dimension to save memory
	real* 8 Vm(ndf*(nn_local+nn_PML_local), inner+1) ! Krylov space matrix
!        real* 8 Vm(ndf*nn, inner+1) ! Krylov space matrix
!=========================================================
    integer i,j,iouter,icount,INFO
	real* 8 e1(inner+1)
	real* 8 x0(ndf*(nn+nn_PML))
	real* 8 beta(inner+1)
	real* 8 eps
	real* 8 r0(ndf*(nn+nn_PML))
	real* 8 rnorm, rnorm0, err
    real* 8 dv(ndf*(nn+nn_PML))
	real* 8 Vy(ndf*(nn+nn_PML))
    real* 8 vloc(ndf,nn+nn_PML), avloc(ndf*(nn+nn_PML))
!==============================
! reduce dimension save memory
	real* 8 temp
!==============================
    character(1) TRAN
	real* 8 workls(2*inner)
	real* 8 av_tmp(ndf,nn+nn_PML)
    integer rngface(neface,ne)
!---------------------------------------
    real(8) fden(nn)
    real(8) fvis(nn)
    real(8) I_fluid(nn)
!============================
! MPI varibalbes
    integer ne_local ! # of element on each processor
    integer ien_local(ne_local) ! subregion-->wholeregion element index
    integer ie_local ! loop parameter
    integer node_local(nn_local)
    integer nn_local
    integer jcount
    integer kcount
    integer node
    integer nn_global_com
    integer global_com(nn_global_com)  ! global node index for communication
    integer nn_local_com
    integer local_com(nn_local_com)  ! local index in the communication region on each processor
    integer ad_length
    integer send_address(ad_length,2)
    integer flag
    real(8) dg_sent(ndf*(nn+nn_PML))
    real(8) space1(inner)
    real(8) space2(inner)
    integer time_arrary_0(8)
    integer time_arrary_1(8)
    real(8) start_time
    real(8) end_time
    real(8) linerr
!---------------------------------------------
    eps = 1.0e-6
    linerr = 1.0e-6
    e1(:) = 0.0
    e1(1) = 1.0
    x0(:) = 0.0
    iouter = 1
    r0(:) = bg(:)
    TRAN = 'N'
    av_tmp(:,:) = 0
    avloc(:) = 0
!	w(:) = 1
    call getnormpml_pa(r0,ndf,nn,node_local,nn_local,rnorm0)
    rnorm = sqrt(rnorm0)

!!!!!!!!!!!!!!!start outer loop!!!!!!!!!
do 111, while((iouter .le. outer) .and. (rnorm .ge. linerr))

!==============================
    Vm(:,:) = 0.0d0
! Already is a local variable defined on each proc
!===============================
    do icount = 1, nn_local
        node=node_local(icount)
        do jcount=1,ndf
            Vm((icount-1)*ndf+jcount,1) = r0((node-1)*ndf+jcount)/rnorm
        enddo
        if (seqcPML(node) > 0) then
            do jcount=1,ndf
                Vm((nn_local+seqcPMLlocal(icount)-1)*ndf+jcount,1) = r0((nn+seqcPML(node)-1)*ndf+jcount)/rnorm
            enddo
        endif
    enddo ! get V1

    beta(:) = rnorm*e1(:) ! get beta*e1
    Hm(:,:) = 0.0
!!!!!!!!!!!!!!!!start inner loop!!!!!!!!!!!!!
    do j=1,inner
        do icount=1, nn_local
            node=node_local(icount)
            do jcount=1,ndf
                dv((node-1)*ndf+jcount) = eps/w((node-1)*ndf+jcount)*Vm((icount-1)*ndf+jcount,j)
            enddo
            if (seqcPML(node) > 0) then
                do jcount=1,ndf
                    dv((nn+seqcPML(node)-1)*ndf+jcount) = eps*Vm((nn_local+seqcPMLlocal(icount)-1)*ndf+jcount,j) &
                                                          /w((nn+seqcPML(node)-1)*ndf+jcount)
                enddo
            endif
        enddo!!!!!!!!!!calcule eps*inv(P)*V1
!===============================================
!		vloc(:,:) = 0.0d0
! Clear the matrix locally to avoid too many loops 
        do icount=1,nn_local
            node=node_local(icount)
            vloc(1:ndf,node)=0.0d0
            if (seqcPML(node) > 0) then
                vloc(1:ndf,nn+seqcPML(node))=0.0d0
            endif
        enddo
! clear all the processor internal nodes
        do icount=1,nn_local_com
            node=global_com(local_com(icount))
            vloc(1:ndf,node)=0.0d0
            if (seqcPML(node) > 0) then
                vloc(1:ndf,nn+seqcPML(node))=0.0d0
            endif
        enddo
! clear the processor boundary nodes
!===============================================
        call equalpml_pa(dv,vloc,ndf,nn,node_local,nn_local)
!============================
        do icount=1,nn_local
            node=node_local(icount)
            vloc(1:ndf,node)=vloc(1:ndf,node)+d(1:ndf,node)
            if (seqcPML(node) > 0) then
                vloc(1:ndf,nn+seqcPML(node))=vloc(1:ndf,nn+seqcPML(node))+qv(1:ndf,node)
            endif
        enddo
! Let vloc=vloc+d first then communicate, and then it should same # of loop (avoiding loop at the whole domain)
!=============================
        call communicate_respml_ad_sub(vloc,ndf,nn,send_address,ad_length)
!----------------------------------------------------------------------------------------------
!vloc(:,:)=vloc(:,:)+d(:,:)
!----------------------------------------------------------------------------------------------
! Clear matrix local first interal then boundary processor nodes
        do icount=1,nn_local
            node=node_local(icount)
            av_tmp(1:ndf,node)=0.0d0
            if (seqcPML(node) > 0) then
                av_tmp(1:ndf,nn+seqcPML(node))=0.0d0
            endif
        enddo
        do icount=1,nn_local_com
            node=global_com(local_com(icount))
            av_tmp(1:ndf,node)=0.0d0
            if (seqcPML(node) > 0) then
                av_tmp(1:ndf,nn+seqcPML(node))=0.0d0
            endif
        enddo
!		av_tmp(:,:) = 0.0d0
!======================================
!!!!!! av_tmp is a temporary residual
        call blockgmresnew(x,vloc,dold,av_tmp,hg,ien,fext,ne_local,ien_local,node_local,nn_local,&
                           fden,fvis,I_fluid,rngface)
        call communicate_respml_ad_sub(av_tmp,ndf,nn,send_address,ad_length)
!===================
! avloc(:)=0.0d0
!==================
        call equalpml_pa(av_tmp,avloc,ndf,nn,node_local,nn_local)
        do icount=1, nn_local
            node=node_local(icount)
            do jcount=1,ndf
                kcount=(node-1)*ndf+jcount
                avloc(kcount) = (-avloc(kcount)+bg(kcount))/eps ! get Av,bg=-r(u)
            enddo
            if (seqcPML(node) > 0) then
                do jcount=1,ndf
                    kcount=(nn+seqcPML(node)-1)*ndf+jcount
                    avloc(kcount) = (-avloc(kcount)+bg(kcount))/eps ! get Av,bg=-r(u)
                enddo
            endif
        enddo
        call setidpml_pa(avloc,ndf,nn,nn_PML,id,node_local,nn_local)
        space1(:)=0.0d0
        space2(:)=0.0d0
        end_time=mpi_wtime()
        do i=1,j
            do icount=1,nn_local
                node=node_local(icount)
                do jcount=1,ndf
                    kcount=(node-1)*ndf+jcount
                    space1(i)=space1(i)+avloc(kcount)*Vm((icount-1)*ndf+jcount,i)
                enddo
                if (seqcPML(node) > 0) then
                    do jcount=1,ndf
                        kcount=(nn+seqcPML(node)-1)*ndf+jcount
                        space1(i)=space1(i)+avloc(kcount)*Vm((nn_local+seqcPMLlocal(icount)-1)*ndf+jcount,i)
                    enddo
                endif
            enddo
! Do the vector product of v_i * v_i set it to be h_i,j
        enddo  ! construct AVj and hi,j
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_allreduce(space1(1),space2(1),j,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        end_time=mpi_wtime()-end_time
        Hm(1:j,j)=space2(1:j)
        do icount = 1, nn_local
            node=node_local(icount)
            do jcount=1,ndf
                kcount=(icount-1)*ndf+jcount
                do i=1,j
                    Vm(kcount,j+1) = Vm(kcount,j+1)-Hm(i,j)*Vm(kcount,i)
                enddo
                Vm(kcount,j+1)=Vm(kcount,j+1)+avloc((node-1)*ndf+jcount)
            enddo
            if (seqcPML(node) > 0) then
                do jcount=1,ndf
                    kcount=(nn_local+seqcPMLlocal(icount)-1)*ndf+jcount
                    do i=1,j
                        Vm(kcount,j+1) = Vm(kcount,j+1)-Hm(i,j)*Vm(kcount,i)
                    enddo
                    Vm(kcount,j+1)=Vm(kcount,j+1)+avloc((nn+seqcPML(node)-1)*ndf+jcount)
                enddo
            endif
        enddo  ! construct v(j+1)
        temp=0.0d0
        do icount = 1, nn_local
            node=node_local(icount)
            do jcount=1,ndf
                kcount=(icount-1)*ndf+jcount
                temp=temp+Vm(kcount,j+1)*Vm(kcount,j+1)
            enddo
            if (seqcPML(node) > 0) then
                do jcount=1,ndf
                    kcount=(nn_local+seqcPMLlocal(icount)-1)*ndf+jcount
                    temp=temp+Vm(kcount,j+1)*Vm(kcount,j+1)
                enddo
            endif
        enddo
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_allreduce(temp,rnorm0,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

        Hm(j+1,j) = sqrt(rnorm0)
        do icount = 1, nn_local
            node=node_local(icount)
            do jcount=1,ndf
                kcount=(icount-1)*ndf+jcount
                Vm(kcount,j+1)=Vm(kcount,j+1)/Hm(j+1,j)
            enddo
            if (seqcPML(node) > 0) then
                do jcount=1,ndf
                    kcount=(nn_local+seqcPMLlocal(icount)-1)*ndf+jcount
                    Vm(kcount,j+1)=Vm(kcount,j+1)/Hm(j+1,j)
                enddo
            endif
        enddo
    enddo  ! end inner loop
!====================================================================================		
! Use lapack to solve the LS problem
!	call DGELS(TRAN,inner+1,inner,1,Hm,inner+1,beta,inner+1,workls,2*inner,INFO)
!===================================================================================
! Use Givens rotation to solve the LS problem
    call mpi_barrier(mpi_comm_world,ierror)
    call givens(Hm,inner,beta)
!!!!!!!!!!!!beta(1:inner) is ym, the solution!!!!!!!!!!!!!!!!!!!!!!!!!
    Vy(:) = 0
    do icount=1, nn_local
        node=node_local(icount)
        do jcount=1,ndf
            kcount=(node-1)*ndf+jcount
            do i=1,inner
                Vy(kcount)=Vy(kcount)+Vm((icount-1)*ndf+jcount,i)*beta(i)
            enddo
            x0(kcount)=x0(kcount)+Vy(kcount)
        enddo
        if (seqcPML(node) > 0) then
            do jcount=1,ndf
                kcount=(nn+seqcPML(node)-1)*ndf+jcount
                do i=1,inner
                    Vy(kcount)=Vy(kcount)+Vm((nn_local+seqcPMLlocal(icount)-1)*ndf+jcount,i)*beta(i)
                enddo
                x0(kcount)=x0(kcount)+Vy(kcount)
            enddo
        endif
    enddo ! calculate Xm
    do icount = 1, nn_local
        node=node_local(icount)
        do jcount=1,ndf
            kcount=(node-1)*ndf+jcount
            dv(kcount) = eps/w(kcount)*x0(kcount)
        enddo
        if (seqcPML(node) > 0) then
            do jcount=1,ndf
                kcount=(nn+seqcPML(node)-1)*ndf+jcount
                dv(kcount) = eps/w(kcount)*x0(kcount)
            enddo
        endif
    enddo
!============================
! Clear matirx locally 
    do icount=1,nn_local
        node=node_local(icount)
        vloc(1:ndf,node)=0.0d0
        if (seqcPML(node) > 0) then
            vloc(1:ndf,nn+seqcPML(node))=0.0d0
        endif
    enddo
! clear all the processor internal nodes
    do icount=1,nn_local_com
        node=global_com(local_com(icount))
        vloc(1:ndf,node)=0.0d0
        if (seqcPML(node) > 0) then
            vloc(1:ndf,nn+seqcPML(node))=0.0d0
        endif
    enddo
! clear the processor boundary nodes
!	vloc(:,:) = 0
!==============================
    call equalpml_pa(dv,vloc,ndf,nn,node_local,nn_local)
    do icount=1,nn_local
        node=node_local(icount)
        vloc(1:ndf,node)=vloc(1:ndf,node)+d(1:ndf,node)
        if (seqcPML(node) > 0) then
            vloc(1:ndf,nn+seqcPML(node))=vloc(1:ndf,nn+seqcPML(node))+qv(1:ndf,node)
        endif
    enddo
    call communicate_respml_ad_sub(vloc,ndf,nn,send_address,ad_length)
!==============================
!	vloc(:,:) = vloc(:,:)+d(:,:)
!===============================
! Clear matrix local first interal then boundary processor nodes
    do icount=1,nn_local
        node=node_local(icount)
        av_tmp(1:ndf,node)=0.0d0
        if (seqcPML(node) > 0) then
            av_tmp(1:ndf,nn+seqcPML(node))=0.0d0
        endif
    enddo
    do icount=1,nn_local_com
        node=global_com(local_com(icount))
        av_tmp(1:ndf,node)=0.0d0
        if (seqcPML(node) > 0) then
            av_tmp(1:ndf,nn+seqcPML(node))=0.0d0
        endif
    enddo
!	av_tmp(:,:) = 0.0d0
!================================
    call blockgmresnew(x,vloc,dold,av_tmp,hg,ien,fext,ne_local,ien_local,node_local,nn_local,&
                       fden,fvis,I_fluid,rngface)
    call communicate_respml_ad_sub(av_tmp,ndf,nn,send_address,ad_length)
!==================
!avloc(:)=0.0d0
!==================
    call equalpml_pa(av_tmp,avloc,ndf,nn,node_local,nn_local)
    call setidpml_pa(avloc,ndf,nn,nn_PML,id,node_local,nn_local)
    do icount=1, nn_local
        node=node_local(icount)
        do jcount=1,ndf
            kcount=(node-1)*ndf+jcount
            avloc(kcount) = (-avloc(kcount)+bg(kcount))/eps
        enddo
        if (seqcPML(node) > 0) then
            do jcount=1,ndf
                kcount=(nn+seqcPML(node)-1)*ndf+jcount
                avloc(kcount) = (-avloc(kcount)+bg(kcount))/eps
            enddo
        endif
    enddo
!!!!!!!!!!calculate AXm
    do icount=1,nn_local
        node=node_local(icount)
        do jcount=1,ndf
            kcount=(node-1)*ndf+jcount
            r0(kcount) = bg(kcount)-avloc(kcount)
        enddo
        if (seqcPML(node) > 0) then
            do jcount=1,ndf
                kcount=(nn+seqcPML(node)-1)*ndf+jcount
                r0(kcount) = bg(kcount)-avloc(kcount)
            enddo
        endif
    enddo !update r0=f-AX0

    call getnormpml_pa(r0,ndf,nn,node_local,nn_local,rnorm0)
    err = sqrt(rnorm0)
    rnorm = sqrt(rnorm0)
    iouter = iouter + 1        
    if (myid == 0) write(*,*) '||r0||=||f-AX0||=',err

111 continue  ! end outer loop

    dg_sent(1:ndf*(nn+nn_PML))=0.0d0
    do icount=1, nn_local
        node=node_local(icount)
        do jcount=1,ndf
            kcount=(node-1)*ndf+jcount
            dg_sent(kcount) = 1/w(kcount)*x0(kcount) ! unscaled x0
        enddo
        if (seqcPML(node) > 0) then
            do jcount=1,ndf
                kcount=(nn+seqcPML(node)-1)*ndf+jcount
                dg_sent(kcount) = 1/w(kcount)*x0(kcount) ! unscaled x0
            enddo
        endif
    enddo
    dg(1:ndf*(nn+nn_PML))=0.0d0
    call mpi_barrier(mpi_comm_world,ierror)
    call mpi_allreduce(dg_sent(1),dg(1),ndf*(nn+nn_PML),mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
!	call mpi_bcast(dg(1),ndf*nn,mpi_double_precision,0,mpi_comm_world,ierror)
    call mpi_barrier(mpi_comm_world,ierror)

return
end subroutine gmres