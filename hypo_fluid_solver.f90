!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   FEM Navier Stokes fluid solver
!     compressible, implicit
! Latest update: Jack, 06/05/2014
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if (myid == 0) then
    write(*,*) "*** Solving Fluids: calculate fluid field ***"
endif
!...store the variables from the previous time step,dold=d 
dold = d
!--------- update PML auxiliary variables -----------
call evaluateQvPML(qv,d)   ! !!!!!!!!!!!!!!!!!!!!!! evaluate for all nodes: remember to CHANGE!!!!!!!!!!
call formid(id,rng,ien)
!----------------------------------------------------
! time varying BC for velocity
call formd_time(d,rng,ien)
!call formd(d,rng,ien)
if (myid == 0) write(*,*) "********** Applied Time Changing BC ***************"

!===========================================================================
! Applying natural(pressure) BC
if (edge_inflow .ne. 0) then
    if (myid == 0) then
        write (*,*) '********** Apply nature BC ***********'
        if (ptotflag == 1) write(*,*) '======== Total Pressure B.C. is used ==========='
    endif
    if (nsd == 2) then
        if (ptotflag == 0) then
            call nature_pre(x,dold,ien,rng,bc4el,ne_inflow,edge_inflow,pin,res_bc)
        else
            call nature_totpre(x,dold,ien,rng,bc4el,ne_inflow,edge_inflow,pin,res_bc,pin_s)
            call nature_pre(x,dold,ien,rng,bc4el,ne_inflow,edge_inflow,pin_s,res_bc)
        endif

        f_fluids(1,:) =  f_fluids(1,:) + res_bc(1,:)
    elseif (nsd == 3) then
        if (ptotflag == 0) then
            call nature_pre_3d(x,dold,ien,rng,bc4el,ne_inflow,edge_inflow,pin,res_bc)
        else
            call nature_totpre_3d(x,dold,ien,rng,bc4el,ne_inflow,edge_inflow,pin,res_bc,pin_s)
            call nature_pre_3d(x,dold,ien,rng,bc4el,ne_inflow,edge_inflow,pin_s,res_bc)
        endif
        f_fluids(1:nsd,:) =  f_fluids(1:nsd,:) + res_bc(1:nsd,:)
    endif
endif
!===========================================================================
do iit=1,nit
    p(:,:) = 0.0d0
    w(:,:) = 0.0d0
    time_com=mpi_wtime()
    call block(x,d,dold,p,w,hg,ien,f_fluids,rng,f_stress,ne_intlocal,ien_intlocal,node_local,nn_local,&
               fden,fvis,I_fluid,qv,sigmaPML)
    time_com=mpi_wtime()-time_com
    if (myid == 0) write(*,*) 'Time for evaluate block', time_com

    time_com=mpi_wtime()
    call communicate_res_ad_sub(p,ndf,nn,send_address,ad_length)
    call communicate_res_ad_sub(w,ndf,nn,send_address,ad_length)
    time_com=mpi_wtime()-time_com
    if (myid == 0) write(*,*) 'Time for communication twice', time_com

    !...apply boundary conditions, calculate residual with corrected BC, i.e. where id=0 -> p=0
!    call setid(p,id,ndf)
    time_com=mpi_wtime()
    call setid_pa(p,ndf,nn,id,node_local,nn_local)
    time_com=mpi_wtime()-time_com
    if (myid == 0) write(*,*) 'Time for set BC', time_com

    !...calculate the normalized error, res_l
    call getnorm_pa(p,ndf,nn,node_local,nn_local,res_l)
    res_l= sqrt(res_l)
    !...calculate the increment(delta_d) for each degree of freedom, dg
    dg(:,:) = 0.0d0
    call gmres(x,d,dold,w,p,dg,hg,ien,f_fluids,id, &
               ne_intlocal,ien_intlocal,node_local,nn_local, &
               global_com,nn_global_com,local_com,nn_local_com,send_address,ad_length,&
               fden,fvis,I_fluid,rng,qv,sigmaPML)
    !...calculate the normalized increment, del_l
    call getnorm_pa(dg,ndf,nn,node_local,nn_local,del_l)
    del_l = sqrt(del_l)
    !...update the calculated variables, d=d+dg
    call mpi_barrier(mpi_comm_world,ierror)

    call update(p, d, dg, ndf)
    if (myid == 0 ) then
        !...output the errors in each iteration
        write(6,102) iit,res_l,del_l
        write(7,102) iit,res_l,del_l
    endif
enddo

if (myid ==0) then
    write(*,'("maximum fluid velocity = ",f13.5)') maxval(d(1:nsd,:))
endif
