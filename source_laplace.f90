subroutine source_laplace(x_fluid,nn_fluid,nsd,x_solid,nn_solid,ien_solid,&
		    ne_solid,nen_solid,ien_sbc,ne_sbc,node_sbc,nn_sbc,data_fluid)
! get source term on fluid node for laplace equation
    use delta_nonuniform, only: cnn, ncnn, shrknode
    use fluid_variables, only: maxconn
    use mpi_variables
    implicit none
    !--------------------------
    real(8) x_fluid(nsd,nn_fluid)
    integer nn_fluid
    integer nsd
    real(8) x_solid(nsd,nn_solid)
    integer nn_solid
    integer ien_solid(ne_solid,nen_solid)
    integer ne_solid
    integer nen_solid
    integer ien_sbc(ne_sbc,nen_solid+2)
    integer ne_sbc
    integer node_sbc(nn_sbc)
    integer nn_sbc
    ! Geometry information
    !------------------------
    !real(8) snvolume(nn_solid) ! solid nodal volume 
    !real(8) adist(nsd,nn_solid) ! region radius
    ! for RKPM interpolation from solid to fluid
    !------------------------
    real(8) x(nsd)
    real(8) rksh(maxconn)
    integer ninf
    integer inf(maxconn)
    real(8) data_fluid(nsd,nn_fluid)
    real(8) data_solid(nsd,nn_solid)

!------------------------
    integer inn
    integer imm
    integer icnn
    integer pt
    integer snode
    integer fnode

    data_fluid(:,:)=0.0
    !data_solid(:,:)=0.0
    if (nsd == 2) then
        call solid_normint(x_solid,nsd,nn_solid,ien_sbc,ne_sbc,nen_solid,&
    			             ien_solid,ne_solid,data_solid)
    else
        if(myid == 0) write(*,*) 'call 3D outward norm calculation'
        call solid_normint_3d(x_solid,nsd,nn_solid,ien_sbc,ne_sbc,nen_solid,&
                            ien_solid,ne_solid,data_solid)
    endif
    
    !call solid_node_volume(x_solid,nsd,nn_solid,ien_solid,ne_solid,nen_solid,snvolume,adist)
    do imm=1,nn_sbc
    	inn=node_sbc(imm)
        do icnn=1,ncnn(inn)
            pt=cnn(icnn,inn)
            data_fluid(1:nsd,pt) = data_fluid(1:nsd,pt) + data_solid(1:nsd,inn) * shrknode(icnn,inn)
        enddo
    enddo
return
end subroutine source_laplace