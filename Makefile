OBJ=global_constants.o global_simulation_parameter.o\
    run_variables.o r_common.o fluid_variables.o solid_variables.o\
    deformation_gradient.o delta_nonuniform.o \
    ensight_output.o form.o formm.o\
    gmres_fluid_ale.o gmres_fluid_mesh.o\
    meshgen_solid.o meshgen_fluid.o\
    solid_fem_BC.o\
    read.o\
    parseinput.o\
    correct3dl.o\
    echoinput.o\
    equal.o\
    error.o\
    facemap.o\
    fclear.o\
    gaussj.o\
    gjinv.o\
    hg.o\
    hydro.o\
    hypo.o\
    initialize.o\
    lenght.o\
    main.o\
    nondimension.o\
    norm.o\
    quad3d4n.o quad3d8n.o quad2d3n.o quad2d4n.o\
    r_bdpd_curr.o r_bdpd_init.o\
    r_element.o\
    r_sbpress.o\
    r_jacob.o\
    r_load.o\
    r_nodalf.o\
    r_sboc.o\
    r_sbpress.o\
    r_scal.o\
    r_scauchy.o\
    r_smaterj.o\
    r_spiola.o\
    r_spiola_viscous.o\
    r_spress.o\
    r_sreadinit.o\
    r_sstif.o\
    r_sstrain.o\
    r_stang.o\
    r_stoxc.o\
    r_timefun.o\
    rkpmshape3d.o\
    set.o\
    shape.o\
    sharp.o\
    solid_solver.o\
    solid_update.o\
    update.o\
    velocity.o\
    vol.o
    
MKLPATH=/opt/intel/mkl60/lib/32/
FC=ifc
#compile option
#COPT=-i4 -r8 -axK -implicitnone -FR -Vaxlib
COPT=-i4 -r8 -implicitnone -Vaxlib -O2
#link option for dynamic link
#LOPT=-L$(MKLPATH) -lmkl  -lmkl_lapack32 -lmkl_lapack64 -lguide -lm -lpthread
#link option for static link
LOPT=$(MKLPATH)libmkl_lapack.a $(MKLPATH)libmkl_ia32.a $(MKLPATH)libguide.a -L$(MKLPATH) -lm -lpthread
%.o : %.f90
	$(FC) $(COPT) -c $< -o $@
ifem: $(OBJ)
	mpif90 $(LOPT) $(OBJ) -o $@
clean:
	rm -f $(OBJ) main
	rm -f *~
	rm -f *.bak
	rm -f *.pc *.pcl
	rm -f *.d
	rm -f ifem aleifem

