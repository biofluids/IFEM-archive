OBJ= global_constants.o global_simulation_parameter.o run_variables.o \
r_common.o fluid_variables.o solid_variables.o deformation_gradient.o \
delta_nonuniform.o ensight_output.o form.o formm.o gmres_fluid_ale.o \
gmres_fluid_mesh.o meshgen_solid.o meshgen_fluid.o solid_fem_BC.o \
read.o parseinput.o correct3dl.o echoinput.o equal.o error.o \
facemap.o fclear.o gaussj.o gjinv.o hg.o hydro.o hypo.o initialize.o \
length.o main.o nondimension.o norm.o quad3d4n.o quad3d8n.o \
quad2d3n.o quad2d4n.o r_bdpd_curr.o r_bdpd_init.o r_element.o \
r_sbpress.o r_jacob.o r_load.o r_nodalf.o r_sboc.o r_sbpress.o \
r_scal.o r_scauchy.o r_smaterj.o r_spiola.o r_spiola_viscous.o \
r_spress.o r_sreadinit.o r_sstif.o r_sstrain.o r_stang.o r_stoxc.o \
r_timefun.o rkpmshape3d.o set.o shape.o sharp.o solid_solver.o \
solid_update.o update.o velocity.o vol.o 

IFEM: $(OBJ)
	f90 $(OBJ) -o IFEM
.o: .f90
	f90 -c $@

clean:
	rm *.o IFEM
