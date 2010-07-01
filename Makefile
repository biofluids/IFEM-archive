.SUFFIXES: .f90

LIBS =/usr/lib64/liblapack.so \
      /usr/lib64/libblas.so
OBJ= global_constants.o global_simulation_parameter.o run_variables.o \
r_common.o fluid_variables.o solid_variables.o delta_nonuniform.o \
ensight_output.o form.o blockdiagstable.o blockgmresnew.o gmresnew.o meshgen_solid.o \
meshgen_fluid.o read.o parseinput.o correct.o  echoinput.o equal.o error.o \
facemap.o gaussj.o gjinv.o hg.o hypo.o initialize.o length.o main.o \
nondimension.o norm.o quad3d4n.o quad3d8n.o quad2d3n.o quad2d4n.o \
r_bdpd_curr.o r_bdpd_init.o r_element.o r_sbpress.o r_jacob.o r_load.o \
r_nodalf.o r_sboc.o r_scal.o r_scauchy.o r_smaterj.o  r_spiola.o \
r_spiola_viscous.o r_spiola_elastic.o r_spress.o r_sreadinit.o r_sstif.o \
r_sstrain.o r_stang.o r_stoxc.o r_timefun.o rkpmshape2d.o rkpmshape3d.o \
set.o shape.o solid_solver.o solid_update.o update.o vol.o \
data_exchange_FEM.o getinf_el_3d.o determinant.o inverse.o search_3d.o \
migs.o search_inf.o shx_tets.o energy_solid.o energy_fluid.o volcorr.o \
cg.o \
givens.o


IFEM: $(OBJ)
	ifort -g -O0 -o IFEM $(OBJ) $(LIBS)
.f90.o:
	ifort -c -g $<
clean:
	rm -rf *.o *.mod IFEM
