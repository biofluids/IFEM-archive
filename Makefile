.SUFFIXES: .f90
LIBS    = /opt/intel/composerxe-2011.3.167/mkl/lib/libmkl_lapack95.a
OBJ= global_constants.o global_simulation_parameter.o run_variables.o \
r_common.o fluid_variables.o solid_variables.o mpi_variables.o\
interface_variables.o denmesh_variables.o\
delta_nonuniform.o \
ensight_output.o form.o \
meshgen_solid.o meshgen_fluid.o meshgen_interface.o\
read.o parseinput.o correct.o  echoinput.o equal.o error.o \
facemap.o gaussj.o gjinv.o hg.o hypo.o initialize.o length.o main.o \
nondimension.o norm.o quad3d4n.o quad3d8n.o quad2d3n.o quad2d4n.o \
r_bdpd_curr.o r_bdpd_init.o r_element.o r_sbpress.o r_jacob.o r_load.o \
r_nodalf.o r_sboc.o r_scal.o r_scauchy.o r_smaterj.o  r_spiola.o \
r_spiola_viscous.o r_spiola_elastic.o r_spress.o r_sreadinit.o r_sstif.o \
r_sstrain.o r_stang.o r_stoxc.o r_timefun.o rkpmshape2d.o rkpmshape3d.o \
set.o shape.o solid_solver.o solid_update.o update.o vol.o \
data_exchange_FEM.o getinf_el_3d.o determinant.o inverse.o search_3d.o \
migs.o search_inf.o shx_tets.o energy_solid.o energy_fluid.o volcorr.o \
cg.o mergefinf.o readpart.o setnqloc.o search_inf_pa.o getinf_el_3d_pa.o \
edgeele.o nature_pre.o \
givens.o \
communicate_res.o getnorm_pa.o equal_pa.o vector_dot_pa.o \
blockdiagstable.o gmresnew.o blockgmresnew.o \
setnei_new.o communicate_res_ad.o setid_pa.o \
scale_shift_inter.o get_submesh_info.o search_inf_den.o search_inf_inter.o\
get_inter_ele.o form_inter_ele.o form_inter_bc.o\
block_Laplace.o gmres_Laplace.o blockgmres_Laplace.o\
set_element_index.o get_interpoint_Ia.o correct_Ip.o\
get_normal_Bspline.o get_curv_Bspline.o B_Spline.o B_Spline_0order.o B_Spline_1order.o B_Spline_2order.o\
points_regen.o point_projection.o get_arc_Bspline.o get_surten_Bspline.o get_intervel_Bspline.o\

IFEM: $(OBJ)
	mpif90 -g -O0 -o IFEM $(OBJ) $(LIBS)
.f90.o:
	mpif90 -c -g $<
clean:
	rm -rf *.o *.mod IFEM
