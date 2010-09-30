.SUFFIXES: .f90
LIBS    = /usr/opt/lapack-3.3.1/lapack_LINUX.a \
          /usr/opt/BLAS/blas_LINUX.a
OBJ= global_constants.o global_simulation_parameter.o run_variables.o interface_variables.o centermesh_variables.o\
r_common.o fluid_variables.o solid_variables.o delta_nonuniform.o \
ensight_output.o form.o blockdiagstable.o blockgmresnew.o gmresnew.o meshgen_solid.o \
meshgen_fluid.o meshgen_interface.o read.o parseinput.o correct.o  echoinput.o equal.o error.o \
facemap.o gaussj.o gjinv.o hg.o hypo.o initialize.o length.o main.o \
nondimension.o norm.o quad3d4n.o quad3d8n.o quad2d3n.o quad2d4n.o \
r_bdpd_curr.o r_bdpd_init.o r_element.o r_sbpress.o r_jacob.o r_load.o \
r_nodalf.o r_sboc.o r_scal.o r_scauchy.o r_smaterj.o  r_spiola.o \
r_spiola_viscous.o r_spiola_elastic.o r_spress.o r_sreadinit.o r_sstif.o \
r_sstrain.o r_stang.o r_stoxc.o r_timefun.o rkpmshape2d.o rkpmshape3d.o \
set.o shape.o solid_solver.o solid_update.o update.o vol.o \
data_exchange_FEM.o getinf_el_3d.o determinant.o inverse.o search_3d.o \
migs.o search_inf.o shx_tets.o energy_solid.o energy_fluid.o volcorr.o \
cg.o get_inter_ele.o form_inter_ele.o form_inter_bc.o block_Laplace.o\
scale_shift_inter.o gmres_Laplace.o blockgmres_Laplace.o \
set_element_index.o get_interpoint_Ia.o B_Spline.o correct_Ip.o \
search_inf_inter.o get_id_inter.o block_norm.o set_id_inter.o\
gmres_normal.o blockgmres_norm.o points_regen.o block_curv.o\
set_id_curv.o gmres_curv.o blockgmres_curv.o get_inter_normal.o\
get_normal_Bspline.o B_Spline_0order.o B_Spline_1order.o \
get_curv_Bspline.o B_Spline_2order.o get_center_coor.o point_projection.o\
get_arc_Bspline.o get_surten_Bspline.o get_intervel_fem.o reset_Icenter.o\
get_intervel_Bspline.o get_centerpoint_Ia.o solve_Laplace_v2.o gmres_Laplace_v2.o\
B_Spline_inter.o corr_I_fluid_center.o\

IFEM: $(OBJ)
	gfortran IFEM $(OBJ) $(LIBS)
.f90.o:
	gfortran -c -g $<
clean:
	rm -rf *.o *.mod IFEM
