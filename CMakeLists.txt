Project(ifem Fortran)
set(ifem_sources global_constants.f90 global_simulation_parameter.f90 run_variables.f90
r_common.f90 fluid_variables.f90 solid_variables.f90 deformation_gradient.f90
delta_nonuniform.f90 ensight_output.f90 form.f90 formm.f90 block.f90 blockgmres.f90 gmres.f90
gmres_fluid_ale.f90 gmres_fluid_mesh.f90 meshgen_solid.f90 meshgen_fluid.f90 solid_fem_BC.f90
read.f90 parseinput.f90 correct.f90  echoinput.f90 equal.f90 error.f90
facemap.f90 fclear.f90 gaussj.f90 gjinv.f90 hg.f90 hydro.f90 hypo.f90 initialize.f90
length.f90 main.f90 nondimension.f90 norm.f90 quad3d4n.f90 quad3d8n.f90
quad2d3n.f90 quad2d4n.f90 r_bdpd_curr.f90 r_bdpd_init.f90 r_element.f90
r_sbpress.f90 r_jacob.f90 r_load.f90 r_nodalf.f90 r_sboc.f90
r_scal.f90 r_scauchy.f90 r_smaterj.f90 r_spiola.f90 r_spiola_viscous.f90 r_spiola_elastic.f90
r_spress.f90 r_sreadinit.f90 r_sstif.f90 r_sstrain.f90 r_stang.f90 r_stoxc.f90
r_timefun.f90 rkpmshape2d.f90 rkpmshape3d.f90 set.f90 shape.f90 sharp.f90 solid_solver.f90
solid_update.f90 update.f90 velocity.f90 vol.f90) 
add_executable(IFEM ${ifem_sources})
