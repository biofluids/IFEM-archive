.SUFFIXES: .f90
OBJ= global_constants.o global_simulation_parameter.o \
run_variables.o solid_variables.o \
fluid_variables.o block.o gmres.o blockgmres.o\
ensight_output.o form.o formm.o \
meshgen_fluid.o \
read.o parseinput.o echoinput.o equal.o error.o \
facemap.o fclear.o hg.o hydro.o hypo.o initialize.o \
lenght.o main.o nondimension.o norm.o quad3d4n.o quad3d8n.o \
quad2d3n.o quad2d4n.o \
set.o shape.o sharp.o \
update.o velocity.o vol.o 

IFEM: $(OBJ)
	g95 -o IFEM $(OBJ)
.f90.o:
	g95 -c $<

clean:
	rm *.o *.mod IFEM
