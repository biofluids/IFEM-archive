OBJ=main.o initialize.o parseinput_fluid.o echoinput.o hypo.o\
    error.o norm.o hg.o facemap.o meshgen.o\
    nondimension.o\
    form.o equal.o set.o shape.o quad3d4n.o quad3d8n.o fclear.o\
    update.o block.o blockgmres.o hydro.o vol.o gjinv.o\
    gmres.o quad2d3n.o quad2d4n.o lenght.o sharp.o\
    r_input.o fiber5.o link2.o link3.o link4.o io8.o fiber1.o\
    fiber8.o fiber2.o link6.o r_main.o r_sinit.o r_sreadinit.o r_stang.o\
    r_element.o r_jacob.o r_bdpd.o r_stoxc.o r_smaterj.o r_spress.o\
    r_sbpress.o r_sboc.o r_sstrain.o r_spiola.o r_sstif.o r_scauchy.o\
    gaussj.o r_scalfu.o r_scalkup.o r_scalfp.o r_scalkpp.o f_fiber2.o\
    f_fiber1.o fiber3.o disturbance.o fiber13.o fiber10.o fiber9.o fiber11.o\
    r_timefun.o r_load.o r_nodalf.o zibm_tec.f zibm_ensGeo.o zibm_ensFluid.o\
    zfem_tec.o zfem_ensGeo.o zfem_ensFluid.o delta.o calcaccel.o\
    zfem_ensCase.o zibm_ensCase.o movepoints.o rkpmshape3d.o correct3dl.o\
    delta_nonuniform.o r_print.o r_main2.o
POSTOBJ=mainpost.o tecplotpost.o dataout.o ewdmem.o error.o ewdio.o fclear.o iei2cray.o getkey.o getstr.o getint.o getreal.o initialize.o meshgen.o parseinput.o postin.o

FC=pgf90
CC=pgcc
LDR=pgf90
COMMAND=ifem 
LIB=

$(COMMAND):	$(OBJ)
	$(LDR) -o $(COMMAND) $(OBJ) $(LIB)

clean:
	rm -f *.o *.l *.T

.SUFFIXES:      .fcm    .pe.o

