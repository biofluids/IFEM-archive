OBJ=main.o initialize.o parseinput.o parseinput_fluid.o echoinput.o hypo.o\
    ewdio.o ewdmem.o error.o norm.o hg.o facemap.o meshgen.o\
    isatty.o getkey.o getstr.o getint.o getreal.o nondimension.o\
    form.o equal.o set.o shape.o disk.o quad3d4n.o quad3d8n.o fclear.o\
    update.o block.o blockgmres.o hydro.o vol.o gjinv.o\
    gmres.o quad2d3n.o quad2d4n.o lenght.o sharp.o cut.o get.o\
    fcc.o readdata.o \

OBJ=block.o blockgmres.o cut.o delta_nonuniform.o disk.o disturbance.o echoinput.o 
    ewdio.o ewdmem.o equal.o error.o \
    facemap.o fclear.o fcc.o form.o fiber1.o fiber8.o fiber2.o fiber5.o fiber13.o 
    fiber10.o fiber9.o fiber11.o f_fiber1.o fiber3.o f_fiber2.o gaussj.o get.o 
    getkey.o getstr.o getint.o getreal.o gjinv.o gmres.o meshgen.o\
    hypo.o hg.o hydro.o link2.o link3.o link4.o link6.o initialize.o io8.o 
    isatty.o lenght.o main.o norm.o nondimension.o\
    parseinput.o parseinput_fluid.o quad3d4n.o quad3d8n.o quad2d3n.o quad2d4n.o 
    r_element.o r_jacob.o r_bdpd.o r_stoxc.o r_smaterj.o r_spress.o\
    r_sbpress.o r_sboc.o r_sstrain.o r_spiola.o r_sstif.o r_scauchy.o r_main.o 
    r_sinit.o r_sreadinit.o r_stang.o\
    r_scalfu.o r_scalkup.o r_scalfp.o r_scalkpp.o r_input.o readmain.o \
    r_timefun.o r_load.o r_nodalf.o\ readdata.o set.o shape.o \
    sharp.o update.o vol.o \
    

    readmain.o r_input.o fiber5.o link2.o link3.o link4.o io8.o fiber1.o\
    fiber8.o fiber2.o link6.o r_main.o r_sinit.o r_sreadinit.o r_stang.o\
    r_element.o r_jacob.o r_bdpd.o r_stoxc.o r_smaterj.o r_spress.o\
    r_sbpress.o r_sboc.o r_sstrain.o r_spiola.o r_sstif.o r_scauchy.o\
    gaussj.o r_scalfu.o r_scalkup.o r_scalfp.o r_scalkpp.o f_fiber2.o\
    f_fiber1.o fiber3.o disturbance.o fiber13.o fiber10.o fiber9.o fiber11.o\
    r_timefun.o r_load.o r_nodalf.o\
    delta_nonuniform.o
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

