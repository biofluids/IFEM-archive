OBJ=main.o initialize.o parseinput.o echoinput.o hypo.o meshgen.o\
    ewdio.o ewdmem.o error.o norm.o connection.o connectionnode.o dyno.o hg.o\
    isatty.o getkey.o getstr.o getint.o getreal.o iei2cray.o\
    form.o equal.o set.o shape.o disk.o quad3d4n.o quad3d8n.o fclear.o\
    vol.o update.o block.o blockf.o blockgmres.o blockgmresf.o hydro.o\
    gmres.o gmresf.o quad2d3n.o quad2d4n.o lenght.o sharp.o cut.o get.o\
    locate.o dim.o gmresm.o blockm.o formm.o blockgmresm.o updatex.o\
    lump.o solerror.o add.o defgrad.o velocity.o liftdrag.o readdata.o\
    fcc.o
POSTOBJ=mainpost.o tecplotpost.o dataout.o ewdmem.o error.o ewdio.o fclear.o iei2cray.o getkey.o getstr.o getint.o getreal.o initialize.o meshgen.o parseinput.o postin.o

FC=pgf90 
CC=pgcc
LDR=pgf90 -Mextend
COMMAND=xts-mpi 
LIB=-lmpich
.f.o:
	$(FC) -g -c -Mextend $<
.c.o:
	$(CC) -g -c $<

$(COMMAND):	$(OBJ)
	$(LDR) -o $(COMMAND) $(OBJ) $(LIB)

clean:
	rm -f *.o *.l *.T

.SUFFIXES:      .fcm    .pe.o

