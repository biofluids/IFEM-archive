SHELL=/bin/sh
OBJ=main.o initialize.o parseinput.o echoinput.o hypo.o meshgen.o\
    ewdio.o ewdmem.o error.o norm.o connection.o connectionnode.o dyno.o hg.o\
    isatty.o getkey.o getstr.o getstrn.o getint.o getreal.o iei2cray.o\
    form.o equal.o set.o shape.o disk.o quad3d4n.o quad3d8n.o fclear.o\
    vol.o update.o block.o blockf.o blockgmres.o blockgmresf.o hydro.o\
    gmres.o gmresf.o quad2d3n.o quad2d4n.o lenght.o sharp.o cut.o get.o\
    locate.o dim.o liftdrag.o
POSTOBJ=mainpost.o tecplotpost.o dataout.o ewdmem.o error.o ewdio.o fclear.o iei2cray.o getkey.o getstr.o getint.o getreal.o initialize.o meshgen.o parseinput.o postin.o

FFLAGS= -N 80 -O2
CFLAGS= -O2
CF77=f90 $(FFLAGS)
CC=cc
LDR=f90 $(FFLAGS)
COMMAND=xts-mpi 
LIB=

$(COMMAND):	$(OBJ)
	$(LDR) $(LDFLAGS) -o $(COMMAND) $(OBJ) $(LIB)

post: $(POSTOBJ)
	$(LDR) $(LDFLAGS) -o post $(POSTOBJ) $(LIB)

debug:
	$(MAKE) $(COMMAND) FFLAGS='$(FFLAGS) -g' CFLAGS='$(CFLAGS) -g'

app:
	$(MAKE) $(COMMAND) FFLAGS='$(FFLAGS) -Wf"-Ta"' LIB='$(LIB) -lapp'

profile:
	$(MAKE) $(COMMAND) FFLAGS='$(FFLAGS) -g -prof' CFLAGS='$(CFLAGS) -g'

clean:
	rm -f *.o *.l *.T

.SUFFIXES:      .fcm    .pe.o

.f.o:
	$(CF77) -c $<
