#include <fortran.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void GETSTRN (str)
_fcd str;
{
      int istat,myid,numproc,ierr;
	extern char cbuffer[];
	char *cstr = _fcdtocp(str), *cstrinb; int i, ilen = _fcdlen(str);
	ierr = MPI_Comm_rank( MPI_COMM_WORLD, &myid);
	ierr = MPI_Comm_size( MPI_COMM_WORLD, &numproc);

	if (myid==0) {
		if (sscanf(cbuffer," %[^\n]",cstr)!=EOF) {
			cstrinb = strpbrk(cbuffer, cstr);
			for (i=0;i<strlen(cstr);i++) *(cstrinb+i) = ' ';
		}
		for (i=strlen(cstr);i<ilen;i++) *(cstr+i) = ' ';
	}

      ierr = MPI_Bcast(cstr,ilen,MPI_CHAR,0,MPI_COMM_WORLD);
	ierr = MPI_Barrier(MPI_COMM_WORLD);

	return;
}
