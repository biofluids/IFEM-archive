/* modified from CRAY to IEEE, by LUCY ZHANG, 1/21/01 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define NUMSIZE 80

void getint_ (pint)
int *pint;
{

	int istat,myid,numproc,ierr;
	extern char cbuffer[];
	char *cnum = calloc(NUMSIZE,sizeof(char)), *cnuminb; int i;

	ierr = MPI_Comm_rank( MPI_COMM_WORLD, &myid);
	ierr = MPI_Comm_size( MPI_COMM_WORLD, &numproc);

	if (myid==0) {
		if (sscanf(cbuffer," %s",cnum)!=EOF) {
			cnuminb = strpbrk(cbuffer, cnum);
			for (i=0;i<strlen(cnum);i++) *(cnuminb+i) = ' ';
			*pint = atoi(cnum);
		}
	} 
	ierr = MPI_Bcast(pint,1,MPI_INT,0,MPI_COMM_WORLD);
	ierr = MPI_Barrier(MPI_COMM_WORLD);

	free(cnum);
	return;
}
