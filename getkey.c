/* modified from CRAY to IEEE, by LUCY ZHANG, 1/21/01 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#define LINESIZE 80 
/*#define TRUE 0xff
  #define FALSE 0 */

int nstack = -1;
FILE *fstack[10];
char cbuffer[LINESIZE];
/* get a line from stdin and read first string as key */
int getkey_(ckey)
     char *ckey;
{
      int istat,myid,numproc,ierr;

      char *ckeyinb;
      int i, noeof = 0;
      int ilen =32;
      /*int ilen = strlen(ckey);*/

      	ierr = MPI_Comm_rank( MPI_COMM_WORLD, &myid);
	ierr = MPI_Comm_size( MPI_COMM_WORLD, &numproc);

	if (myid==0) {
	  if (nstack == -1) { fstack[0] = stdin; nstack = 0; }
	  /* clean the buffer */
	  for (i=0;i<LINESIZE;i++) *(cbuffer+i) = ' ';
	  for (i=0;i<ilen;i++) *(ckey+i) = ' ';
	  /* prompt */
	  if (isatty(fileno(fstack[nstack]))) printf("> ");
	  /* get a line */
	  if (fgets(cbuffer,LINESIZE,fstack[nstack])!=0) {
	    if (sscanf(cbuffer,"%s",ckey)!=EOF) {
				/* find key in buffer and erase it */
	      ckeyinb = strpbrk(cbuffer, ckey);
	      for (i=0;i<strlen(ckey);i++) *(ckeyinb+i) = ' ';
	    }
	    /* fill rest of key with spaces, so that .eq. works */
	    for (i=strlen(ckey);i<ilen;i++) *(ckey+i) = ' ';
	    noeof = 1;
	  }
	}

	/*printf("length of key= %d\n",strlen(ckey));*/
	ierr = MPI_Bcast(ckey,ilen,MPI_CHAR,0,MPI_COMM_WORLD);
	ierr = MPI_Bcast(&noeof,1,MPI_INT,0,MPI_COMM_WORLD);
	ierr = MPI_Barrier(MPI_COMM_WORLD);

	return(bconvert_(&noeof));
}

