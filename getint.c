/* modified from CRAY to IEEE, by LUCY ZHANG, 1/21/01 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NUMSIZE 80

void getint_ (pint)
int *pint;
{

	/*int istat,myid,numproc,ierr;*/
	extern char cbuffer[];
	char *cnum = calloc(NUMSIZE,sizeof(char)), *cnuminb; int i;

		if (sscanf(cbuffer," %s",cnum)!=EOF) {
			cnuminb = strpbrk(cbuffer, cnum);
			for (i=0;i<strlen(cnum);i++) *(cnuminb+i) = ' ';
			*pint = atoi(cnum);
		}


	free(cnum);
	return;
}
