/* modified from CRAY to IEEE, by LUCY ZHANG, 1/21/01 */
/* modified for single processor */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  /* int istat,myid,numproc,ierr;*/

      char *ckeyinb;
      int i, noeof = 0;
      int ilen =32;
      /*int ilen = strlen(ckey);*/

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

	return(bconvert_(&noeof));
}

