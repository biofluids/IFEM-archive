/* modified from CRAY to IEEE, by LUCY ZHANG, 1/21/01 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void getstr_ (cstr)
     char *cstr;
{
	/*int istat,myid,numproc,ierr;*/
	extern char cbuffer[];
	/*char *cstr = _fcdtocp(str), *cstrinb; */
	char *cstrinb;
	int i, ilen = 8;


	  *cstr = (char) 0;
	  if (sscanf(cbuffer,"%s",cstr)!=EOF) {
	    cstrinb = strpbrk(cbuffer, cstr);
	    for (i=0;i<strlen(cstr);i++) *(cstrinb+i) = ' ';
	  }
	  for (i=strlen(cstr);i<ilen;i++) *(cstr+i) = ' ';

	return;
}
