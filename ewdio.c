#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
/* open a file and return file descriptor */
void ewd_open_ (cstr,fd)

char *cstr;
int *fd;
{

  int i;
  int clen;
  char *temp;
  clen=strlen(cstr);
  clen=clen-1;
  for (i=0;i<clen+5;i++){
    if (*(cstr+i)=='.') {clen=clen+1;}
  }

  /* *(cstr+clen)='\0';
  temp="mxyz";
  printf("%s",temp);*/
  *fd=open(cstr,O_RDWR | O_CREAT | O_SYNC ,384);
  if (*fd<=0) {
    printf("ewd_open: opening file as read-only\n");
    *fd = open(cstr, O_RDONLY | O_SYNC );
    perror(cstr);
    exit(1);
  }
  /* if ((*fd=open(cstr,O_RDWR | O_CREAT | O_SYNC ,0))<0) {
    printf("ewd_open: can't open file %s ",cstr);
    exit(1);
    }*/
  
  return;
}

/********************************************************************/

/* close a file given a file descriptor */
void ewd_close_ (fd)
int *fd;
{
  if (*fd != 0) close(*fd);
  return;
}

/******************************************************************/
/* write to a file given a file descriptor */
void ewd_write_ (fd, array, nbytes)
int *fd;
int *array;
int *nbytes;
{
	int ntest;
	if (*fd == 0) return;
	ntest = write(*fd, array, *nbytes);
	if (ntest != *nbytes) { perror("ewd_write"); exit(1); }
	return;
}

/******************************************************************/
/* read from a file given a file descriptor */
void ewd_read_ (fd, array, nbytes)
int *fd;
int *array;
int *nbytes;
{
	int ntest;
	if (*fd == 0) return;

	ntest = read(*fd, array, *nbytes);
	
	if (ntest != *nbytes) { perror("ewd_read"); exit(1); }
	return;
}

/*********************************************************************/
/* position a file given a file descriptor */
void ewd_lseek_ (fd, offset, origin)
int *fd, *offset, *origin;
{
	int otest;
	if (*fd == 0) return;
	otest = (int) lseek(*fd, *offset, *origin);
	if (otest != *offset) { perror("ewd_lseek"); exit(1); }
	return;
}
