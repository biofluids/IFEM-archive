#include <fortran.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* open a file and return file descriptor */
void EWD_OPEN(str, fd)
_fcd str;
int *fd;
{
	char *cstr; int clen = _fcdlen(str) + 1;
	cstr = malloc(clen);
	cstr = strncpy(cstr, _fcdtocp(str), clen);
	while (clen > 1) {
		clen--;
		if (*(cstr+clen-1) != ' ') {
			if (*(cstr+clen-1) != (char) 0) *(cstr+clen) = (char) 0;
			break;
		}
	}
	*fd = open(cstr,O_RDWR | O_CREAT | O_SYNC | O_RAW, 384);
	/* flags should be passed from the calling program - instead we do this */
	if (*fd <= 0) { 
		printf("ewd_open: opening file as read-only\n");
		*fd = open(cstr, O_RDONLY | O_SYNC | O_RAW);
	}
	if (*fd <= 0) { perror(cstr); exit(1); }
	return;
}

/* close a file given a file descriptor */
void EWD_CLOSE (fd)
int *fd;
{
	if (*fd != 0) close(*fd);
	return;
}

/* write to a file given a file descriptor */
void EWD_WRITE (fd, array, nbytes)
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

/* read from a file given a file descriptor */
void EWD_READ (fd, array, nbytes)
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

/* position a file given a file descriptor */
void EWD_LSEEK (fd, offset, origin)
int *fd, *offset, *origin;
{
	int otest;
	if (*fd == 0) return;
	otest = (int) lseek(*fd, *offset, *origin);
	if (otest != *offset) { perror("ewd_lseek"); exit(1); }
	return;
}
