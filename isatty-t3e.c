#include <fortran.h>

/* isatty(3C) wrapper */
int ISATTY(int *fd)
{
	return(_btol(isatty(*fd)));
}
