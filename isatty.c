/* #include <fortran.h> */

/* isatty(3C) wrapper */
int ISATTY(int *fd)
{
	return(isatty(*fd));
}
