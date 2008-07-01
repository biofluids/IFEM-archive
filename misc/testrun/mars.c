/*#include <DTRA/fc.h>*/
extern void mars3d_RPICFD ( int *i4opt, double *dt, int *m4wn, int *n4wn, double *c8wn, double *v8wn, double *a8wn, double *f8wn);
void mars3d_rpicfd_( int *i4opt, double *dt, int *m4wn, int *n4wn, double *c8wn, double *v8wn, double *a8wn, double *f8wn)
{
  mars3d_RPICFD(i4opt, dt, m4wn, n4wn, c8wn, v8wn, a8wn, f8wn);
}
