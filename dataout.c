/***************/
/* dataout.c   */
/* G. Wagner   */
/***************/

#include <fortran.h>
#include <stdio.h>

void HEADEROUT(_fcd title, int *nsd, int *ndf) {
  char *ctitle = _fcdtocp(title);
  int i, titlelen = _fcdlen(title);
  FILE *f = fopen("tecout.dat", "w");
  printf("title length = %d\n", strlen(ctitle));
  fprintf(f, "TITLE = \"%s\"\n", ctitle);
  fprintf(f, "VARIABLES = ");
  if(*nsd == 2) {
	fprintf (f, "\"X\" \"Y\" ");
  }
  else if(*nsd == 3) {
    fprintf (f, "\"X\" \"Y\" \"Z\" "); 
  }
  if(*ndf - *nsd == 1) {
	fprintf (f, "\"PHI\" ");
  }
  fprintf(f, "\n");
  fclose(f);
  return;
}

void POSTOUT(int *count, int *xn, int *dd, int *nsd, 
             int *ndf, int *nn, int *ne, int *nen) {
  int i,j;
  int *xnout = xn, *ddout = dd;
  FILE *f = fopen("tecout.dat", "a");
  *ndf = 0;
  printf("Writing file %d...\n", *count);
  if ((*count) == 0) {
    fprintf(f, "ZONE T=\"%d\", N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON\n",
            *count, *nn, *ne);
  }
  else {
    fprintf(f, "ZONE T=\"%d\", N=%d, E=%d, F=FEPOINT, ", *count, *nn, *ne);
    fprintf(f,"ET=TETRAHEDRON, D=(1,2,3,FECONNECT)\n");
  }
  for (j = 0; j < *nn; j++) {
    if ((*count) == 0) {
      for (i = 0; i < *nsd; i++) {
        fprintf(f, "%9.5lf ", *xnout);
        xnout++;
      }
    }
	for (i = 0; i < *ndf+1; i++) {
	  fprintf(f, "%9.5lf ", *ddout);
	  ddout++;
	}
    fprintf (f, "\n");
  }
  fprintf (f, "\n");
  fclose(f);
  return;
}

void MESHOUT(int *ien, int *nen, int *ne) {
  int i,j;
  int *iout = ien;
  FILE *f = fopen("tecout.dat", "a");
  for(j = 0; j < *ne; j++) {
    for(i = 0; i < *nen; i++) {
      fprintf(f, "%d ", *iout);
	  iout++;
    }
    fprintf(f, "\n");
  }
  fclose(f);
  return;
}







