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
	fprintf (f, "\"X\" \"Y\" \"U\" \"V\" ");
  }
  else if(*nsd == 3) {
    fprintf (f, "\"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" "); 
  }
  if(*ndf - *nsd == 1) {
	fprintf (f, "\"P\" \"PHI\" ");
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
  printf("Writing file %d...\n", *count);
  if ((*count) == 0) {
   // fprintf(f, "ZONE T=\"%d\", N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON\n",
   //         *count, *nn, *ne);
    fprintf(f, "ZONE T=\"%d\", N=%d, E=%d, F=FEPOINT, ET=TETRAHEDRON\n",
            *count, *nn, *ne);
  }

  else {
    fprintf(f, "ZONE T=\"%d\", N=%d, E=%d, F=FEPOINT, ", *count, *nn, *ne);
    //fprintf(f,"ET=TETRAHEDRON, D=(1,2,3,FECONNECT)\n");
    fprintf(f,"ET=TETRAHEDRON, D=(FECONNECT)\n");
  }
  for (j = 0; j < *nn; j++) {
   // if ((*count) == 0) {
   //   for (i = 0; i < *nsd; i++) {
   //     fprintf(f, "%9.5lf ", *xnout);
   //     xnout++;
    //  }
   // }
	for (i = 0; i < *ndf+4; i++) {
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

void HEADEROUT2(_fcd title, int *nsd, int *ndf) {
  char *ctitle = _fcdtocp(title);
  int i, titlelen = _fcdlen(title);
  FILE *f = fopen("tecoutmesh.dat", "w");
  printf("title length = %d\n", strlen(ctitle));
  fprintf(f, "TITLE = \"%s\"\n", ctitle);
  fprintf(f, "VARIABLES = ");
  if(*nsd == 2) {
	fprintf (f, "\"X\" \"Y\" ");
  }
  else if(*nsd == 3) {
    fprintf (f, "\"X\" \"Y\" \"Z\" "); 
  }
  
  fprintf(f, "\n");
  fclose(f);
  return;
}


void POSTOUT2(int *count, int *xn, int *nsd, int *ndf, 
	     int *nn, int *ne, int *nen) {

  int i,j;
  int *xnout = xn;
  //int *dnout = dn;
  FILE *f = fopen("tecoutmesh.dat", "a");
  printf("Writing file %d...\n", *count);
  if ((*count) == 0) {
    fprintf(f, "ZONE T=\"%d\", N=%d, E=%d, F=FEPOINT, ET=brick\n",
            *count, *nn, *ne);
  }
  else {
    fprintf(f, "ZONE T=\"%d\", N=%d, E=%d, F=FEPOINT, ", *count, *nn, *ne);
    //    fprintf(f,"ET=brick, D=(1,2,3,FECONNECT)\n");
    fprintf(f, "ET=brick, D=(FECONNECT)\n");
  }
  for (j = 0; j < *nn; j++) {
    //if ((*count) == 0) {
      for (i = 0; i < *nsd; i++) {
	fprintf(f, "%9.5lf ", *xnout);
	//printf("%9.5lf ",*xnout);
	xnout++;
      }
      //printf("\n");
      fprintf(f, "\n");
      //    }
    // else{
      //  for (i = 0; i < *nsd; i++) {
	//	fprintf(f, "%9.5lf ", *dnout);
	//	dnout++;
	//      }
      //    fprintf (f, "\n");
      //    }  
  }
  fprintf (f, "\n");
  fclose(f);
  return;
}

void MESHOUT2(int *ien, int *nen, int *ne) {
  int i,j;
  int *iout = ien;
  FILE *f = fopen("tecoutmesh.dat", "a");
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







