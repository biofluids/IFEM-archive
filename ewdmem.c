#include <malloc.h>

/* malloc(3C) wrapper */
int *MALLOC(int *size)
{
	int *newptr;
	newptr = (int *) malloc((unsigned) *size);
	if (newptr == 0) {
		printf("malloc: null pointer\n");
		/* printf("malloc: _memerror is %d\n", _memerror);*/
		printf("malloc: see /usr/include/malloc.h for more info\n");
		printf("malloc: called with size = %d\n", *size);
		/* tracebk();*/
		exit(1);
	}
	return(newptr);
}

/* realloc(3C) wrapper */
int *REALLOC(int *ptr, int *size)
{
	int *newptr;
	if (*ptr == 0) newptr = (int *) malloc((unsigned) *size);
	else newptr = (int *) realloc((char *) *ptr, (unsigned) *size);
	if (newptr == 0) {
		printf("realloc: null pointer\n");
		/*printf("realloc: _memerror is %d\n", _memerror);*/
		printf("realloc: see /usr/include/malloc.h for more info\n");
		printf("realloc: called with ptr = %d, size = %d\n", *ptr, *size);
		/*tracebk();*/
		exit(1);
	}
	return(newptr);
}

/* free(3C) wrapper                     */
/* 940316 - it finally works...         */
/* used to be free((void *) ptr);       */
/* which frees the cft77 pointer itself */
void FREE(int *ptr)
{
	free((void *) *ptr);
	*ptr = 0;
	return;
}

/* calloc(3C) wrapper */
int *CALLOC(int *nelem, int *elsize)
{
	int *newptr;
	newptr = (int *) calloc((unsigned) *nelem, (unsigned) *elsize);
	if (newptr == 0) { printf("calloc: null pointer\n"); exit(1); }
	return(newptr);
}

/* void *EWDLOC(int *ptr)
{
	printf("ewdloc: pointer value is %d\n", ptr);
	return;
}*/
