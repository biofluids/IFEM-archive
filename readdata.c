#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void readdata_(cstr,offset,origin,array,nbytes,myid)
     int *offset, *origin, *array, *nbytes,*myid;
     char *cstr;
{   
  int fp,i;
  int clen;
    char temp[5];
  int otest,ntest;

  clen=strlen(cstr);

  clen=clen-1;
  printf("clen=%d",clen);
  for (i=0;i<clen+5;i++){
    if (*(cstr+i)=='.') {clen=clen+1;}
  }
  *(cstr+clen)=' ';
  *(cstr+clen+1)=' ';
  
  i=0;
  if (myid==0) {
  do  {
    *(temp+i)=*(cstr+i);
    printf("temp=%s,cstr=%s,\n",temp,cstr);
    i++;
  } while (*(cstr+i)!=' ');
  }
   printf("temp=%s,",temp);
  /* open file */
  
  if ((fp=open("mxyz",O_RDONLY ,0))==-1){
    printf("ewd_open:can't open file mxyz\n");
    exit(1);
  }

  /* lseek */
  otest = (int) lseek(fp,*offset, *origin);
  if (otest!=*offset)  { perror("ewd_lseek"); exit(1); }

  /* read file*/
  ntest = read(fp, array, *nbytes);

  if (ntest != *nbytes) { perror("ewd_read"); exit(1); }  
  close(fp);
  printf("okay");
  return;
}
/***********************************************/
void readmien_(offset,origin,array,nbytes)
     int *offset, *origin, *array, *nbytes;
{
  int i;
  int fp;
  int otest,ntest;

  /* open file */
  if ((fp=open("mien",O_RDONLY ,0))==-1){
    printf("ewd_open:can't open file mien\n");
    exit(1);
  }

  /* lseek */
  otest = (int) lseek(fp,*offset, *origin);
  if (otest!=*offset)  { perror("ewd_lseek"); exit(1); }

  /* read file*/
  ntest = read(fp, array, *nbytes);

  if (ntest != *nbytes) { perror("ewd_read"); exit(1); }  
  close(fp);
  
  return;
}
/***********************************************/
void readmxyz_(offset,origin,array,nbytes)
     int *offset, *origin, *nbytes;
     double *array;
{
  int i;
  int fp;
  int otest,ntest;

  /* open file */
  if ((fp=open("mxyz",O_RDONLY ,0))==-1){
    printf("ewd_open:can't open file mxyz\n");
    exit(1);
  }

  /* lseek */
  otest = (int) lseek(fp,*offset, *origin);
  if (otest!=*offset)  { perror("ewd_lseek"); exit(1); }

  /* read file*/
  ntest = read(fp, array, *nbytes);

  if (ntest != *nbytes) { perror("ewd_read"); exit(1); }  
  close(fp);
  
  return;
}
/***********************************************/
void readmrng_(offset,origin,array,nbytes)
     int *offset, *origin, *array, *nbytes;
{
  int i;
  int fp;
  int otest,ntest;

  /* open file */
  if ((fp=open("mrng",O_RDONLY ,0))==-1){
    printf("ewd_open:can't open file mrng\n");
    exit(1);
  }

  /* lseek */
  otest = (int) lseek(fp,*offset, *origin);
  if (otest!=*offset)  { perror("ewd_lseek"); exit(1); }

  /* read file*/
  ntest = read(fp, array, *nbytes);

  if (ntest != *nbytes) { perror("ewd_read"); exit(1); }  
  close(fp);
  
  return;
}
/***********************************************/
void readrestart_(offset,origin,array,nbytes)
     int *offset, *origin, *nbytes;
     double *array;
{
  int i;
  int fp;
  int otest,ntest;

  /* open file */
  if ((fp=open("data",O_RDONLY ,0))==-1){
    printf("ewd_open:can't open file data\n");
    exit(1);
  }

  /* lseek */
  otest = (int) lseek(fp,*offset, *origin);
  if (otest!=*offset)  { perror("ewd_lseek"); exit(1); }

  /* read file*/
  ntest = read(fp, array, *nbytes);

  if (ntest != *nbytes) { perror("ewd_read"); exit(1); }  
  close(fp);
  
  return;
}
/********************************************/
void datawrite_(i4,i3,i2,i1,offset,origin,array,nbytes)
     char *i4,*i3,*i2,*i1;
     int *offset, *origin, *nbytes;
     double *array;
{
  char name[10];
  char head[10];
  int i;
  int fp;
  int otest,ntest;

  /*head="data.";*/
  strcpy(name,"data.");
  strcat(name,i4);
  strcat(name,i3);
  strcat(name,i2);
  strcat(name,i1);
  printf("writing to %s;\n",name);

  /* open file */
   if ((fp=open(name,O_RDWR | O_CREAT ,384))==-1){
    printf("ewd_open:can't open file data\n");
    exit(1);
    }
  /*fp=open(name,O_RDWR| O_CREAT ,384);*/

  /* lseek */
  otest = (int) lseek(fp,*offset, *origin);
  if (otest!=*offset)  { perror("ewd_lseek"); exit(1); }

  /* read file*/
  ntest = write(fp, array, *nbytes);

  if (ntest != *nbytes) { perror("ewd_read"); exit(1); }  
  close(fp);
  
  return;
}


     
