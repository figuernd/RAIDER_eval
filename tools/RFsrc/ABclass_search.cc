#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define  TRUE  1
#define  FALSE  0
#ifndef  EXIT_FAILURE
  #define  EXIT_FAILURE  -1
#endif
#ifndef  EXIT_SUCCESS
  #define  EXIT_SUCCESS  0
#endif



typedef struct bacs{
  int name;
  int id;
  char title[20];
  int signaln;
} BAC_Type, * BAC_Ptr;


typedef struct inMergleng *inMLeng_Ptr;

typedef struct inMergleng{

  int length;
  int copies;
  inMLeng_Ptr next;
  inMLeng_Ptr prev;
  
} inMLeng_Type;


typedef struct inMergname *inMerg_Ptr;

typedef struct inMergname{

  int name;
  inMLeng_Ptr Length;
  int length;
  int copies;
  int depend;
  inMerg_Ptr next;
  inMerg_Ptr prev;
  
} inMerg_Type;

typedef struct inpname *inName_Ptr;

typedef struct inpname{

  int name;
  int length;
  int copies;
  inName_Ptr next;
  inName_Ptr prev;
  
} inName_Type;

typedef struct rcls{
  
  int name;
  int clname;
  
} Rep_Type, * Rep_Ptr;

typedef struct rcela *Exc_Ptr;

typedef struct rcela{
  int name;
  int length;
  int copies;
  inMerg_Ptr Ereps;
  inName_Ptr clname;
} Exc_Type;

typedef struct buf{
  int clcoor;
 } Buf_Type,* Buf_Ptr;

typedef struct inpcls *inPool_Ptr;

typedef struct inpcls{

  int name;
  int length;
  int copies;
  int num;
  inPool_Ptr next;
  inPool_Ptr prev;
  
}inPool_Type;


typedef struct pcls * Pool_Ptr;

typedef struct pcls{

  int name;
  int length;
  int copies;
  int min_name;
  int min_length;
  int num;
  inMerg_Ptr Ereps;
  inPool_Ptr pool;
  
}Pool_Type;


Pool_Ptr  Clist;
Buf_Type Blist[200];
bacs Baclist[1000000];
int  Baclen;

FILE *  File_Open  (const char *, const char *);
void *  Safe_calloc  (size_t,size_t);
inPool_Ptr inPalloc(void);
int addInPool(int ,int ,int , int, int , int);
void Read_Rep_map(const char * );
void classPrint(int);
void RepeaPrint(int);
inName_Ptr inNalloc(void);
void ABRepeaPrint(int , const char *, const char *);
void BlastaRep_map(const char *);



int  Binary_Bsearch ( bacs [], int , int , int );
int  Binary_Bsearch2 ( bacs [], int , int , int );

main(int argc, char * argv [])
{
   FILE  * fp,*fp1,*fp2,*Q;
   FILE * f1, *f2,*f3,*f5,*f6,*f8,*f9;

   int Data_Len;
   int Data_Len2;
   int Data_Len3;
   
      
   Data_Len3= 300000;
    
   Clist=(Pool_Ptr) Safe_calloc (Data_Len3, sizeof (Pool_Type)); 
     
   Read_Rep_map(argv [argc - 4]);
  
   BlastaRep_map(argv [argc - 3]);
   
   ABRepeaPrint(Data_Len3,argv [argc - 1],argv [argc - 2]);

   free(Clist);
   
}



inPool_Ptr inPalloc(void)
{
  return (inPool_Ptr)malloc(sizeof(inPool_Type));
}


int addInPool(int i1,int i, int name, int length,int copies,int es)
{
  inPool_Ptr Temp,Tpool,T1pool, Prev, Npool;
  
  inName_Ptr Tname;
  int j1,stops;
  if(Clist[i1].pool==NULL)
    {
   
      Npool=inPalloc();
      Npool->name=name;
      Npool->length=length;
      Npool->copies=copies;
      Npool->prev=NULL;
      Npool->next=NULL;
      Clist[i1].pool=Npool;
    
      Prev= Npool;
           
    }
  else 
    if (es)
      {
		
        stops=1;
         
	  for (Tpool=Clist[i1].pool;(Tpool!=NULL && stops); Tpool=Tpool->next){
   
            Prev=Tpool->prev;
	    if (Tpool->name==name) stops=0;
            if (Tpool->name> name) {
	      Npool=inPalloc();
              if(Prev!=NULL){
              Prev->next=Npool;

	      Npool->name= name;
              Npool->length=length;
	      Npool->copies=copies;
	      Npool->prev=Tpool->prev;
	      Npool->next=Tpool;
	      Tpool->prev=Npool;
	      }
              else{
		Npool->prev=NULL;
                Npool->name= name;
                Npool->length=length;
		Npool->copies=copies;
                Npool->next=Tpool;
                Tpool->prev=Npool;
                Clist[i1].pool=Npool;
	      }
              
	      stops=0;
	    }
           Prev=Tpool;
	  }
	  if (Tpool==NULL && stops) {
	    Npool=inPalloc();
            if(Prev!=NULL)
            Prev->next=Npool;
	    Npool->name=name;
            Npool->length=length;
	    Npool->copies=copies;
	    Npool->prev=Prev; 
            Npool->next=NULL; 
	  }
  	
      }
    else {
       
   
      for(T1pool=Clist[i].pool;(T1pool !=NULL ); T1pool=T1pool->next){
	stops=1;
	
	for (Tpool=Clist[i1].pool;(Tpool!=NULL && stops); Tpool=Tpool->next){
	  Prev=Tpool->prev;
	  if (Tpool->name==T1pool->name) stops=0;
	  if (Tpool->name>T1pool->name) {
	    Npool=inPalloc();
            if(Prev!=NULL){
            Prev->next=Npool;
	    Npool->name=T1pool->name;
            Npool->length=T1pool->length;
            Npool->copies=T1pool->copies;
	    Npool->prev=Tpool->prev;
	    Npool->next=Tpool;
	    Tpool->prev=Npool;
	    } 
	    else{
		Npool->prev=NULL;
                Npool->name=T1pool->name;
                Npool->copies=T1pool->copies;
                Npool->length=T1pool->length;
		Npool->next=Tpool;
                Tpool->prev=Npool;
                Clist[i1].pool=Npool;
	      }
	    stops=0;
           
	  }
          Prev=Tpool;
	}
	if (Tpool==NULL && stops) {
	  Npool=inPalloc();
          if(Prev!=NULL)
           Prev->next=Npool;
	  Npool->name=T1pool->name;
	  Npool->length=T1pool->length;
	  Npool->copies=T1pool->copies;
	  Npool->prev=Prev; 
	  Npool->next=NULL; 
	}
      }
    }
  /*   fprintf (stderr," %10d ",i1);
  for (Tpool=Clist[i1].pool;(Tpool!=NULL); Tpool=Tpool->next){
  fprintf (stderr," %10d ",Tpool->name);
  } 
  fprintf (stderr," \n ");
  */
}

void Read_Rep_map(const char * Filename)
{ 
  FILE *fp;
  char c[6];
  int name, coord, length, scores;
 
  fp = File_Open (Filename, "r");

   while (fscanf(fp," %d %d %d %d \n",&name,&coord,&length,&scores)!=EOF ){
    
    addInPool(name,0,coord,length,scores,1);
    if(Clist[name].name==0) Clist[name].name=coord;
    
    
  }
  fclose(fp);
}
  
/*Blasta searching!*/

void BlastaRep_map(const char * Filename)
{ 
  FILE *Q;
  int i,j,ib,pleng,i1,i2,j1,cb,cb1,es,stopsgn,bsign;
  int stops;
  inPool_Ptr Tpool, Tname;;
  Q = File_Open(Filename,"r") ;
  i=0;j=0;
  pleng=0;
  cb=0;
  es=0;
  /*fprintf (stderr, " We are in BlastaRep_map\n");*/
 while (fscanf(Q,"%d %d\n",&i,&j)!=EOF ){
    while ( Clist[i].name <0 && Clist[i].name!=-1)
      i=- (Clist[i].name/10 +1);
    
    while(Clist[j].name <0 && Clist[j].name!=-1)
      j=-(Clist[j].name/10 +1);

    /*fprintf (stderr, " \n %10d %10d \n",i,j);
    for(Tpool=Clist[j].pool;(Tpool!=NULL); Tpool=Tpool->next){
      fprintf (stderr, " %10d ",Tpool->name);
      }
    */

    if ( i!=j ){
      /*fprintf (stderr, " \n %10d %10d \n",i,j);*/
      if (i<j){
	addInPool(i,j,0,0,0,es);
  
	Clist[j].name=-(i+1)*10;
      }
      if (i>j){
	addInPool(j,i,0,0,0,es);
  
	Clist[i].name=-(j+1)*10;
      }
    }
 }
 fclose(Q);
}

void ABRepeaPrint( int Data_Len3,const char * Filename1, const char * Filename2) 
{
  FILE *f3,*f8 ,*f9;
  int i,j,k,z,j4,j1,ci,ci1,cr1,cr2;
  inPool_Ptr T1pool,Tpool,Prev;
  int bac_num;
  char tc1[20];
  int loc_coord;
  int zero=0;
 

 f9=File_Open(Filename2,"r");
  i=0;
  while (fscanf(f9,"%s %d\n",&tc1,&cr2)  != EOF) {
    strcpy(Baclist[i].title, tc1);
    Baclist[i].name=cr2;
    Baclist[i].signaln=1; 
    i++;
   }
   fclose(f9);
   Baclen=i-1;
  i=0;
 i=0;j=0;
  while ( i<Data_Len3)
    {
      if(Clist[i].name >0){
	j++;
		
	  ci=0;
	 for(Tpool=Clist[i].pool;(Tpool!=NULL); Tpool=Tpool->next){
	   ci++;
           Prev=Tpool;
	   
           bac_num= Binary_Bsearch( Baclist, zero, Baclen,Tpool->name );
           loc_coord=Tpool->name-Baclist[bac_num].name +1;

           printf( ">Class %10d ",j);
           printf( "%10d %10d %10d %20s %10d\n",Tpool->name,Tpool->length,Tpool->copies,Baclist[bac_num].title,loc_coord);

	  
	 }
	 if (ci==1) {
           fprintf(f3, ">Class %10d ",j);
	   /*  fprintf(f3, "%10d ",j);
	       fprintf(f3, "%10d %10d \n",Prev->name+Clist[i].length,Clist[i].length);*/
	  
	}
	
      }
      
      i++;
      }  
  /*fprintf(f3, ">The end. ");*/
  /*   printf ("The number of the repeat classes  %10d \n",j); 
       fclose(f3);  */   
}

int  Binary_Bsearch ( bacs Blist[], int Lo, int Hi, int Val)

/* Find and return  i  such that  Lo <= i <= Hi  and
*  A [i].name <= Val < A [i + 1].name .  If no such value return  -1 . */

  {
   int  i;
  
   /*printf ("%10d %10d %10d\n ", Lo, Hi,i);*/
   while  (Lo <= Hi)
     {
       
      i = (Lo + Hi) / 2;
      /*printf ("%10d %10d %10d %10d\n ", Lo, Hi,i, Blist[i].name);*/
      if  (Blist[i].name-10 <= Val && Val < Blist[i + 1].name-10)
          return  i;
      else if  (Val < Blist [i].name)
          Hi = i - 1;
        else
          Lo = i + 1;
     }

   return  -1;
  }

int  Binary_Bsearch2 ( bacs Blist[], int Lo, int Hi, int Val)

/* Find and return  i  such that  Lo <= i <= Hi  and
*  A [i].name <= Val < A [i + 1].name .  If no such value return  -1 . */

  {
   int  i;
  
   /*printf ("%10d %10d %10d\n ", Lo, Hi,i);*/
   while  (Lo <= Hi)
     {
       
      i = (Lo + Hi) / 2;
      /*printf ("%10d %10d %10d %10d\n ", Lo, Hi,i, Blist[i].name);*/
      if  (Blist[i].name+10 <= Val && Val < Blist[i + 1].name+10)
          return  i;
      else if  (Val < Blist [i].name+10)
          Hi = i - 1;
        else
          Lo = i + 1;
      }

   return  -1;
  }

FILE *  File_Open  (const char * Filename, const char * Mode)

/* Open  Filename  in  Mode  and return a pointer to its control
*  block.  If fail, print a message and exit. */

  {
   FILE  *  fp;

   fp = fopen (Filename, Mode);
   if  (fp == NULL)
       {
        fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
        exit (EXIT_FAILURE);
       }

   return  fp;
  }
void *  Safe_calloc  (size_t N, size_t Len)

/* Allocate and return a pointer to an array of  N  elements of
*   Len  bytes each.  All set to 0.  Exit if fail. */

  {
   void  * P;

   P = calloc (N, Len);
   if  (P == NULL)
       {
        fprintf (stderr, "ERROR:  calloc failed\n");
        exit (EXIT_FAILURE);
       }

   return  P;
  }






















 



















