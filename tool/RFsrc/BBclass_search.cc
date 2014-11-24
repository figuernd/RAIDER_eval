#include  <stdio.h>
#include  <stdlib.h>
#include  <iostream>
#include  <iomanip>
#include  <math.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>
#include  <unistd.h>


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


Exc_Ptr Elist;
Rep_Ptr Rlist;
Pool_Ptr  Clist;
Buf_Type Blist[200];
bacs Baclist[1000000];
int  Baclen;

FILE *  File_Open  (const char *, const char *);
void *  Safe_calloc  (size_t,size_t);
void reverse (char *);
void itoa(int , char *);
inPool_Ptr inPalloc(void);
int addInPool(int ,int ,int );
void Rep_map(int, int);
void classPrint(int);
/*void RepeaPrint(int,const char *);*/
void RepeaPrint(int);
inName_Ptr inNalloc(void);
int addInName(int ,int);

int  Binary_Esearch2 ( int , int , int );
int  Binary_Esearch ( int , int , int );
int  Binary_Search ( int , int , int );
int  Binary_Bsearch ( bacs [], int , int , int );
int  Binary_Bsearch2 ( bacs [], int , int , int );

main(int argc, char * argv [])
{
  FILE  * fp,*fp1,*fp2,*Q;
  FILE * f1, *f2,*f3,*f5,*f6,*f8,*f9;
   double  Max_rep;
   int  i,j, n, ch,k,z, errflg,i1,i2,pleng;
   int  cr1,cr2, num1, leng1,new_cur;
   int  fr1,fr2, num2, leng2,buff,zero,cur_cr2,ib;
   int diff1, diff2, cur_gap,tandem1,res_temp,long_temp;
   int diff12,diff3,Ident;
   int fsign,csign,long_leng,merge,classr1,classr2,classrn;
   int classr1_leng,cur_cr1,cur_leng;
   int exact, subs, tandem,res_leng,Idelta,Igap;
   int Cur_class1,Cur_class2;
   int Cur_cbac2,Cur_cbacsof2,right_cr2;
   int Cur_rep1,Cur_rep2, Cur_rep1_cl,cur_e;
   int ci, ci1, ci2,csign1,csign2;
   int cb, cb1,cb2,rsign,bsign,stopsgn;
   int Cur_fbac1,Cur_fbacsof,Cur_cbac1, Cur_cbacsof,right_fr1,right_cr1;
   char tc1[20];
   double delta;
   char s;
   inName_Ptr Tname;
   int Lo,Hi,Lo2,Hi2;
   int  rep_len;
   int e_num,e_cur;
   inMerg_Ptr ETname;
   inMLeng_Ptr ETname1;
   int Data_Len;
   int Data_Len2;
   int Data_Len3;
   char *Nstr;
   char *c;
   char  * p;
   Lo=0;

  Data_Len=1000000;
  Data_Len2= Data_Len;
  Data_Len3= Data_Len/3;
  
  
  Rlist = (Rep_Ptr) calloc (Data_Len, sizeof (Rep_Type));
  Elist= (Exc_Ptr) calloc (Data_Len2, sizeof (Exc_Type));
  Clist=(Pool_Ptr) Safe_calloc (Data_Len3, sizeof (Pool_Type)); 
  /*printf ("Allocation completed\n");*/
   
  
   //  Parse the argument list using "man 3 getopt"
   errflg=0;
   optarg = NULL;
   Max_rep=2.;
   zero=0;
  while  (! errflg && ((ch = getopt (argc, argv, "r:o:")) != EOF))
    switch  (ch)
        {
	case 'r' :
         Igap=(int) strtol (optarg, & p, 10);
	 if  (p == optarg )
              {
               fprintf (stderr, "ERROR:  Illegal min repeat length \"%s\"\n",
                        optarg);
               exit (-1);
              }
              optarg = NULL;
	      break;
        case 'o' :
         Ident=(int) strtol (optarg, & p, 10);
	 if  (p == optarg )
              {
               fprintf (stderr, "ERROR1:  Illegal identity value \"%s\"\n",
                        optarg);
               exit (-1);
              }
	 if(Ident<0 || Ident>100)
             {
               fprintf (stderr, "ERROR2:  Illegal identity value \"%s\"\n",
                        optarg);
               exit (-1);
              }
	      optarg = NULL;
	      break;
      case  '?' :
	fprintf (stderr, "Unrecognized option -%c", optopt);
      default :
	errflg++;
      }
  

  if  (argc - optind != 3)
    {
      fprintf (stderr, "USAGE: \n", argv [0]);
      fprintf (stderr, "  # \n");
      exit (EXIT_FAILURE);
    }
  /*fprintf (stderr, "Open all.coord.nodu\n" );*/
  f1 = File_Open (argv [argc - 2],"r") ;
  i=0;
  while (fscanf(f1,"%d\n",&cr1)!=EOF)
    {
      Rlist[i].name=cr1;
      Rlist[i].clname=0;
      i++;
    }
  Hi=i-1;
  /*fprintf (stderr, "Open seq.coord= OK\n" );*/
  fp2=File_Open(argv [argc - 1],"r");
  i=0;
  while (fscanf(fp2,"%s %d\n",&tc1,&cr2)  != EOF) {
    strcpy(Baclist[i].title, tc1);
    Baclist[i].name=cr2;
    Baclist[i].signaln=1; 
    i++;
   }
   fclose(fp2);
   Baclen=i-1;
  
  i=0;
  
  fp = File_Open (argv [argc - 3], "r");

  /*f2 = fopen ("SPLT.merge","a") ;*/
  fr1=0;
  fr2=0;
  num1=0;
  leng1=0;
  fsign=1;
 
  exact=0;
  subs=0;
  res_leng=0;
  merge=0;
  e_num=1;
  while (fscanf(fp,"%d %d %d %d \n",&cr1,&cr2,&leng2,&num2) != EOF) {
       Cur_fbac1= Binary_Bsearch( Baclist, zero, Baclen, fr1);
    
       Cur_cbac2= Binary_Bsearch( Baclist, zero, Baclen, cr2);
       Cur_cbac1= Binary_Bsearch( Baclist, zero, Baclen, cr1);
       Cur_cbacsof= Binary_Bsearch2( Baclist, zero, Baclen, cr1+leng2);
     
     if(subs || merge ){
       Cur_fbac1=Binary_Bsearch( Baclist, zero, Baclen, classr1);
      }
     diff1=cr1-fr1;
     
     diff2=-cr1-leng2+fr1 +leng1;
     diff12=cr2 -fr2>0? cr2 -fr2-leng1:fr2-cr2-leng2;
     diff3=fr1+leng1-cr1;
  
    if ( (diff1>=0 && (diff2>=0 || Ident*diff1<=(100-Ident)*diff3 || -Ident*diff2<=(100-Ident)*diff3  )) &&((num1 !=num2) &&((cr2<cr1 &&cr2<classr1 && cr2>0 && !fsign) ||(classr1_leng+classr1>cr1+leng2 && !fsign)||(cr2<0 && !fsign)||(cr2>cr1 && !fsign)||(cr1 -classr1 +1<Igap && !fsign) || (fsign))) &&(Cur_fbac1==Cur_cbacsof)){
      
      if ((fsign && fr1)||merge){
         
         
	 if (merge){
	    subs=merge;
            merge=0;
          
                      
	 }
	 else{
	  
	   if (diff1==0) {
	     if(num2<num1){ 
	       num1=num2;
	     }
	   }
           
	   classr1=fr1;
	   classr2=fr2;
       
           classr1_leng=leng1;

           Cur_rep1_cl=Binary_Search (Lo,Hi,fr1);

	   Cur_rep2=Binary_Search (Lo,Hi,fr2);
           

	   Elist[e_num].name= classr1;
	  
           addInName(e_num,classr1);
           
	    if (Rlist[Cur_rep1_cl].clname==0){
	       Rlist[Cur_rep1_cl].clname=e_num;
               e_cur=e_num;
	      
	     }
             else {
	       e_cur= Rlist[Cur_rep1_cl].clname;

              addInName(e_cur,classr1);

	       

	       }
                      
	       ci2=0; csign2=0;
  
	    cur_e=Rlist[Cur_rep2].clname;
            if (cur_e !=0){
	      new_cur=Binary_Search(Lo,Hi,Elist[cur_e].name) ;
	                
	      if (Elist[cur_e].name !=0){
		e_cur=Rlist[new_cur].clname;
                 addInName(e_cur,classr1);

		
	      }
	    }
            else Rlist[Cur_rep2].clname=e_num;
          
	  
	    classrn=num1;
	    subs=1;
	    
	 }
       }
       res_temp= leng1 > diff1+leng2  ? leng1 : diff1 + leng2 ;
       if (classr1 + res_temp <classr2 ){
	 
	 res_leng=res_temp;
       }
       else{
	 if (classr2<classr1){
           res_leng=res_temp;
	 }
         else res_leng=leng1;
       }
       cur_cr1=cr1;
       cur_cr2=cr2;
       cur_leng=leng2;
       
       Cur_rep1=Binary_Search (Lo,Hi,cr1);
      
        if (diff1 ==0 && leng1!=leng2){
	 subs=subs+1;}
        else if (diff1 !=0)  subs=subs+1;
	cr1=fr1;
	cr2=fr2;
	leng2=res_leng;
	
	buff=num2;
	num2=num1;
	fsign=0;
	exact=0; 
	
    }
    
    
    
    else {
      if(fsign && fr1){
	
	/*fprintf(f2, "here1 %10d %10d %10d %10d\n", fr1,fr2,leng1,fsign); */
	/*fprintf(stderr, "here1 %10d %10d %10d %10d\n", fr1,fr2,leng1,fsign); */
	Cur_rep1=Binary_Search (Lo,Hi,fr1);
	Cur_rep2=Binary_Search (Lo,Hi,fr2);
	
	
	addInName(e_num,fr1);
	
	Elist[e_num].name= fr1; 
	Elist[e_num].length=leng1;
	Elist[e_num].copies=fsign;
	
	
	if (Rlist[Cur_rep1].clname==0){
	  Rlist[Cur_rep1].clname=e_num;
	  e_cur=e_num;
                
	}
	else {           
	  ci2=0; csign2=0;
	  cur_e=Rlist[Cur_rep1].clname;
	  new_cur=Binary_Search(Lo,Hi,Elist[cur_e].name) ;
	  
	  e_cur=Rlist[new_cur].clname;
	  
	  addInName(e_cur,fr1);
	  
	  
	}
	
          
	if (Rlist[Cur_rep2].clname==0) Rlist[Cur_rep2].clname=e_num;
	else {
	  
	  ci2=0; csign2=0;
	  cur_e=Rlist[Cur_rep2].clname;
	  new_cur=Binary_Search(Lo,Hi,Elist[cur_e].name) ;
	  
	  e_cur=Rlist[new_cur].clname;
	  
	  addInName(e_cur,fr1);
	  
	  
	}
	e_num++;
	
      }
      
      if (res_leng){
	if (cr2>classr1 && classr1 + res_leng> cr1 && cr2<cr1 && cr1 -classr1 +1>=Igap){
	  res_leng=cr1 -classr1 +1;
	}
	/*
	  fprintf(f2, "here3 %10d %10d %10d %10d\n",classr1,classr2, res_leng,subs);
	*/
	/*  fprintf(stderr, "here3 %10d %10d %10d %10d\n",classr1,classr2, res_leng,subs);*/
	
	
	Elist[e_num].length=res_leng;
	Elist[e_num].copies=subs;
	
	if (Rlist[Cur_rep1].clname==0){
	  Rlist[Cur_rep1].clname=e_num;
	  e_cur=e_num;
	  
	}
	else{
	  if (Rlist[Cur_rep1].clname!=e_num){
	    ci2=0; csign2=0;
	    cur_e=Rlist[Cur_rep1].clname;
	    new_cur=Binary_Search(Lo,Hi,Elist[cur_e].name) ;
	    e_cur=Rlist[new_cur].clname;
	    addInName(e_cur,classr1);
	    
	  }
	}
	
	
	ci2=0; csign2=0;
	new_cur=Binary_Search(Lo,Hi,cur_cr2) ;
	
	e_cur=Rlist[new_cur].clname;
	if (e_cur !=0){
	  addInName(e_cur,classr1);
	  
	  
	} else Rlist[new_cur].clname=e_num; 
	e_num++;
	res_leng=0;
	subs=0;
      }
      fsign =1;   
      
    }
     
     if (res_leng){
       if (Rlist[Cur_rep1].clname==0){
	 Rlist[Cur_rep1].clname=e_num;
	 e_cur=e_num;
	 
       }
       else{
	 
	 ci2=0; csign2=0;
	 cur_e=Rlist[Cur_rep1].clname;
	 new_cur=Binary_Search(Lo,Hi,Elist[cur_e].name) ;
	 
	 
	 e_cur=Rlist[new_cur].clname;
	 addInName(e_cur,classr1); 
	 
       }
       
       
  
       ci2=0; csign2=0;
       new_cur=Binary_Search(Lo,Hi,cur_cr2) ;
       
       
       e_cur=Rlist[new_cur].clname;
       if (e_cur!=0){
	 addInName(e_cur,classr1);
	      
       }
       else Rlist[new_cur].clname=e_num;
       
       
       /*fprintf(f2, "here5 %10d %10d %10d %10d \n",classr1,cur_cr2,zero,zero);*/
       /*fprintf(stderr, "here5 %10d %10d %10d %10d \n",classr1,cur_cr2,zero,zero);*/
       
     }
     
     fr1=cr1;
     fr2=cr2;
     num1=num2;
     leng1=leng2;
     
 }
 
 if(fsign && fr1){
   
    /*fprintf(f2, "here6 %10d %10d %10d %10d\n", fr1,fr2,leng1,fsign); */
   /*fprintf(stderr, "here6 %10d %10d %10d %10d\n", fr1,fr2,leng1,fsign);*/
   Cur_rep1=Binary_Search (Lo,Hi,fr1);
   Cur_rep2=Binary_Search (Lo,Hi,fr2);
   
   Elist[e_num].name= fr1; 
   Elist[e_num].length=leng1;
   Elist[e_num].copies=fsign;
   addInName(e_num,fr1);
   
   if (Rlist[Cur_rep1].clname==0){
     Rlist[Cur_rep1].clname=e_num;
     e_cur=e_num;
     
   }
   else {           
     ci2=0; csign2=0;
     cur_e=Rlist[Cur_rep1].clname;
     new_cur=Binary_Search(Lo,Hi,Elist[cur_e].name) ;
     
     e_cur=Rlist[new_cur].clname;
     addInName(e_cur,fr1);
     
     
     
   }
   if (Rlist[Cur_rep2].clname==0) Rlist[Cur_rep2].clname=e_num;
   else {
     
     ci2=0; csign2=0;
     cur_e=Rlist[Cur_rep2].clname;
     new_cur=Binary_Search(Lo,Hi,Elist[cur_e].name) ;
     
     e_cur=Rlist[new_cur].clname;
     addInName(e_cur,fr1);
     
     
   }
   e_num++;       
   
 } 
 if (res_leng){
   if (cr2>classr1 && classr1 + res_leng> cr1 && cr2<cr1 && cr1 -classr1 +1>=Igap){
     res_leng=cr1 -classr1 +1;
   }
   /*fprintf(f2, "here7 %10d %10d %10d %10d\n",classr1,classr2, res_leng,subs);*/
   /*fprintf(stderr, "here7 %10d %10d %10d %10d\n",classr1,classr2, res_leng,subs);*/
   
   
   Elist[e_num].length=res_leng;
   Elist[e_num].copies=subs;
	 
   if (Rlist[Cur_rep1].clname==0){
     Rlist[Cur_rep1].clname=e_num;
     e_cur=e_num;
     
   }
   else{
     if (Rlist[Cur_rep1].clname!=e_num){
       ci2=0; csign2=0;
       cur_e=Rlist[Cur_rep1].clname;
       new_cur=Binary_Search(Lo,Hi,Elist[cur_e].name) ;
       
       e_cur=Rlist[new_cur].clname;
       addInName(e_cur,classr1 );
       
     }
   }
   
   
   ci2=0; csign2=0;
   new_cur=Binary_Search(Lo,Hi,cur_cr2) ;
   
   
   e_cur=Rlist[new_cur].clname;
   if (e_cur !=0){
     addInName(e_cur,classr1); 
     
     
   } else Rlist[new_cur].clname=e_num; 
   
   e_num++;
   
   res_leng=0;
   subs=0;
 }
 
 
 
 fclose(fp);
 
 
 /*printf ("The end of preprocessing\n"); */  
 /*fprintf (stderr, " We are before Rep_map\n");*/
 Rep_map(Data_Len, Data_Len3);
 /*fprintf (stderr, " We are after Rep_map\n");*/
 
 classPrint(Data_Len3);
 /*RepeaPrint(Data_Len3,argv [argc - 1]);*/
 RepeaPrint(Data_Len3);
 free(Rlist);
 free(Elist);
 free(Clist);
 
}

  int  Binary_Esearch2 ( int Lo, int Hi, int Val)

/* Find and return  i  such that  Lo <= i <= Hi  and
*  A [i].name <= Val < A [i + 1].name .  If no such value return  -1 . */

  {
   int  i;
   
   /*printf ("%10d %10d %10d\n ", Lo, Hi,i);*/
   while  (Lo <= Hi)
     {
       
      i = (Lo + Hi) / 2;
      /*printf ("%10d %10d %10d %10d\n ", Lo, Hi,i, Rlist[i].name);*/
      if  (Elist[i].name == Val)
          return  i;
      else if  (Val > Elist [i].name)
          Hi = i - 1;
        else
          Lo = i + 1;
     }

   return  -1;
  }
  int  Binary_Esearch ( int Lo, int Hi, int Val)

/* Find and return  i  such that  Lo <= i <= Hi  and
*  A [i].name <= Val < A [i + 1].name .  If no such value return  -1 . */

  {
   int  i;
   
   /*printf ("%10d %10d %10d\n ", Lo, Hi,i);*/
   while  (Lo <= Hi)
     {
       
      i = (Lo + Hi) / 2;
      /*printf ("%10d %10d %10d %10d\n ", Lo, Hi,i, Rlist[i].name);*/
      if  (Elist[i].name == Val)
          return  i;
      else if  (Val < Elist [i].name)
          Hi = i - 1;
        else
          Lo = i + 1;
     }

   return  -1;
  }
  int  Binary_Search ( int Lo, int Hi, int Val)

/* Find and return  i  such that  Lo <= i <= Hi  and
*  A [i].name <= Val < A [i + 1].name .  If no such value return  -1 . */

  {
   int  i;
   
   /*printf ("%10d %10d %10d\n ", Lo, Hi,i);*/
   while  (Lo <= Hi)
     {
       
      i = (Lo + Hi) / 2;
      /*printf ("%10d %10d %10d %10d\n ", Lo, Hi,i, Rlist[i].name);*/
      if  (Rlist[i].name == Val)
          return  i;
      else if  (Val < Rlist [i].name)
          Hi = i - 1;
        else
          Lo = i + 1;
     }

   return  -1;
  }


inPool_Ptr inPalloc(void)
{
  return (inPool_Ptr)malloc(sizeof(inPool_Type));
}

inName_Ptr inNalloc(void)
{
  return (inName_Ptr)malloc(sizeof(inName_Type));
}

int addInName (int e_num ,int clNname)
{
  inName_Ptr Temp, Tname, Nname,Prev;
  int stops;

  if ( Elist[e_num].clname==NULL )
    {
      Nname=inNalloc();
      Nname->name= clNname;
      Nname->next=NULL;
      Nname->prev=NULL;
      Elist[e_num].clname=Nname;  
      Prev=Nname;
      
    }
  else {
    stops=1;
     for (Tname=Elist[e_num].clname;(Tname!=NULL && stops); Tname=Tname->next){     
      Prev=Tname->prev;
      if (Tname->name==clNname) stops=0;
      if (Tname->name>clNname) {
	Nname=inNalloc();
	
	if(Prev!=NULL){
	  Prev->next=Tname;
	 
	  Nname->name=clNname;
	  Nname->prev=Tname->prev;
	  Nname->next=Tname;
	  Tname->prev=Nname;
	}
	else{
          
	  Nname->prev=NULL;
	  Nname->name=clNname;
	  Nname->next=Tname;
	  Tname->prev=Nname;
          Elist[e_num].clname=Nname;  
	}
	stops=0;
      }
        Prev=Tname;
     }
     if (Tname==NULL && stops) {
   
       Nname=inNalloc();
       if(Prev!=NULL){
	 Prev->next=Nname;
         /*printf ( "%10d %10d %10d\n",e_num,Prev->next->name, clNname); */
       }
       Nname->name=clNname;
       Nname->prev=Prev; 
       Nname->next=NULL; 
     }
  }
}

int addInPool(int i1,int i,int es)
{
  inPool_Ptr Temp,Tpool,T1pool, Prev, Npool;
  
  inName_Ptr Tname;
  int j1,stops;
  if(Clist[i1].pool==NULL)
    {
   
      Npool=inPalloc();
      Npool->name=Elist[i].clname->name;
      Npool->prev=NULL;
      Npool->next=NULL;
      Clist[i1].pool=Npool;
    
      Prev= Npool;
      
      for (Tname=Elist[i].clname->next;(Tname!=NULL); Tname=Tname->next){
              
	Npool=inPalloc();
        if(Prev!=NULL)
        Prev->next=Npool;
	Npool->name=Tname->name;
	Npool->next=NULL;
	Npool->prev=Prev;
	Prev=Npool;
	
      }
     
    }
  else 
    if (es)
      {
	for (Tname=Elist[i].clname;(Tname!=NULL); Tname=Tname->next){
        

	
        stops=1;
         
	  for (Tpool=Clist[i1].pool;(Tpool!=NULL && stops); Tpool=Tpool->next){
   
            Prev=Tpool->prev;
	    if (Tpool->name==Tname->name) stops=0;
            if (Tpool->name>Tname->name) {
	      Npool=inPalloc();
              if(Prev!=NULL){
              Prev->next=Npool;

	      Npool->name=Tname->name;
	      Npool->prev=Tpool->prev;
	      Npool->next=Tpool;
	      Tpool->prev=Npool;
	      }
              else{
		Npool->prev=NULL;
                Npool->name=Tname->name;
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
	    Npool->name=Tname->name;
	    Npool->prev=Prev; 
            Npool->next=NULL; 
	  }
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
  /*  fprintf (stderr," %10d ",i1);
  for (Tpool=Clist[i1].pool;(Tpool!=NULL); Tpool=Tpool->next){
  fprintf (stderr," %10d ",Tpool->name);
  } 
  fprintf (stderr," \n ");*/
}

void Rep_map(int Data_Len2,int Data_Len3)
{ 
  int i,j,ib,pleng,i1,i2,j1,cb,cb1,es,stopsgn,bsign;
  int stops;
  inPool_Ptr Tpool;
  inName_Ptr Tname;
  i=1;j=0;
  pleng=0;
  cb=0;
  /*fprintf (stderr, " We are in Rep_map\n");!!!!*/
  while(i<Data_Len2 && pleng<Data_Len3){
    /*printf("Ei =%10d %10d\n",i,Elist[i].name);*/
    
    if (Elist[i].length!=0){
      cb=0;for(ib=0;Blist[ib].clcoor!=0;ib++)Blist[ib].clcoor=0;
     
      for (i1=0;(Clist[i1].name !=0 && i1<Data_Len3 );i1++){
	
	if (Clist[i1].name!=-1){
	  stopsgn=1;
           for (Tname=Elist[i].clname;(Tname!=NULL && stopsgn); Tname=Tname->next){ 
	  
            stops=1;
	    for (Tpool=Clist[i1].pool;(Tpool!=NULL && stops); Tpool=Tpool->next){
	      if (Tpool->name==Tname->name) stops=0;
            
	    }
	    if (!stops)
	      {     
		bsign=0;stopsgn=0;
		for(ib=0;Blist[ib].clcoor!=0;ib++)
		  {
		    if(Blist[ib].clcoor==i1+1) bsign=1;
		  }
		if(!bsign){
		  Blist[cb].clcoor=i1+1;
		  cb++;
		}
	      }
	  }
	  
	}
      }
      /*printf("here2 cb= %10d \n",cb);*/
      /*printf ("Blist :");
      for(ib=0;Blist[ib].clcoor!=0;ib++) printf (" %10d ",Blist[ib].clcoor-1);
      printf( "\n");*/
      if(cb==0){
	Clist[pleng].name=Elist[i].name;
        es=1;
	addInPool(pleng,i,es);
      
	Clist[pleng].num=j1;
	pleng++;
	/*printf("\n pleng %10d \n", pleng);*/
      }
      if(cb==1){
	es=1;
	i1=Blist[0].clcoor-1;
	addInPool(i1,i,es);
	
      }
      if (cb >1){
		i1=Blist[0].clcoor-1;
	es=0;
	for(cb1=1;cb1<cb;cb1++){
	  i2=Blist[cb1].clcoor-1;
	  addInPool(i1,i2,es);
        
	  /* printf("Cnames :%10d %10d \n", Clist[i2].name, Clist[i1].name);*/
	  Clist[i2].name=-1;
	}
	es=1;
	addInPool(i1,i,es);
        
      }
    }
    i++;
  }
}
void classPrint(int Data_Len3)
{
  FILE  *f4, *f7;
  int i,j,ib,pleng,i1,i2,j1,cb,cb1,es,stopsgn,bsign,ci;
  int Lo, Hi;
  int stops,new_cur;
  inPool_Ptr T1pool;
  /*f4 = fopen ("Repeat_classes","a");
  f7= fopen ("Repeat_classes_min","a");
  fprintf (f4,"    Class    MaxLength   Copies \n");*/  
  i=1;
  Lo=1;
  while (Elist[i].name>0)
    i++;
  Hi=i-1;
  i=0;j=0;
  while (Clist[i].name!=0 && i<Data_Len3)
    {
      if(Clist[i].name !=-1){
	j++;
        Clist[i].length=0;
        Clist[i].copies=0;
        Clist[i].min_name= Clist[i].name;
        Clist[i].min_length=10000;
       
	 for(T1pool=Clist[i].pool;(T1pool!=NULL); T1pool=T1pool->next){
	  new_cur=Binary_Esearch(Lo,Hi,T1pool->name) ;
	   T1pool->length=Elist[new_cur].length;
	   T1pool->copies=Elist[new_cur].copies;
	   
	  if (Elist[new_cur].length <Clist[i].min_length){
	    Clist[i].min_length=Elist[new_cur].length;
	    Clist[i].min_name=Elist[new_cur].name;
	    
	  }
	  
	  if (Elist[new_cur].length >Clist[i].length && Elist[new_cur].name>0){
	    Clist[i].length=Elist[new_cur].length;
	    Clist[i].name=Elist[new_cur].name;
	    
	  }
	  
	  
	  Clist[i].copies=Clist[i].copies+Elist[new_cur].copies;
	  Elist[new_cur].length=-1;
	}
            
	 /*fprintf (f4," %10d %10d %10d %10d %10d \n", Clist[i].name, Clist[i].length, Clist[i].copies, Clist[i].min_name, Clist[i].min_length);
	   fprintf (f7," %10d %10d %10d \n", Clist[i].name, Clist[i].min_name, Clist[i].min_length);*/
      }
      i++;
    }
  /* printf ("The number of the repeat classes  %10d \n",j);
  fclose(f4);
  fclose(f7);*/
}


/****Preparing for multyfasta file******/

/*void RepeaPrint(int Data_Len3,const char * Filename)*/
void RepeaPrint(int Data_Len3)
{
  FILE *f3,*f8 ,*f9;
  int i,j,k,z,j4,j1,ci,ci1;
 inPool_Ptr T1pool,Tpool,Prev,Npool;
 char c[100];
 /*f3 =  File_Open(Filename,"a") ;*/
 /*f9 = fopen("gnu.com6","a");*/
 i=0;j=0;
  while (Clist[i].name!=0 && i<Data_Len3)
    {
      if(Clist[i].name !=-1){
	j++;
		
	  ci=0;
	 for(Tpool=Clist[i].pool;(Tpool!=NULL); Tpool=Tpool->next){
	   ci++;
           Prev=Tpool;
	   
           /*fprintf(f3, ">Class %10d ",i);
	     fprintf(f3, "%10d %10d %10d \n",Tpool->name,Tpool->length,Tpool->copies); */
	   printf(">Class %10d ",i);
	   printf("%10d %10d %10d \n",Tpool->name,Tpool->length,Tpool->copies);
	 }
	 if (ci==1) {
	   /* fprintf(f3, ">Class %10d ",i);
	  fprintf(f3, "%10d %10d %10d \n",Prev->name+Clist[i].length,Clist[i].length,Clist[i].copies);
	   */
          printf( ">Class %10d ",i);
	  printf( "%10d %10d %10d \n",Prev->name+Clist[i].length,Clist[i].length,Clist[i].copies);
	    Npool=inPalloc();
          if(Prev!=NULL)
           Prev->next=Npool;
	  Npool->name=Prev->name+Clist[i].length;
	  Npool->length=Clist[i].length;
	  Npool->copies=1;
	  Npool->prev=Prev; 
	  Npool->next=NULL; 
	}
	
      }
      
      i++;
      }  
  /*fprintf(f3, ">The end. ");*/
  /*printf ("The number of the repeat classes  %10d \n",j); */
  /*fclose(f3);     */
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





















 
























