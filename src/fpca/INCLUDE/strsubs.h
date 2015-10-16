#include <stdlib.h>

int splitup (char *strin,  char *strpt[],int maxpt) ;
int oldsplitup (char *strin,  char *strpt[],int maxpt) ;
void freeup (char *strpt[],int numpt) ;
int split1 (char *strin,  char *strpt[], char splitc);
int first_word(char *string, char *word, char *rest) ;
char * fnwhite (char * ss) ;
char * fwhite (char * ss) ;
int isnumber (char c) ;
int isnumword (char *str)  ;
void fatalx( char *fmt, ...) ;
long seednum() ;
void printbl(int n) ;
void striptrail(char *sss, char c) ;
void catx(char *sout, char **spt, int n) ;
void catxx(char *sout, char **spt, int n) ;
void makedfn(char *dirname, char *fname, 
  char *outname, int maxstr) ;
int substring (char **ap, char *inx, char *outx) ;


#define ZALLOC(item,n,type)      if ((item = (type *)calloc((n),sizeof(type))) == NULL) \
                                        fatalx("Unable to allocate %d unit(s) for item \n",n)

#undef MAX
#undef MIN

#define MAX(a,b)   ( (a) < (b) ?  (b) : (a) ) 
#define MIN(a,b)   ( (a) < (b) ?  (a) : (b) ) 
#define YES  1
#define NO  0

