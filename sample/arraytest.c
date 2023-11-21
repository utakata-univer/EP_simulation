#include "stdio.h"
#include "ranlxd.h"

int main(void)
{
    int i;
    double a,b[2];	

    for (i=0;i<2;i++)
   {
      ranlxd(&a,1);
      b[i]=a;
   }
  printf("%e\n",b[1]); 

  return 0;
}  
