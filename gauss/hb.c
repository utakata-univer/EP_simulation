#include "stdio.h"
#include "ranlxs.h"
#include "math.h"

const int PI = 3.141592;
static int seed;

int main(void)
{   
    int i,n;
    float x[2*n],y[2*n];

    rlxs_init(1,seed);
    
    ranlxs(x,2*n);

    double gauss(void)
    {
        for (i=0;i<n;i++)
	
	y[2*i]=sqrt(-2*log(1-x[2*i]))*cos(2*PI*x[2*i+1]);
        y[2*i+1]=sqrt(-2*log(1-x[2*i]))*sin(2*PI*x[2*i+1]);

	return y[i];
    }
    
    printf("gauss = %17.7e\n",gauss);

    return 0;
}


