#include "phi42.h"

void hopping(int hop[V][2*D] ){ 
    int x, y, Lk;
    int xk, k, dxk ;


    for (x=0; x < V ; x++){
        Lk = V;
        y  = x;

        for (k=D-1; k >= 0; k--){
		   
            Lk/=L;
	    xk =y/Lk;
    	    y  =y-xk*Lk;

            if (xk<L-1) dxk = Lk;
	    else        dxk = Lk*(1-L);
	    hop[x][k] = x + dxk;

            if (xk>0)   dxk = -Lk;
	    else        dxk = Lk*(L-1);
            hop[x][k+D] = x + dxk;

       	}
    }
}
