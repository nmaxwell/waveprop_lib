
#include <fstream>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#include "/workspace/waveprop_lib/waveprop.h"


int main()
{
    int n1 = 256;
    int n2 = 256;
    double dx1 = 1.0/n1;
    double dx2 = 1.0/n2;
    double tf = 1.0;
    double dt = 0.1;
    
    int hdaf_order = 8;
    double hdaf_gamma = 0.8;
    int exp_order = 7;
    
    double * vel = (double *)malloc( n1*n2*8 );
    
    for (int k=0; k<n1*n2; k++)
        vel[k] = 1.0;
    
    double * damp = (double *)malloc( n1*n2*8 );
    
    for (int k=0; k<n1*n2; k++)
        damp[k] = 1.0;
    
    double * u = (double *)malloc( n1*n2*8 );
    double * v = (double *)malloc( n1*n2*8 );
    
    for (int k=0; k<n1*n2; k++)
        u[k] = 1.0;
    for (int k=0; k<n1*n2; k++)
        v[k] = 1.0;
    
    void *data=0;
    method1_init( &data, n1, n2, dx1, dx2, vel, damp, exp_order, hdaf_order, hdaf_order, hdaf_gamma, hdaf_gamma );
    
    double t = 0.0;
    
    while (t <= tf)
    {
        method1_execute( data, dt, u,v,u,v );
        
        t += dt;
    }
}
