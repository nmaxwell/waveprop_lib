#ifndef WAVEPROP_CPP
#define WAVEPROP_CPP

#include "waveprop.h"

//#include <iostream>
#include <fstream>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#define N_FFT_THREADS  3

#include <mathlib/math/std_math.h>
#include <mathlib/math/PDE/laplacian_hdaf.h>
#include <mathlib/link.cpp>


using namespace std;

 #ifdef __cplusplus
 extern "C" {
 #endif

#define _debug_here( pos ) ( printf( "debug: file %s; \n line %d; \t code %d\n", __FILE__ , __LINE__, pos ) );



int method1_free( method1_data *data )
{
    if ((*data).velocity != NULL) {
        free( (*data).velocity );
        (*data).velocity = 0;
    }
    
    if ((*data).damping != NULL) {
        free( (*data).damping );
        (*data).damping = 0;
    }
    
    if ((*data).U != NULL && data->expansion_order>0 ) {
        ml_free( (*data).U, data->expansion_order+2 );
        (*data).U = 0;
    }
    
    if ((*data).V != NULL && data->expansion_order>0 ) {
        ml_free( (*data).V, data->expansion_order+1 );
        (*data).V = 0;
    }
    
    if ((*data).U != NULL ) {
        free( (*data).U );
    }
    
    if ((*data).V != NULL ) {
        free( (*data).V );
    }
    
    if ((*data).Del2 != NULL) {
        delete (laplacian_2d_hdaf *)((*data).Del2);
        (*data).Del2 = 0;
    }
    
    return 0;
}


int method1_init( method1_data *data, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_gamma_1, double hdaf_gamma_2 )
{
    (*data).n1 = n1;
    (*data).n2 = n2;
    (*data).dx1 = dx1;
    (*data).dx2 = dx2;
    (*data).expansion_order = expansion_order;
    
    int N = n1*n2;
    
    if (velocity == NULL) return 1;
    if (damping == NULL) return 1;
    
    method1_free( data );
    
    (*data).velocity = ml_alloc<double> ( N*8 );
    (*data).damping = ml_alloc<double> ( N*8 );
    (*data).U = ml_alloc<double> (expansion_order+2, N);
    (*data).V = ml_alloc<double> (expansion_order+1, N);
    
    if ((*data).velocity == NULL) return 2;   
    if ((*data).damping == NULL) return 2;
    if ((*data).U == NULL) return 2;
    if ((*data).V == NULL) return 2;
    
    for (int k=0; k<N; k++) {
        (*data).damping[k] = damping[k];
        (*data).velocity[k] = velocity[k];
    }
    
    (*data).Del2 = (void *)( new laplacian_2d_hdaf );
    ((laplacian_2d_hdaf *)(data->Del2))->init( n1, n2, n1*dx1, n2*dx2, hdaf_order_1, hdaf_order_2, hdaf_gamma_1, hdaf_gamma_2 );
    
    return 0;
}


int method1_execute( method1_data *data, double t, double *ui, double *vi, double *uf, double *vf )
{
    int N= (*data).n1*(*data).n2;
    double **U = (*data).U; 
    double **V = (*data).V;
    double *velocity = (*data).velocity;
    laplacian_2d_hdaf *Del2 = (laplacian_2d_hdaf *)(data->Del2);
    int exp_order = (*data).expansion_order;
    double * damping = (*data).damping;
    
    
    for (int k=0; k<N; k++)
        U[0][k] = ui[k];
    
    for (int k=0; k<N; k++)
        V[0][k] = vi[k];
    
    for (int n=1; n<=exp_order+1; n++)
    {
        Del2->execute( U[n-1], U[n] );
        
        for (int k=0; k<N; k++)
        {
            U[n][k] *= velocity[k];
            U[n][k] *= t*t/(4*n*n-2*n);
        }
    }
    
    for (int n = 1; n<=exp_order; n++)
    {
        Del2->execute( V[n-1], V[n] );
        
        for (int k=0; k<N; k++)
        {
            V[n][k] *= velocity[k];
            V[n][k] *= t*t/(4*n*n-2*n);
        }
    }
    
    for (int k=0; k<N; k++)
    {
        uf[k] = 0;
        vf[k] = 0;
    }
    
    for (int n = 0; n<=exp_order; n++)
    {
        double a_n = t/(2.0*n+1.0);
        double b_n = ((4.0*n+6.0)*n+2.0)/(t*(2.0*n+1.0));
        
        for (int k=0; k<N; k++)
        {
            uf[k] += U[n][k];
            uf[k] += V[n][k]*a_n;
            
            vf[k] += U[n+1][k]*b_n;
            vf[k] += V[n][k];
        }
    }
    
    for (int k=0; k<N; k++)
        uf[k] *= exp(-t*damping[k]);
    
    return 0;
}



 #ifdef __cplusplus
 }
 #endif
















#endif
