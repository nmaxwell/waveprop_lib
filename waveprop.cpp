#ifndef WAVEPROP_CPP
#define WAVEPROP_CPP

#include "waveprop.h"

#include <fstream>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

#define N_FFT_THREADS 1
#include <mathlib/math/std_math.h>
#include <mathlib/math/PDE/laplacian_hdaf.h>
#include <mathlib/link.cpp>


#include <mathlib/tools/std_tools.h>

#include "pngwriter.h"





using namespace std;

 #ifdef __cplusplus
 extern "C" {
 #endif

#define _debug_here( pos ) ( printf( "debug: file %s; \n line %d; \t code %d\n", __FILE__ , __LINE__, pos ) );




int render_png_scalar(
    char *fname, int nx, int ny, double *data, double center, double major_scale,
    double red_midpoint, double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_midpoint, double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_midpoint, double blue_leftvalue, double blue_midvalue, double blue_rightvalue )
{
    //cout << "writing " << fname << endl;
    
	pngwriter png(nx, ny, 1.0, fname);
    
 	for (int i = 0; i<nx;i++)
 	for (int j = 0; j<ny;j++)
    {
        double r = (data[i*ny+j]-center)/major_scale;
        double s = atan(r)/ml_pi+0.5;
        double red,green,blue;
        
        if ( s <= red_midpoint )
            red = (s/red_midpoint)*(red_midvalue-red_leftvalue)+red_leftvalue;
        else
            red = (s-red_midpoint)*(red_rightvalue-red_midvalue)/(1.0-red_midpoint)+red_midvalue;
        if ( s <= green_midpoint )
            green = (s/green_midpoint)*(green_midvalue-green_leftvalue)+green_leftvalue;
        else
            green = (s-green_midpoint)*(green_rightvalue-green_midvalue)/(1.0-green_midpoint)+green_midvalue;
        if ( s <= blue_midpoint )
            blue = (s/blue_midpoint)*(blue_midvalue-blue_leftvalue)+blue_leftvalue;
        else
            blue = (s-blue_midpoint)*(blue_rightvalue-blue_midvalue)/(1.0-blue_midpoint)+blue_midvalue;
        
 		png.plot(i+1,j+1, red, green, blue );
	}
    
	png.close();
    
    return 0;
}



int method1_free( void **pdata )
{
    method1_data *&data = *((method1_data **)pdata);
    
    if ( data != NULL)
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
        
        free(data);
        data = 0;
    }
    
    return 0;
}


int method1_init( void **pdata, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_gamma_1, double hdaf_gamma_2 )
{
    method1_free( pdata );
    
    if (*((method1_data **)pdata)!=NULL)
        return 9;
    
    *((method1_data **)pdata) = (method1_data *)malloc(sizeof(method1_data));
    
    if (*((method1_data **)pdata)==NULL)
        return 15;
    
    method1_data *data = *((method1_data **)pdata);  
    
    (*data).n1 = n1;
    (*data).n2 = n2;
    (*data).dx1 = dx1;
    (*data).dx2 = dx2;
    (*data).expansion_order = expansion_order;
    
    int N = n1*n2;
    
    if (velocity == NULL) return 1;
    if (damping == NULL) return 1;
    
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


int method1_execute_nodamping( void *pdata, double t, double *ui, double *vi, double *uf, double *vf )
{
    // split dapming and rest into two steps, usefull for other methods to reuse code, i.e. forcing termp happens before damping.
    
    method1_data *data = (method1_data *)pdata;
    
    if (data==NULL) return 10;
    
    int N= (*data).n1*(*data).n2;
    double **U = (*data).U; 
    double **V = (*data).V;
    double *velocity = (*data).velocity;
    laplacian_2d_hdaf *Del2 = (laplacian_2d_hdaf *)(data->Del2);
    int exp_order = (*data).expansion_order;
    
    for (int k=0; k<N; k++)
        U[0][k] = ui[k];
    
    for (int k=0; k<N; k++)
        V[0][k] = vi[k];
    
    for (int n=1; n<=exp_order+1; n++)
    {
        Del2->execute( U[n-1], U[n] );
        double c = t*t/((4*n-2)*n);
        
        for (int k=0; k<N; k++)
        {
            U[n][k] *= velocity[k];
            U[n][k] *= c;
        }
    }
    
    for (int n = 1; n<=exp_order; n++)
    {
        Del2->execute( V[n-1], V[n] );
        double c= t*t/((4*n-2)*n);
        
        for (int k=0; k<N; k++)
        {
            V[n][k] *= velocity[k];
            V[n][k] *= c;
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
    
    
    return 0;
}


int method1_execute( void *pdata, double t, double *ui, double *vi, double *uf, double *vf )
{
    int err = method1_execute_nodamping( pdata, t, ui, vi, uf, vf );
    
    if (err == 0)
    {
        method1_data *data = (method1_data *)pdata;
        if (data==NULL) return 10;
        int N= (*data).n1*(*data).n2;
        double * damping = (*data).damping;
        
        for (int k=0; k<N; k++)
            uf[k] *= exp(-t*damping[k]);
    }
    
    return err;
}


int method2_free( void **pdata )
{
    return method1_free( pdata );
}


int method2_init( void **pdata, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_gamma_1, double hdaf_gamma_2 )
{
    return method1_init( pdata, n1, n2, dx1, dx2, velocity, damping, expansion_order, hdaf_order_1, hdaf_order_2, hdaf_gamma_1, hdaf_gamma_2 );
}


int method2_execute( void *pdata, double t, double *ui, double *vi, double *uf, double *vf, double *forcing_1, double *forcing_2 )
{
    int err = method1_execute_nodamping( pdata, t, ui, vi, uf, vf );
    
    if (err == 0)
    {
        method1_data *data = (method1_data *)pdata;
        if (data==NULL) return 10;
        int N= (*data).n1*(*data).n2;
        double **U = (*data).U;
        double *velocity = (*data).velocity;
        double * damping = (*data).damping;
        laplacian_2d_hdaf *Del2 = (laplacian_2d_hdaf *)(data->Del2);
        int exp_order = (*data).expansion_order;
        
        for (int k=0; k<N; k++)
            U[0][k] = forcing_1[k];
        
        for (int n=1; n<=exp_order+1; n++)
        {
            Del2->execute( U[n-1], U[n] );
            double c = t*t/((4*n-2)*n);
            
            for (int k=0; k<N; k++)
            {
                U[n][k] *= velocity[k];
                U[n][k] *= c;
            }
        }
        
        for (int n = 0; n<=exp_order; n++)
        {
            double a_n = t/2.0;
            double b_n = ((4.0*n+6.0)*n+2.0)/(2*(2*n+1));
            
            for (int k=0; k<N; k++)
            {
                uf[k] += U[n][k]*a_n;
                vf[k] += U[n+1][k]*b_n;
            }
        }
        
        for (int k=0; k<N; k++)
            uf[k] += (t/2.0)*forcing_2[k];
        
        // apply damping after the forcing term, hence the use of method1_execute_nodamping
        for (int k=0; k<N; k++)
            uf[k] *= exp(-t*damping[k]);
    }
    
    return err;
}




















 #ifdef __cplusplus
 }
 #endif









#endif
