#ifndef WAVEPROP_CPP
#define WAVEPROP_CPP

#include "waveprop.h"

#include <fstream>
#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>

//#define FFTW_PLAN_MODE FFTW_PATIENT
#define N_FFT_THREADS 1
#include <mathlib/math/std_math.h>
#include <mathlib/math/PDE/laplacian_hdaf.h>
#include <mathlib/link.cpp>


#include <mathlib/tools/std_tools.h>

#include "pngwriter.h"



#define wp_pi 3.14159265358979323846264338327950

using namespace std;

 #ifdef __cplusplus
 extern "C" {
 #endif

#define _debug_here( pos ) ( printf( "debug: file %s; \n line %d; \t code %d\n", __FILE__ , __LINE__, pos ) );


// ordering: 0 if row major, 1 if column major

int render_png_scalar(
    char *fname, int nx, int ny, double *data, int ordering,
    double center, double major_scale,
    double red_midpoint, double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_midpoint, double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_midpoint, double blue_leftvalue, double blue_midvalue, double blue_rightvalue )
{
    //cout << "writing " << fname << endl;
    
	pngwriter png(nx, ny, 1.0, fname);
    
 	for (int i = 0; i<nx;i++)
 	for (int j = 0; j<ny;j++)
    {
        double r=0;
        if (ordering==0)
            r = (data[i*ny+j]-center)/major_scale;
        else
            r = (data[j*nx+i]-center)/major_scale;
        
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


int render_png_blend_background(
    char *fname, int nx, int ny, double *data, double *background, int ordering,
    double center, double major_scale, double midpoint
    double red_leftvalue, double red_rightvalue,
    double green_leftvalue, double green_rightvalue,
    double blue_leftvalue,  double blue_rightvalue )
{
	pngwriter png(nx, ny, 1.0, fname);
    
    if (midpoint <= 1E-10) return 10;
    
    for (int i = 0; i<nx;i++)
 	for (int j = 0; j<ny;j++)
    {
        double r=0;
        double backgroud_red=
        if (ordering==0)
            r = (data[i*ny+j]-center)/major_scale;
        else
            r = (data[j*nx+i]-center)/major_scale;
        
        
        double s = atan(r)/ml_pi+0.5;
        double red=0.0, blue=0.0, green=0.0;
        
        if ( s <= midpoint )
        {
            s /= midpoint;
            red += (1.0-s)*red_leftvalue;
            red += s*red_leftvalue;
        }
        
        
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


int method1_init( void **pdata, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_sigma_1, double hdaf_sigma_2 )
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
    
    double gamma1 = 2.0*dx1*sqrt(2*hdaf_order_1+1)/hdaf_sigma_1;
    double gamma2 = 2.0*dx2*sqrt(2*hdaf_order_2+1)/hdaf_sigma_2;
    
    ((laplacian_2d_hdaf *)(data->Del2))->init( n1, n2, n1*dx1, n2*dx2, hdaf_order_1, hdaf_order_2, gamma1, gamma2 );
    
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


int method2_init( void **pdata, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_sigma_1, double hdaf_sigma_2 )
{
    return method1_init( pdata, n1, n2, dx1, dx2, velocity, damping, expansion_order, hdaf_order_1, hdaf_order_2, hdaf_sigma_1, hdaf_sigma_2 );
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

int resample2d_linear( double *in, double *out, int in_n1, int in_n2, int out_n1, int out_n2, double in_dx1, double in_dx2, double out_dx1, double out_dx2, double in_off1, double in_off2, double out_off1, double out_off2 )
{
    if ( in==NULL)  return 1;
    if (out==NULL)  return 2;
    
    for (int k=0; k<out_n1*out_n2; k++)
        out[k] = 0;
    
    for (int i=0; i<out_n1; i++)
    for (int j=0; j<out_n2; j++)
    {
        //cout << i << "\t" << j << endl;
        
        double x=out_off1+out_dx1*i;
        double y=out_off2+out_dx2*j;
        
        int u= (int) floor( (x-in_off1)/in_dx1 );
        int v= (int) floor( (y-in_off2)/in_dx2 );
        
        if ( 0 <= u and u < in_n1 and 0 <= v and v < in_n2 )
        {
            double x1 = in_off1+in_dx1*u;
            double x2 = in_off1+in_dx1*(u+1);
            double y1 = in_off2+in_dx2*v;
            double y2 = in_off2+in_dx2*(v+1);
            
            double F11=in[u*in_n2+v];
            double F12=in[u*in_n2+v+1];
            double F21=in[(u+1)*in_n2+v];
            double F22=in[(u+1)*in_n2+v+1];
            
            out[i*out_n2+j] = ( F11*(x2-x)*(y2-y) + F21*(x-x1)*(y2-y) + F12*(x2-x)*(y-y1) + F22*(x-x1)*(y-y1) )/(in_dx1*in_dx2);
        }
        else
            out[i*out_n2+j] = 0;
    }
    
    return 0;
}



int cosh_damping( double *data, int nx, int ny, double ax, double ay, double dx, double dy, double x1, double x2, double y1, double y2, double xs, double ys, double amp )
{
    if (data==NULL)  return 1;
    
    for (int i=0; i<nx; i++)
    for (int j=0; j<ny; j++)
    {
        double x= dx*i+ax;
        double y= dy*j+ay;
        
        double q = cosh((y1-y)/ys);
        data[i*ny+j]  = ( 1.0/(q*q) );
        
        q = cosh((y2-y)/ys);
        data[i*ny+j] += ( 1.0/(q*q) );
        
        q = cosh((x1-x)/xs);
        data[i*ny+j] += ( 1.0/(q*q) );
        
        q = cosh((x2-x)/xs);
        data[i*ny+j] += ( 1.0/(q*q) );
        
        data[i*ny+j] *= amp;
    }
    
    return 0;
}







 #ifdef __cplusplus
 }
 #endif









#endif
