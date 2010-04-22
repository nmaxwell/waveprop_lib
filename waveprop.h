#ifndef WAVEPROP_H
#define WAVEPROP_H

/*

documentation:

    nicholas.maxwell@gmail.com for questions.


methods:
    
    method1:
        2D, acoustic, regular cartesian grid, non-homogenous velocity, hdaf for laplacian, taylor expansion of acoustic propagator for time, exponential damping.
    
    method2:
        method1 + forcing term with linear approx of integral
    
    
    
    
    
    
compile fftw as

./configure --enable-sse2 CPPFLAGS=-fPIC CC=g++ -enable-threads

compile arprec as

./configure CXXFLAGS="-fPIC -O3" CXX=g++

also pngwriter with fPIC


*/




 #ifdef __cplusplus
 extern "C" {
 #endif



int render_png_scalar(
    char *fname, int nx, int ny, double *data, double center, double major_scale,
    double red_midpoint, double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_midpoint, double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_midpoint, double blue_leftvalue, double blue_midvalue, double blue_rightvalue );


int render_png_scalar_resample(
    char *fname, double *data,
    double a1, double b1, double a2, double b2, int nx, int ny,
    double a1_new, double b1_new, double a2_new, double b2_new, int nx_new, int ny_new,
    double center, double major_scale,
    double red_midpoint, double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_midpoint, double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_midpoint, double blue_leftvalue, double blue_midvalue, double blue_rightvalue,
    double red_default, double green_default, double blue_default );

struct method1_data {
    
    int n1,n2;
    double dx1,dx2;
    
    double *velocity;
    double *damping;
    int expansion_order;
    
    double **U;
    double **V;
    void *Del2;
    
    method1_data():n1(0),n2(0),dx1(0),dx2(0),velocity(0),damping(0),expansion_order(0),U(0),V(0),Del2(0) {}
};

int method1_init( void **data, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_gamma_1, double hdaf_gamma_2 );

int method1_free( void **data );

int method1_execute_nodamping( void *data, double t, double *ui, double *vi, double *uf, double *vf );

int method1_execute( void *data, double t, double *ui, double *vi, double *uf, double *vf );


int method2_init( void **data, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_gamma_1, double hdaf_gamma_2 );

int method2_free( void **data );

int method2_execute( void *data, double t, double *ui, double *vi, double *uf, double *vf, double *forcing_1, double *forcing_2 );





 #ifdef __cplusplus
 }
 #endif


#endif
















