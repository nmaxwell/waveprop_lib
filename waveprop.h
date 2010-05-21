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


// ordering: 0 if row major, 1 if column major

int render_png_scalar(
    char *fname, int nx, int ny, double *data, int ordering,
    double center, double major_scale,
    double red_midpoint, double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_midpoint, double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_midpoint, double blue_leftvalue, double blue_midvalue, double blue_rightvalue );

int render_png_blend_background(
    char *fname, int nx, int ny, double *data, double *background, int ordering,
    double center, double major_scale, double midpoint
    double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_leftvalue, double blue_midvalue, double blue_rightvalue );



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

int method1_init( void **data, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_sigma_1, double hdaf_sigma_2 );

int method1_free( void **data );

int method1_execute_nodamping( void *data, double t, double *ui, double *vi, double *uf, double *vf );

int method1_execute( void *data, double t, double *ui, double *vi, double *uf, double *vf );


int method2_init( void **data, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_sigma_1, double hdaf_sigma_2 );

int method2_free( void **data );

int method2_execute( void *data, double t, double *ui, double *vi, double *uf, double *vf, double *forcing_1, double *forcing_2 );


int resample2d_linear( double *in, double *out, int in_n1, int in_n2, int out_n1, int out_n2, double in_dx1, double in_dx2, double out_dx1, double out_dx2, double in_off1, double in_off2, double out_off1, double out_off2 );

int cosh_damping( double *data, int nx, int ny, double ax, double ay, double dx, double dy, double x1, double x2, double y1, double y2, double xs, double ys, double amp );



 #ifdef __cplusplus
 }
 #endif


#endif
















