#ifndef WAVEPROP_H
#define WAVEPROP_H

/*

documentation:

    nicholas.maxwell@gmail.com for questions.


methods:
    
    method1:
        2D, acoustic, regular cartesian grid, non-homogenous velocity, hdaf for laplacian, taylor expansion of acoustic propagator for time, absorbing boundaries.
    
    
    
*/


 #ifdef __cplusplus
 extern "C" {
 #endif


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

int method1_init( method1_data *data, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_gamma_1, double hdaf_gamma_2 );

int method1_free( method1_data *data );

int method1_execute( method1_data *data, double t, double *ui, double *vi, double *uf, double *vf );








 #ifdef __cplusplus
 }
 #endif


#endif
















