



from ctypes import *


c_waveprop = cdll.LoadLibrary('/usr/lib/libwaveprop.so.1')






#import numpy


#c_waveprop = cdll.LoadLibrary('/usr/lib/libhdaf.so.1')

"""
int method1_init( method1_data *data, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_gamma_1, double hdaf_gamma_2 );

int method1_free( method1_data *data );

int method1_execute( method1_data *data, double t, double *ui, double *vi, double *uf, double *vf );
"""


"""
class py_method1_data(Structure):
    _fields_ = [ (n1, c_int),(n2, c_int),(dx1, c_double),(dx2, c_double),  ]

    int n1,n2;
    double dx1,dx2;
    
    double *velocity;
    double *damping;
    int expansion_order;
    
    double **U;
    double **V;
    void *Del2;


c_waveprop.method1_init.argtyprs = [ c_double, c_int ]


"""













"""









def freq_to_sigma( f, m ):
    return c_hdaf.sigma_from_cutoff_frequency( c_double(f), c_int(m) )



c_hdaf.hdaf_equate_arrays.argtypes = [ c_void_p, c_void_p, c_ulong ]

c_hdaf.hdaf_free_array.argtypes = [ c_void_p ]




c_hdaf.get_hdaf_kernel.argtypes=[ c_void_p, c_void_p, c_double, c_int, c_double, c_char_p  ]
c_hdaf.get_hdaf_kernel.restype = c_int


def hdaf_kernel( sampling_period, order, sigma ):
    
    ptr = c_void_p()
    size = c_ulong()

    error = c_hdaf.get_hdaf_kernel( byref(ptr), byref(size), c_double(sampling_period), c_int(order), c_double(sigma), c_char_p(hdaf_data_dir) )
    
    if error != 0:
        print "Error in 'get_hdaf_kernel'; code ", error
    
    
    ar = numpy.zeros(size.value)
    c_hdaf.hdaf_equate_arrays( ar.ctypes.data_as(c_void_p), ptr, size )
    c_hdaf.hdaf_free_array( ptr )
    ptr = c_void_p()
    
    return ar




c_hdaf.get_hdaf_kernel_arbitrary_points.argtypes=[ c_void_p, c_void_p, c_int, c_int, c_double, c_char_p  ]
c_hdaf.get_hdaf_kernel_arbitrary_points.restype = c_int


def hdaf_kernel_bypts( eval_points, order, sigma ):
    
    eval_points = numpy.array(eval_points)
    ar = numpy.zeros(len(eval_points))
    length = c_int(len(eval_points))
    
    error = c_hdaf.get_hdaf_kernel_arbitrary_points( eval_points.ctypes.data_as(c_void_p), ar.ctypes.data_as(c_void_p), length, c_int(order), c_double(sigma), c_char_p(hdaf_data_dir) )
    
    if error != 0:
        print "Error in 'get_hdaf_kernel_arbitrary_points'; code ", error
    
    return ar






c_hdaf.get_hdaf_kernel_lp.argtypes=[ c_void_p, c_void_p, c_double, c_int, c_double, c_char_p  ]
c_hdaf.get_hdaf_kernel_lp.restype = c_int


def lp_hdaf_kernel( sampling_period, order, cutoff_frequency ):
    
    ptr = c_void_p()
    size = c_ulong()

    error = c_hdaf.get_hdaf_kernel_lp( byref(ptr), byref(size), c_double(sampling_period), c_int(order), c_double(cutoff_frequency), c_char_p(hdaf_data_dir) )
    
    if error != 0:
        print "Error in 'get_hdaf_kernel_lp'; code ", error
    
    
    ar = numpy.zeros(size.value)
    c_hdaf.hdaf_equate_arrays( ar.ctypes.data_as(c_void_p), ptr, size )
    c_hdaf.hdaf_free_array( ptr )
    ptr = c_void_p()
    
    return ar




"""





"""

int get_hdaf_kernel(double *&kernel, int & kernel_size, double sampling_period, int order, double sigma, const char *hdaf_data_dir );

int get_hdaf_kernel_lp(double *&kernel, int & kernel_size, double sampling_period, int order, double cutoff_frequency, const char *hdaf_data_dir );

int get_hdaf_kernel_bp(double *&kernel, int & kernel_size, double sampling_period, int low_pass_order, double low_pass_frequency, int high_pass_order, double high_pass_frequency, const char *hdaf_data_dir );
"""



        





