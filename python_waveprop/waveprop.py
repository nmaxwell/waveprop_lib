

from ctypes import *
from ctypes.util import *
import numpy

"""
int method1_init( method1_data *data, int n1, int n2, double dx1, double dx2, double *velocity, double *damping, int expansion_order, int hdaf_order_1, int hdaf_order_2, double hdaf_gamma_1, double hdaf_gamma_2 );

int method1_free( method1_data *data );

int method1_execute( method1_data *data, double t, double *ui, double *vi, double *uf, double *vf );
"""

c_waveprop = cdll.LoadLibrary(find_library('waveprop'))

c_waveprop.method1_init.argtyprs = [ c_void_p, c_int, c_int, c_double, c_double, c_void_p, c_void_p, c_int, c_int, c_int, c_double, c_double ]

c_waveprop.method1_free.argtyprs = [ c_void_p ] 

c_waveprop.method1_execute.argtypes = [ c_void_p, c_double, c_void_p, c_void_p,  c_void_p,  c_void_p, ]

class py_method1_data(Structure):
    _fields_ = [ ("n1", c_int), ("n2", c_int), ("dx1", c_double), ("dx2", c_double), ("velocity", c_void_p), ("damping", c_void_p), ("exp_order", c_int), ("U", c_void_p), ("V", c_void_p), ("Del2", c_void_p) ]



class propagator1:
    
    def __init__( self ):
        self.data = None
    
    def free():
        if self.data != None:
            try:
                c_waveprop.method1_free( self.data )
            except:
                print "propagator1: free error."
        
    
    def set( a1=0.0, b1=None, a2=0.0, b2=None, n1=None, n2=None, dx1=None, velocity=lambda x,y:1.0, damping=lambda x,y:0.0, expansion_order=5, hdaf_m1=12, hdaf_m2=12, hdaf_gamma1=0.5, hdaf_gamma2=0.5 ):
        
        if b1 != None and b2 != None:
            if n1 != None and n2 != None:
                dx1 = float(b1-a1)/n1
                dx2 = float(b2-a2)/n2
            else:
                if dx1 != None and dx2 != None:
                    n1 = float(b1-a1)/dx1
                    n2 = float(b2-a2)/dx2
        else:
            if dx1 != None and dx2 != None:
                if n1 != None and n2 != None:
                    b1 = float(dx1)*n1 + a1
                    b2 = float(dx2)*n2 + a2
                else:
                    if b1 != None and b2 != None:
                        n1 = float(b1-a1)/dx1
                        n2 = float(b2-a2)/dx2
            else:
                if n1 != None and n2 != None:
                    if dx1 != None and dx2 != None:
                        b1 = float(dx1)*n1 + a1
                        b2 = float(dx2)*n2 + a2
                    else:
                        if b1 != None and b2 != None:
                            dx1 = float(b1-a1)/n1
                            dx2 = float(b2-a2)/n2
        
        try:
            C = numpy.zeros((n1,n2))
            i1=0
            while i1<n1:
                i2=0
                while i2<n2:
                    x1 = float(i1)*dx1+a1
                    x2 = float(i2)*dx2+a2
                    C[i1,i2] = velocity(x1,x2)
                    i2 += 1
                i1 += 1
            
            velocity = C
        except:
            try:
                if numpy.shape(velocity) == (n1,n2):
                    velocity = numpy.array(velocity)
                else:
                    print "propagator1: I don't know what to do with 'velocity'"
            except:
                print "propagator1: I don't know what to do with 'velocity'"
        
        try:
            D = numpy.zeros((n1,n2))
            i1=0
            while i1<n1:
                i2=0
                while i2<n2:
                    x1 = float(i1)*dx1+a1
                    x2 = float(i2)*dx2+a2
                    d[i1,i2] = damping(x1,x2)
                    i2 += 1
                i1 += 1
            
            damping = D
        except:
            try:
                if numpy.shape(damping) == (n1,n2):
                    damping = numpy.array(damping)
                else:
                    print "propagator1: I don't know what to do with 'damping'"
            except:
                print "propagator1: I don't know what to do with 'damping'"
        
        velocity = numpy.array(velocity)
        damping = numpy.array(damping)
        
        self.free()
        
        try:
            self.data = c_void_p()
            c_waveprop.method1_init( self.data, c_int(n1), c_int(n2) c_double(dx1), c_double(dx2), velocity.ctypes.data_as(c_void_p), damping.ctypes.data_as(c_void_p), c_int(expansion_order), c_int(hdaf_m1), c_int(hdaf_m2), c_double(hdaf_gamma1), c_double(hdaf_gamma2) )
            
        except:
            print "propagator1: initialization error."
    
    
    def __call__( time_step, ui, vi, uf, vf ):
        
        try:
            c_waveprop.method1_execute( self.data, c_double(time_step), ui.ctypes.data_as(c_void_p), vi.ctypes.data_as(c_void_p), uf.ctypes.data_as(c_void_p), vf.ctypes.data_as(c_void_p) )
        except:
            print "propagator1: propagation error."
    
    
    
    
if __name__ == "__main__":
    
    print c_waveprop
    
    
    
    
    












"""

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
2



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



        





