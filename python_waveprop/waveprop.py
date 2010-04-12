
from types import *
from ctypes import *
from ctypes.util import *
import numpy


c_waveprop = cdll.LoadLibrary(find_library('waveprop'))

c_waveprop.method1_init.argtyprs = [ c_void_p, c_int, c_int, c_double, c_double, c_void_p, c_void_p, c_int, c_int, c_int, c_double, c_double ]

c_waveprop.method1_free.argtyprs = [ c_void_p ] 

c_waveprop.method1_execute.argtypes = [ c_void_p, c_double, c_void_p, c_void_p,  c_void_p,  c_void_p, ]


c_waveprop.method2_init.argtyprs = [ c_void_p, c_int, c_int, c_double, c_double, c_void_p, c_void_p, c_int, c_int, c_int, c_double, c_double ]

c_waveprop.method2_free.argtyprs = [ c_void_p ] 

c_waveprop.method2_execute.argtypes = [ c_void_p, c_double, c_void_p, c_void_p,  c_void_p,  c_void_p, c_void_p,  c_void_p ]


c_waveprop.render_png_scalar.argtypes = [ c_char_p, c_int, c_int, c_void_p,  c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double  ]




"""
int render_png_scalar(
    char *fname, int nx, int ny, double *data, double center, double major_scale,
    double red_midpoint, double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_midpoint, double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_midpoint, double blue_leftvalue, double blue_midvalue, double blue_rightvalue );
"""


def write_png( data, fname, major_scale=1.0, center=0., red_params=(0.5, 0.,0.5,1. ), green_params=(0.5, 0.,0.5,1. ), blue_params=(0.5, 0.,0.5,1. ) ):
    
    nx,ny = numpy.shape(data)
    
    try:
        err = c_waveprop.render_png_scalar( c_char_p(fname), c_int(nx), c_int(ny), data.ctypes.data_as(c_void_p), center, major_scale, c_double(red_params[0]),c_double(red_params[1]),c_double(red_params[2]),c_double(red_params[3]), c_double(green_params[0]),c_double(green_params[1]),c_double(green_params[2]),c_double(green_params[3]), c_double(blue_params[0]),c_double(blue_params[1]),c_double(blue_params[2]),c_double(blue_params[3]) )
    except:
        print "write_png error, c_waveprop.write_png, exception."



class grid2d:
    
    def __init__(self, a1=0.0, b1=None, a2=0.0, b2=None, n1=None, n2=None, dx1=None, dx2=None, ):
        
        if b1 != None and b2 != None:
            if n1 != None and n2 != None:
                dx1 = float(b1-a1)/n1
                dx2 = float(b2-a2)/n2
            elif dx1 != None and dx2 != None:
                n1 = float(b1-a1)/dx1
                n2 = float(b2-a2)/dx2
        elif dx1 != None and dx2 != None:
            if n1 != None and n2 != None:
                b1 = dx1*n1 + a1
                b2 = dx2*n2 + a2
            elif b1 != None and b2 != None:
                n1 = float(b1-a1)/dx1
                n2 = float(b2-a2)/dx2
        elif n1 != None and n2 != None:
            if dx1 != None and dx2 != None:
                b1 = float(dx1)*n1 + a1
                b2 = float(dx2)*n2 + a2
            elif b1 != None and b2 != None:
                dx1 = float(b1-a1)/n1
                dx2 = float(b2-a2)/n2
        
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.n1 = int(n1)
        self.n2 = int(n2)
        self.dx1 = float(dx1)
        self.dx2 = float(dx2)
    
    def coord_mesh(self ):
        return numpy.mgrid[ self.a1: self.a1+self.n1*self.dx1: self.dx1, self.a2: self.a2+self.n2*self.dx2: self.dx2 ]
    
    def zeros(self ):
        return numpy.zeros((self.n1,self.n2))
    
    def evaluate(self, f):
        F = self.zeros()
        i1=0
        while i1<self.n1:
            i2=0
            while i2<self.n2:
                x1 = self.dx1*i1+self.a1
                x2 = self.dx2*i2+self.a2
                F[i1,i2] = f(x1,x2)
                i2 += 1
            i1 += 1
        return F
    
    def shape(self ):
        return (self.n1, self.n2 )





class propagator1:
    
    def __init__(self ):
        self.data = None
        self.grid = None
    
    def free(self ):
        if self.data != None:
            if self.data.value != 0:
                try:
                    err = c_waveprop.method1_free( byref(self.data) )
                    if err != 0 or self.data.value != None:
                        print "propagator1: free error,", err, "\t", self.data.value
                except:
                    print "propagator1: free error, exception."
    
    
    def set( self, grid=None, a1=0.0, b1=None, a2=0.0, b2=None, n1=None, n2=None, dx1=None, dx2=None, velocity=lambda x,y:1.0, damping=lambda x,y:0.0, expansion_order=5, hdaf_m1=8, hdaf_m2=8, hdaf_gamma1=0.8, hdaf_gamma2=0.8 ):
        
        if grid != None:
            self.grid = grid
        else:
            self.grid = grid2d( a1=a1, b1=b1, a2=a2, b2=b2, n1=n1, n2=n2, dx1=dx1, dx2=dx2 )
        
        try:
            velocity = self.grid.evaluate(velocity)
        except:
            try:
                if numpy.shape(velocity) == self.grid.shape():
                    velocity = numpy.array(velocity)
                else:
                    print "propagator1: I don't know what to do with 'velocity'"
            except:
                print "propagator1: I don't know what to do with 'velocity'"
        
        try:
            damping = self.grid.evaluate(damping)
        except:
            try:
                if numpy.shape(damping) == self.grid.shape():
                    damping = numpy.array(damping)
                else:
                    print "propagator1: I don't know what to do with 'damping'"
            except:
                print "propagator1: I don't know what to do with 'damping'"
        
        velocity = numpy.array(velocity)
        damping = numpy.array(damping)
        
        self.free()
        
        try:
            self.data = c_void_p(0)
            err = c_waveprop.method1_init( byref(self.data), c_int(self.grid.n1), c_int(self.grid.n2), c_double(self.grid.dx1), c_double(self.grid.dx2), velocity.ctypes.data_as(c_void_p), damping.ctypes.data_as(c_void_p), c_int(expansion_order), c_int(hdaf_m1), c_int(hdaf_m2), c_double(hdaf_gamma1), c_double(hdaf_gamma2) )
            
            if err != 0 or self.data.value==0:
                print "propagator1: initialization error: ", err
            
        except:
            print "propagator1: initialization error, exception."
        
    
    def __call__( self, time_step, ui, vi, uf, vf, initial_force=None, final_force=None ):
        
        try:
            if initial_force is None and final_force is None:
                err = c_waveprop.method1_execute( self.data, c_double(time_step), ui.ctypes.data_as(c_void_p), vi.ctypes.data_as(c_void_p), uf.ctypes.data_as(c_void_p), vf.ctypes.data_as(c_void_p) )
                if err != 0:
                    print "propagator1: propagation error: ", err
            else:
                err = c_waveprop.method2_execute( self.data, c_double(time_step), ui.ctypes.data_as(c_void_p), vi.ctypes.data_as(c_void_p), uf.ctypes.data_as(c_void_p), vf.ctypes.data_as(c_void_p), initial_force.ctypes.data_as(c_void_p), final_force.ctypes.data_as(c_void_p) )
                if err != 0:
                    print "propagator1: propagation, force, error: ", err
        except:
            print "propagator1: propagation error, exception."
            quit()





if __name__ == "__main__" and True:
    
    
    #from enthought.mayavi import mlab
    from numpy import *
    import time as time_module
    from math import *
    dir="/workspace/output/scratch/"
    
    
    class interval:
        "interval (a,b] of the real numbers."
        
        def __init__(self, a,b ):
            self.a = a
            self.b = b
     
        def __contains__(self, x):
            return self.a < x and x <= self.b
    
    
    def indicator_function(A):
        return lambda x: float(x in A)
    
    norm = numpy.linalg.norm
    
    print "\n\n"    
    #------------------------------
    
    
    
    
    
    
    def u0(x,y):
        s = norm((x-3.,y))
        a = 0.5
        if s <= .5*a:
            return cos(2.0*pi*s/a)+1.0
        else:
            return 0
    
    def velocity(x,y):
        return 10.0
    
    def damping(x,y):
        ax = 0.4
        ay = 0.4
        return 1500.0*( 1.0/(cosh((x-10)/ax)) +  1.0/(cosh((x+10)/ax)) + 1.0/(cosh((y-10)/ay)) +  1.0/(cosh((y+10)/ax))  )
    
    def forcing_term(t ):
        return lambda x,y: 3.0*u0(x,y)*cos(2.*pi*t/2)
    
    
    grid = grid2d( n1 = 256, n2 = 256, a1 = -10, b1 = 10, a2 = -10, b2 = 10 )
    X,Y = grid.coord_mesh()
    
    velocity = grid.evaluate( velocity )
    damping = grid.evaluate( damping )

    write_png( velocity, dir+ "velocity.png" )
    write_png( damping, dir+ "damping.png" )
    
    P = propagator1()
    P.set( grid, expansion_order=7, velocity=velocity, damping=damping )
    
    final_time = 10.0
    time_step = 0.01
    
    
    u = grid.evaluate( u0 )
    v = grid.zeros()
    time=0
    count = 0
    
    write_png( u, dir+ "%04d.png"%count, major_scale=0.1, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )
    count += 1
    
    b_times = []
    c_times = []
    
    while time<final_time:
    
        f1 = None
        f2 = None
        """
        t1 = time_module.clock()
        f1 = grid.evaluate( forcing_term(time) )
        f2 = grid.evaluate( forcing_term(time+time_step) )
        t2 = time_module.clock()
        c_times.append(t2-t1)
        print "\t", numpy.mean(c_times)
        """
        
        
        
        t1 = time_module.clock()
        
        P( time_step, u,v, u,v, f1,f2 )
        
        t2 = time_module.clock()
        
        b_times.append(t2-t1)
        print numpy.mean(b_times)
        
        time += time_step
        
        write_png( u, dir+ "%04d.png"%count, major_scale=0.05, red_params=(0.5,0,0,1.0), green_params=(0.5,0,.7,0), blue_params=(0.5,1.0,0,0,) )
        count += 1
    
    
    P.free()
    
    
    







"""

if __name__ == "__main__" and True:
    
    from numpy import *
    from math import *
    dir="/workspace/output/scratch/"
    
    
    def color_map( s ):
        s *= 10.0
        r = atan(s)*2.0/pi+0.5;
        return (r,r,r)
    
    grid = grid2d( n1 = 256, n2 = 256, a1 = -10, b1 = 10, a2 = -10, b2 = 10 )
    
    
    def test(x,y):
        s = 0
        
        if 7 < x and x < 9 and -1 < y and y < 1:
            s += 1
        
        if -1 < x and x < 1 and 7 < y and y < 9:
            s -= 1
        
        return s
    
    
    
    
    test = grid.evaluate( test )
    write_png( test, color_map, dir+ "test.png" )

"""
















"""

    from enthought.mayavi import mlab
    from numpy import *
    from time import *
    print "\n\n"
    
    
    grid = grid2d( dx1 = 0.1, dx2 = 0.1, a1 = -3, b1 = 3, a2 = -3, b2 = 3 )
    X,Y = grid.coord_mesh()
    
    vel = grid.evaluate( lambda x,y: x*x*y/12  )
    u0 = grid.evaluate( lambda x,y: exp(-x*x-y*y) )
    v0 = grid.zeros()
    final_time = 1
    time_step = 0.2
    
    P = propagator1()
    P.set( grid, expansion_order=5, velocity=vel )
    
    u = array(u0)
    v = array(v0)
    time=0
    
    plot = mlab.surf(X,Y,vel)
    #mlab.show()
    
    while time<final_time:
        print time
        P( time_step, u,v, u,v )
        time += time_step
        sleep(0.3)
        plot.mlab_source.scalars = vel
    
    
    P.free()
"""





