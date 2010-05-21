
from types import *
from ctypes import *
from ctypes.util import *
from hdaf import freq_to_sigma, gamma_to_sigma
import numpy


c_waveprop = cdll.LoadLibrary(find_library('waveprop'))

c_waveprop.method1_init.argtyprs = [ c_void_p, c_int, c_int, c_double, c_double, c_void_p, c_void_p, c_int, c_int, c_int, c_double, c_double ]

c_waveprop.method1_free.argtyprs = [ c_void_p ]

c_waveprop.method1_execute.argtypes = [ c_void_p, c_double, c_void_p, c_void_p,  c_void_p,  c_void_p, ]


c_waveprop.method2_init.argtyprs = [ c_void_p, c_int, c_int, c_double, c_double, c_void_p, c_void_p, c_int, c_int, c_int, c_double, c_double ]

c_waveprop.method2_free.argtyprs = [ c_void_p ]

c_waveprop.method2_execute.argtypes = [ c_void_p, c_double, c_void_p, c_void_p,  c_void_p,  c_void_p, c_void_p,  c_void_p ]


c_waveprop.render_png_scalar.argtypes = [ c_char_p, c_int, c_int, c_void_p, c_int,  c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double  ]

c_waveprop.render_png_scalar_resample.argtypes = [ c_char_p, c_void_p, c_double, c_double, c_double, c_double, c_int, c_int, c_double, c_double, c_double, c_double, c_int, c_int, c_double, c_double, c_double, c_double, c_double, c_double,  c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double,  c_double, c_double, c_double ]

c_waveprop.resample2d_linear.argtypes = [ c_void_p, c_void_p, c_int, c_int, c_int, c_int, c_double, c_double, c_double, c_double, c_double, c_double, c_double ]

c_waveprop.cosh_damping.argtypes = [ c_void_p, c_int, c_int, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double ]



"""
int render_png_scalar(
    char *fname, int nx, int ny, double *data, double center, double major_scale,
    double red_midpoint, double red_leftvalue, double red_midvalue, double red_rightvalue,
    double green_midpoint, double green_leftvalue, double green_midvalue, double green_rightvalue,
    double blue_midpoint, double blue_leftvalue, double blue_midvalue, double blue_rightvalue );
"""




def hard_copy(ip):
    nx, ny = numpy.shape(ip)
    op = numpy.zeros((nx,ny))
    for i in range(nx):
        for j in range(ny):
            op[i][j] = ip[i][j]
    return op


def rotate_array(ip):
    nx, ny = numpy.shape(ip)
    op = numpy.zeros((ny,nx))
    for i in range(nx):
        for j in range(ny):
            op[j][nx-i-1] = ip[i][j]
    return op




class linearResampler2d:

    def __init__(self,  in_grid=None, out_grid=None, in_n1=None, in_n2=None, out_n1=None, out_n2=None, in_dx1=None, in_dx2=None, out_dx1=None, out_dx2=None, in_off1=None, in_off2=None, out_off1=None, out_off2=None ):

        if in_grid is not None:
            self.in_grid = in_grid
        else:
            self.in_grid = grid2d( a1=in_off1, a2=in_off2, n1=in_n1, n2=in_n2, dx1=in_dx1, dx2=in_dx2 )

        if out_grid is not None:
            self.out_grid = out_grid
        else:
            self.out_grid = grid2d( a1=out_off1, a2=out_off2, n1=out_n1, n2=out_n2, dx1=out_dx1, dx2=out_dx2 )
    
    def __call__(self, input ):

        result = self.out_grid.zeros()

        try:
            err = c_waveprop.resample2d_linear( input.ctypes.data_as(c_void_p), result.ctypes.data_as(c_void_p), c_int(self.in_grid.n1), c_int(self.in_grid.n2), c_int(self.out_grid.n1), c_int(self.out_grid.n2), c_double(self.in_grid.dx1), c_double(self.in_grid.dx2), c_double(self.out_grid.dx1), c_double(self.out_grid.dx2), c_double(self.in_grid.a1), c_double(self.in_grid.a2), c_double(self.out_grid.a1), c_double(self.out_grid.a2) )
            if err != 0:
                print "error: linearResampler2d call error, error: ", err
        except:
            print "error: linearResampler2d call error, exception."
            quit()

        return result





def write_png( data, fname, major_scale=1.0, center=0., red_params=(0.5, 0.,0.5,1. ), green_params=(0.5, 0.,0.5,1. ), blue_params=(0.5, 0.,0.5,1. ), ordering='rm' ):

    nx,ny = numpy.shape(data)
    ordering_int = 0
    if ordering not in [ 'rm', 'row_major', 'row' ]:
        ordering_int = 1

    try:
        err = c_waveprop.render_png_scalar( c_char_p(fname), c_int(nx), c_int(ny), data.ctypes.data_as(c_void_p), c_int(ordering_int), center, major_scale, c_double(red_params[0]),c_double(red_params[1]),c_double(red_params[2]),c_double(red_params[3]), c_double(green_params[0]),c_double(green_params[1]),c_double(green_params[2]),c_double(green_params[3]), c_double(blue_params[0]),c_double(blue_params[1]),c_double(blue_params[2]),c_double(blue_params[3]) )
    except:
        print "write_png error, c_waveprop.write_png, exception."


def write_png_resample( data, fname, grid_new, grid_old,  major_scale=1.0, center=0., red_params=(0.5, 0.,0.5,1. ), green_params=(0.5, 0.,0.5,1. ), blue_params=(0.5, 0.,0.5,1. ), default_color=(0.0,0.0,0.0) ):

    nx,ny = numpy.shape(data)

    try:
        err = c_waveprop.render_png_scalar_resample( c_char_p(fname), data.ctypes.data_as(c_void_p), grid_old.a1, grid_old.b2, grid_old.a2, grid_old.b2, grid_old.n1, grid_old.n2, grid_new.a1, grid_new.b2, grid_new.a2, grid_new.b2, grid_new.n1, grid_new.n2, center, major_scale, c_double(red_params[0]),c_double(red_params[1]),c_double(red_params[2]),c_double(red_params[3]), c_double(green_params[0]),c_double(green_params[1]),c_double(green_params[2]),c_double(green_params[3]), c_double(blue_params[0]),c_double(blue_params[1]),c_double(blue_params[2]),c_double(blue_params[3]), c_double(default_color[0]), c_double(default_color[1]), c_double(default_color[2]) )
    except:
        print "write_png_resample error, c_waveprop.render_png_scalar, exception."


class grid2d:

    def __init__(self, a1=None, b1=None, a2=None, b2=None, n1=None, n2=None, dx1=None, dx2=None, L1=None, L2=None ):
        
        if L1 is not None:
            a1=0.0;
            b1=L1;
        
        if L2 is not None:
            a2=0.0;
            b2=L2;
        
        if a1 is not None and b1 is not None and dx1 is not None:
            n1 = float(b1-a1)/dx1
        elif a1 is not None and b1 is not None and n1 is not None:
            dx1 = float(b1-a1)/n1
        elif a1 is not None and dx1 is not None and n1 is not None:
            b1 = dx1*n1 + a1
        elif b1 is not None and dx1 is not None and n1 is not None:
            a1 = b1 - dx1*n1
        else:
            print "error: grid2d, not enough parameters specified to construct grid, in first coordinates."
            quit()
        
        if a2 is not None and b2 is not None and dx2 is not None:
            n2 = float(b2-a2)/dx2
        elif a2 is not None and b2 is not None and n2 is not None:
            dx2 = float(b2-a2)/n2
        elif a2 is not None and dx2 is not None and n2 is not None:
            b2 = dx2*n2 + a2
        elif b2 is not None and dx2 is not None and n2 is not None:
            a2 = b2 - dx2*n2
        else:
            print "error: grid2d, not enough parameters specified to construct grid, in second coordinates."
            quit()
        
        self.a1 = float(a1)
        self.a2 = float(a2)
        self.b1 = float(n1*dx1+a1)
        self.b2 = float(n2*dx2+a2)
        self.n1 = int(n1)
        self.n2 = int(n2)
        self.dx1 = float(dx1)
        self.dx2 = float(dx2)

        [self.X1, self.X2] = self.coord_mesh()
    
    def debug(self, ):
        print "grid2d instance:", '\n', self
        print "[a1,b1) = ", self.a1, self.b1
        print "[a2,b2) = ", self.a2, self.b2
        print "n1,n2 = ", self.n1, self.n2
        print "dx1,dx2 = ", self.dx1, self.dx1
    
    def __call__(self, i,j ):
        return ( self.a1*self.dx1*i, self.a2*self.dx2*j )
    
    def coord_mesh(self ):
        return numpy.mgrid[ self.a1: self.a1+self.n1*self.dx1: self.dx1, self.a2: self.a2+self.n2*self.dx2: self.dx2 ]
    
    def zeros(self ):
        return numpy.zeros((self.n1,self.n2))
    
    def ones(self ):
        return numpy.ones((self.n1,self.n2))
    
    def evaluate(self, f):
        g = numpy.vectorize(f)
        try:
            return g(self.X1, self.X2)
        except:
            print "grid2d error: evaluate"
            return self.zeros()
    
    def shape(self ):
        return (self.n1, self.n2 )
    
    def index1(self, x ):
        return int((x-self.a1)/self.dx1)
    
    def index2(self, y ):
        return int((y-self.a2)/self.dx2)
    
    def index(self, X):
        return (self.index1(X[0]), self.index2(X[1]) )
    
    def extend_top(self, n_points=None, distance=None ):
        if n_points is not None:
            return grid2d( a1=self.a1, dx1=self.dx1, n1=self.n1,  a2=self.a2, dx2=self.dx2, n2=self.n2+n_points )
        elif distance is not None:
            return grid2d( a1=self.a1, dx1=self.dx1, n1=self.n1,  a2=self.a2, b2=self.b2+distance, dx2=self.dx2 )
    
    def extend_bottom(self, n_points=None, distance=None ):
        if n_points is not None:
            return grid2d( a1=self.a1, dx1=self.dx1, n1=self.n1,   b2=self.b2, dx2=self.dx2, n2=self.n2+n_points )
        elif distance is not None:
            return grid2d( a1=self.a1, dx1=self.dx1, n1=self.n1,   a2=self.a2-distance, b2=self.b2, dx2=self.dx2 )
    
    def extend_left(self, n_points=None, distance=None ):
        if n_points is not None:
            return grid2d( b1=self.b1, dx1=self.dx1, n1=self.n1+n_points,   a2=self.a2, b2=self.b2, dx2=self.dx2 )
        elif distance is not None:
            return grid2d( a1=self.a1-distance, b1=self.b1, dx1=self.dx1,   a2=self.a2, b2=self.b2, dx2=self.dx2 )
    
    def extend_right(self, n_points=None, distance=None ):
        if n_points is not None:
            return grid2d( a1=self.a1, dx1=self.dx1, n1=self.n1+n_points,   a2=self.a2, b2=self.b2, dx2=self.dx2 )
        elif distance is not None:
            return grid2d( a1=self.a1, dx1=self.dx1, b1=self.b1+distance,   a2=self.a2, b2=self.b2, dx2=self.dx2 )
    


class propagator1:

    def __init__(self ):
        self.data = None
        self.grid = None

    def __del__(self ):
        self.free()

    def free(self ):
        if self.data != None:
            if self.data.value != 0:
                try:
                    err = c_waveprop.method1_free( byref(self.data) )
                    if err != 0 or self.data.value != None:
                        print "propagator1: free error,", err, "\t", self.data.value
                except:
                    print "propagator1: free error, exception."


    def set( self, grid=None, a1=0.0, b1=None, a2=0.0, b2=None, n1=None, n2=None, dx1=None, dx2=None, velocity=lambda x,y:1.0, damping=lambda x,y:0.0, expansion_order=5, hdaf_order=12, hdaf_gamma=0.8 ):

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
            
            dx = max([self.grid.dx1, self.grid.dx2])
            sigma = gamma_to_sigma( hdaf_gamma, hdaf_order, dx )
            
            err = c_waveprop.method1_init( byref(self.data), c_int(self.grid.n1), c_int(self.grid.n2), c_double(self.grid.dx1), c_double(self.grid.dx2), velocity.ctypes.data_as(c_void_p), damping.ctypes.data_as(c_void_p), c_int(expansion_order), c_int(hdaf_order), c_int(hdaf_order), c_double(sigma), c_double(sigma) )

            if err != 0 or self.data.value==0:
                print "propagator1: initialization error: ", err

        except:
            print "propagator1: initialization error, exception."


    def __call__( self, time_step, ui, vi, uf, vf, initial_force=None, final_force=None ):
        
        err=None
        
        try:
            if initial_force is None and final_force is None:
                #print ui.ctypes.data_as(c_void_p), vi.ctypes.data_as(c_void_p), uf.ctypes.data_as(c_void_p), vf.ctypes.data_as(c_void_p)
                err = c_waveprop.method1_execute( self.data, c_double(time_step), ui.ctypes.data_as(c_void_p), vi.ctypes.data_as(c_void_p), uf.ctypes.data_as(c_void_p), vf.ctypes.data_as(c_void_p) )
                if err != 0:
                    print "propagator1: propagation error: ", err
                else:
                    err=None
            else:
                err = c_waveprop.method2_execute( self.data, c_double(time_step), ui.ctypes.data_as(c_void_p), vi.ctypes.data_as(c_void_p), uf.ctypes.data_as(c_void_p), vf.ctypes.data_as(c_void_p), initial_force.ctypes.data_as(c_void_p), final_force.ctypes.data_as(c_void_p) )
                if err != 0:
                    print "propagator1: propagation, force, error: ", err
                else:
                    err=None
        except:
            print "propagator1: propagation error, exception."
            err = 'excep'
        
        return err


 

def cosh_damping(grid, x1, x2, y1, y2, xs, ys, amp ):
    
    err=None
    damp=None
    
    try:
        damp=grid.zeros()
        err = c_waveprop.cosh_damping( damp.ctypes.data_as(c_void_p), c_int(grid.n1), c_int(grid.n2), c_double(grid.a1), c_double(grid.a2), c_double(grid.dx1), c_double(grid.dx2), c_double(x1), c_double(x2), c_double(y1), c_double(y2), c_double(xs), c_double(ys), c_double(amp) )
        
        if err != 0:
            print "error: cosh_damping, error: ", err
        else:
            err=None
    except:
        print "error: cosh_damping, exception."
    
    return damp







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

        t1 = time_module.clock()

        f1 = grid.evaluate( lambda x,y: 3.0*exp(-x*x-y*y)*cos(2.*pi*time/2) ) #forcing_term(time)
        f2 = grid.evaluate( lambda x,y: 3.0*exp(-x*x-y*y)*cos(2.*pi*(time+time_step)/2) ) #forcing_term(time+time_step)

        t2 = time_module.clock()
        c_times.append(t2-t1)
        print "\t", numpy.mean(c_times)



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





