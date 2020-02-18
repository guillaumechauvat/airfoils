#===================================================#
# airfoils.py                                       #
#                                                   #
# author: Guillaume Chauvat                         #
# email: chauvat@kth.se                             #
# created 2020-02-10, last modification: 2020-02-18 #
#===================================================#

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import neksuite as nek


"""
Class representing an airfoil as closed 2D curve that can be interpolated in different ways
"""
class Airfoil:

    def __init__(self, x, y, n=10000, tol=1e-13, verbose=False):
        self.n = n
        self.tol = tol
        self.verbose = verbose
        self.x = x
        self.y = y
        (self.t, self.length) = self.compute_abscissas()
        self.x_int = interpolate.CubicSpline(self.t, self.x)
        self.y_int = interpolate.CubicSpline(self.t, self.y)
        self.dx_int = self.x_int.derivative(1)
        self.dy_int = self.y_int.derivative(1)
        self.d2x_int = self.x_int.derivative(2)
        self.d2y_int = self.y_int.derivative(2)
        self.trigonometric_direction = self.is_trigonometric()
        
    
    """
    Computes the normalized curvilinear abscissas t corresponding to the coordinates of the airfoil.
    """
    def compute_abscissas(self):
        # first, approximate the curve by linear interpolation to get an approximation of the curvilinear abscissa.
        lengths = np.sqrt(np.diff(self.x)**2 + np.diff(self.y)**2)
        s = np.concatenate(([0], lengths)).cumsum()  # curvilinear abscissa
        t = s/s[-1]  # normalize
        
        # compute the new curvilinear abscissa and iterate until convergence
        tc1 = np.linspace(0, 1, self.n)
        err = 1.
        i = 0
        if self.verbose:
            print("Computing the curvilinear abscissa...")
        while err > self.tol:
            # get the high resolution points
            x1 = interpolate.CubicSpline(t, self.x)(tc1)
            y1 = interpolate.CubicSpline(t, self.y)(tc1)
            # approximate the curvilinear abscissa based on a linear interpolation between high resolution points
            lengths = np.sqrt(np.diff(x1)**2 + np.diff(y1)**2)
            s1 = np.concatenate(([0], lengths)).cumsum()
            t1 = s1/s1[-1]
            # update the normalized curvilinear abscissa at coordinate points
            t_new = interpolate.CubicSpline(tc1, t1)(t)
            err = np.linalg.norm(t_new-t)
            t = t_new
            if self.verbose:
                print("iteration", i, "err =", err)
                i = i + 1
        return (t, s1[-1])
        

    """
    Interpolates x: t->x(t) and y: t->y(t) at points in the array tc.
    """
    def pos(self, tc):
        return (self.x_int(tc), self.y_int(tc))


    """
    Returns the normalized curvilinear abscissa corresponding to a x/c location on the airfoil (where c is the chord).
    The last argument indicates which side of the airfoil must be returned, since there are always two possibilities for 0 < x < 1.
    
    x0: location where the curvilinear abscissa must be computed
    t0: first guess of the location
    dt0: initial step size for iterating.
    """
    def curvilinear_abscissa(self, x0, t0=0.8, dt0=0.1, tol=1e-13):
        # We want to interpolate t as a function of x, but this only works if the function x: t -> x(t) is bijective, which is not the case for an airfoil.
        # Instead, we can look for x iteratively. We also avoid running into problems where x'(t) is very small (in which case interpolating t as a function of x generates a lot of oscillations).
        # make sure t1 is in [0, 1]
        t1 = t0 + dt0
        if t1>1:
            t1 = t0 - dt0
        err = 1
        
        # initial values
        f0 = self.x_int(t0)-x0
        f1 = self.x_int(t1)-x0
        
        # iterate with a secant method
        i = 0
        if self.verbose:
            print("Computing t using the secant method...")
        while err > tol:
            t2 = t1 - f1*((t1-t0)/(f1-f0))
            f2 = self.x_int(t2)-x0
            f0 = f1
            f1 = f2
            t0 = t1
            t1 = t2
            err = abs(t1-t0)
            if self.verbose:
                print("iteration " + str(i) + ", err = " + str(err))
                i = i+1
        return t1

    
    """
    returns the normal vector to the surface at a given point
    """
    def normal(self, t0):
        dx = self.dx_int(t0)
        dy = self.dy_int(t0)
        l = np.sqrt(dx**2+dy**2)
        nx = dy/l
        ny = -dx/l
        if not self.trigonometric_direction:
            nx = -nx
            ny = -ny
        return (nx, ny)
    
    
    """
    returns the tangent vector to the surface at a given point
    """
    def tangent(self, t0):
        dx = self.dx_int(t0)
        dy = self.dy_int(t0)
        l = np.sqrt(dx**2+dy**2)
        nx = dx/l
        ny = dy/l
        return (nx, ny)

    
    """
    returns the derivative of the tangent vector to the surface at a given point
    """
    def dtangent(self, t0):
        # ||dx**2 + dy**2|| is approximately equal to self.length, but we want to be very accurate and account for small deviations here
        # so speed ~ self.length and fact1 << 1/self.length
        dx = self.dx_int(t0)
        dy = self.dy_int(t0)
        d2x = self.d2x_int(t0)
        d2y = self.d2y_int(t0)
        speed = np.sqrt(dx**2+dy**2)
        fact1 = (dx*d2x+dy*d2y)/speed**3
        return (d2x/speed-fact1*dx, d2y/speed-fact1*dy)
    
    
    """
    returns the closest point (x, y) on the surface to an arbitrary point (x0, y0)
    """
    def closest(self, x0, y0, t0=0.5, tol=1e-13, alpha=1):
        # t0 is the starting point
        # the variable alpha can be chosen between 0 and 1 to sacrifice speed for stability
        # find a zero of the function f: t -> (tangent(t)|x(t)-x0)
        # f': t: -> 1 + (tangent'(t)|x(t)-x0)

        t1 = t0
        (x1, y1) = self.pos(t1)
        (taux, tauy) = self.tangent(t1)
        f1 = taux*(x1-x0)+tauy*(y1-y0)
        err = 1
        i = 0
        if self.verbose:
            print("Projecting point (" + str(x0) + ", " + str(y0) + ") on airfoil...")
        while err > tol:
            (dtaux, dtauy) = self.dtangent(t1)
            dx = self.dx_int(t1)
            dy = self.dy_int(t1)
            df = dx*taux + dy*tauy + dtaux*(x1-x0) + dtauy*(y1-y0)
            t2 = t1 - f1/df
            err = abs(t2-t1)
            t1 = (1-alpha)*t1+alpha*t2
            (x1, y1) = self.pos(t1)
            (taux, tauy) = self.tangent(t1)
            f1 = taux*(x1-x0)+tauy*(y1-y0)
            if self.verbose:
                print("iteration " + str(i) + ", t = " + str(t1) + ", f = " + str(f1) + ", err = " + str(err))
                i = i+1
        return (t1, x1, y1)

    """
    checks whether the airfoil is described in the trigonometric direction
    """
    def is_trigonometric(self):
        dx1 = self.x[1]-self.x[0]
        dx2 = self.x[-2]-self.x[0]
        dy1 = self.y[1]-self.y[0]
        dy2 = self.y[-2]-self.y[0]
        return dx1*dy2-dx2*dy1 > 0
