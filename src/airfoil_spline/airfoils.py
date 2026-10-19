#===================================================#
# airfoils.py                                       #
#                                                   #
# author: Guillaume Chauvat                         #
# created 2020-02-10, last modification: 2026-10-18 #
#===================================================#

import numpy as np
from scipy import interpolate, optimize


class Airfoil:
    """
    Class representing an airfoil as closed 2D curve {t -> (x(t), y(t)) for t in [0, 1]} using cubic splines.
    This allows for smooth interpolation of the coordinates as well as their derivatives.
    """

    def __init__(self, x, y, n=10000, tol=1e-13, verbose=False):
        self.n = n
        self.tol = tol
        self.verbose = verbose
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        (self.t, self.length) = self.compute_abscissas()
        self.x_int = interpolate.CubicSpline(self.t, self.x)
        self.y_int = interpolate.CubicSpline(self.t, self.y)
        self.dx_int = self.x_int.derivative(1)
        self.dy_int = self.y_int.derivative(1)
        self.d2x_int = self.x_int.derivative(2)
        self.d2y_int = self.y_int.derivative(2)
        self.trigonometric_direction = self.is_trigonometric()

    @staticmethod
    def _as_t(t):
        """Converts value to a 1D array."""
        return np.atleast_1d(np.asarray(t, dtype=float))


    def compute_abscissas(self):
        """
        Computes the normalized curvilinear abscissas t corresponding to the coordinates of the airfoil.
        """

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
                print(f"iteration {i}, err = {err}")
                i = i + 1
        return (t, s1[-1])

    def pos(self, tc):
        """
        Interpolates x: t->x(t) and y: t->y(t) at points in the array tc.
        """

        tc = self._as_t(tc)
        return np.stack([self.x_int(tc), self.y_int(tc)], axis=-1)

    def curvilinear_abscissa(self, x0, t0, tol=1e-13):
        """
        Returns the normalized curvilinear abscissa corresponding to a x/c location on the airfoil (where c is the chord).
        The choice of t0 is important since there are always two possibilities for 0 < x < 1.

        x0: location where the curvilinear abscissa must be computed
        t0: first guess of the location
        """

        x0 = self._as_t(x0)
        t0 = self._as_t(t0)

        # initial value
        f0 = self.x_int(t0)-x0

        # iterate with Newton-Raphson
        i = 0
        err = 1.0
        if self.verbose:
            print("Computing t...")
        while np.any(err > tol):
            df = self.dx_int(t0)
            t1 = t0 - f0/df
            f0 = self.x_int(t1)-x0
            err = abs(t1-t0)
            t0 = t1
            if self.verbose:
                print(f"iteration {i}, err = {err}")
                i = i+1
        return t1

    def normal(self, t0):
        """
        returns the normal vector to the surface at a given point
        """

        t0 = self._as_t(t0)
        dx = self.dx_int(t0)
        dy = self.dy_int(t0)
        length = np.sqrt(dx**2+dy**2)
        nx = dy/length
        ny = -dx/length
        if not self.trigonometric_direction:
            nx = -nx
            ny = -ny
        return np.stack([nx, ny], axis=-1)

    def tangent(self, t0):
        """
        returns the tangent vector to the surface at a given point
        """

        t0 = self._as_t(t0)
        dx = self.dx_int(t0)
        dy = self.dy_int(t0)
        length = np.sqrt(dx**2+dy**2)
        nx = dx/length
        ny = dy/length
        return np.stack([nx, ny], axis=-1)

    def dtangent(self, t0):
        """
        returns the derivative of the tangent vector to the surface at a given point
        """

        # ||dx**2 + dy**2|| is approximately equal to self.length, but we want to be very accurate and account for small deviations here
        # so speed ~ self.length and fact1 << 1/self.length
        t0 = self._as_t(t0)
        dx = self.dx_int(t0)
        dy = self.dy_int(t0)
        d2x = self.d2x_int(t0)
        d2y = self.d2y_int(t0)
        speed = np.sqrt(dx**2+dy**2)
        fact1 = (dx*d2x+dy*d2y)/speed**3
        return np.stack([d2x/speed-fact1*dx, d2y/speed-fact1*dy], axis=-1)

    def closest(self, x0, y0, t0, tol=1e-13, alpha=1):
        """
        returns the closest point (x, y) on the surface to an arbitrary point (x0, y0)
        """

        # t0 is the starting point
        # the variable alpha can be chosen between 0 and 1 to sacrifice speed for stability
        # find a zero of the function f: t -> (tangent(t)|x(t)-x0)
        # f': t: -> 1 + (tangent'(t)|x(t)-x0)

        t0 = self._as_t(t0)
        x0 = self._as_t(x0)
        y0 = self._as_t(y0)
        t1 = t0
        p1 = self.pos(t1)
        x1, y1 = p1[:, 0], p1[:, 1]
        tau = self.tangent(t1)
        taux, tauy = tau[:, 0], tau[:, 1]
        f1 = taux * (x1 - x0) + tauy * (y1 - y0)
        err = 1
        i = 0
        if self.verbose:
            print(f"Projecting point ({x0}, {y0}) on airfoil...")
        while np.any(err > tol):
            dtau = self.dtangent(t1)
            dtaux, dtauy = dtau[:, 0], dtau[:, 1]
            dx = self.dx_int(t1)
            dy = self.dy_int(t1)
            df = dx * taux + dy * tauy + dtaux * (x1 - x0) + dtauy * (y1 - y0)
            t2 = t1 - f1/df
            err = np.abs(t2 - t1)
            t1 = (1 - alpha) * t1 + alpha * t2
            p1 = self.pos(t1)
            x1, y1 = p1[:, 0], p1[:, 1]
            tau = self.tangent(t1)
            taux, tauy = tau[:, 0], tau[:, 1]
            f1 = taux * (x1 - x0) + tauy * (y1 - y0)
            if self.verbose:
                print(f"iteration {i}, t = {t1}, f = {f1}, err = {err}")
                i = i+1
        return t1

    def le(self, tol=1e-13, relax=0.9):
        """
        Find the position of the leading-edge, defined as the position where
        the airfoil surface is vertical.
        """

        # The LE should be around 0.5
        t0 = np.array([0.5])
        # we want ny = 0 (horizontal normal)

        def err(t):
            return self.normal(t)[:, 1]

        return optimize.fsolve(err, t0)


    def is_trigonometric(self):
        """
        checks whether the airfoil is described in the trigonometric direction
        """

        x, y = self.x, self.y
        area = 0.5 * np.sum(x[:-1]*y[1:] - x[1:]*y[:-1])
        return area > 0


def rotate(x, y, angle, xr, yr):
    """
    rotates a set of points around a point (xr, yr)
    """

    x1 = xr + (x-xr)*np.cos(angle) + (y-yr)*np.sin(angle)
    y1 = yr - (x-xr)*np.sin(angle) + (y-yr)*np.cos(angle)
    return (x1, y1)
