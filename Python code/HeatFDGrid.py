
# %%
import numpy as np
import matplotlib.pyplot as plt
import warnings


class HeatFDGrid:
    """ 
    A Python class to solve Heat Partial Differenial Equation with Finite Difference method.
    This class could implement all three common Schemes of Finite Difference method: Explicit, Implicit, and Crank-Nicolson method
    Pleae find the the details of the three methods in the readme file
    The main part of this class is a numpy matrix that represents the grid, 
    and solution at every discreted time point could be exacted by 1-demension or 2-demonstion index, but only after the instance has completed calculation.
    With 1-D index x, the HeatFDGRid object returns the solution at time x
    With 2-D index specified as list [x,y], the HeatFDGrid object returns the solution at time y and the xth boundary condition
    This class supports users to input upper, lower and initial boundary conditions with self written function, 
    this file contains a sample initial boundary condition function, and upper and lower boundary condition are 0 if not specified 
    This class could be directly demonstrated with print() function, it would show the entire grid
    With Emplicit method, the solution may have unstable solutions, in this case, there will be a warning
    The file contains no error checking 
    
    Sample code to use these codes:
        a = HeatFDGrid(0, dx, 1, dt, 1000*dt)
        a.SetBoundaryConditions(demo_fun_explict)
        a.CalGrid()
    
    Author: Ruinan Lu
    References:
    [1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.
    [2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.
    [3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.
    """
    typecode = 'd'
    
    def __init__(self, xmin, dx, xmax, dt, tmax):
        """
        The constructor of the class specifies the following variables to initialise a Finite Difference grid:
        1) xmin: the minimum x, which is the number that will appear at the [0,0] location of the grid
        2) dx: boundary condition step scale 
        3) xmax: the maximum x, which is the number that will appear at the [-1,0] location of the grid
        4) dt: time step scale 
        5) tmax: the maximum time in the grid
        """
        self.xmin = xmin
        self.dx = dx
        self.xmax = xmax
        self.dt = dt
        self.tmax = tmax
        self.N = round((self.xmax-self.xmin)/self.dx)
        self.xmax = self.xmin+self.N*self.dx
        self.M = round(self.tmax/self.dt)
        self.__FDGrid = np.zeros((self.N+1, self.M+1))

    def __getitem__(self, idx_list):
        """Special method to realise extract solution by index.
        Supports either 1 or 2-D indexing"""
        if isinstance(idx_list, int):
            x = idx_list
            return self.FDGrid_full[:, x]
        if len(idx_list) == 2:
            x, y = idx_list
            return self.FDGrid_full[x, y]

    def __repr__(self):
        """Special method to output the grid with print() function"""
        self.FDGrid_repr = self.__FDGrid
        return np.array_str(self.FDGrid_repr)

    def SetBoundaryConditions(self, f_init, f_ub=0, f_lb=0):
        """
        Class method to input boundary conditions
        the first input is the initial consition of the equation, and is a must for the grid
        the second input is the upper bound of the grid, could be user specified function or number, would be 0 if not specified
        the third input is the lower bound  of the grid, could be user specified function or number, would be 0 if not specified
        """
        self.rho = self.dt/(self.dx**2)

        if callable(f_ub):
            self.__FDGrid[0, :] = np.array(
                list(map(f_ub, np.arange(0, self.tmax+self.dt, self.dt))))
        else:
            self.__FDGrid[0, :] = f_ub

        if callable(f_lb):
            self.__FDGrid[-1, :] = np.array(list(map(f_lb,
                                                     np.arange(0, self.tmax+self.dt, self.dt))))
        else:
            self.__FDGrid[-1, :] = f_ub

        self.__FDGrid[:, 0] = np.array(
            list(map(f_init, np.arange(self.xmin, self.xmax+self.dx, self.dx))))

        self.FDGrid_ViewBoundary = self.__FDGrid

    def CalGrid(self, method='CrankNicolson'):
        """
        The method to calculate the grid, specify method by 'Explicit' 'Implicit' or 'CrankNicolson', use Crank-Nocolson method by default
        """
        if method == 'Explicit':
            """
            Explicit method calculation process
            """
            di = (1-2*self.rho)*np.ones((1, self.N-1))
            rho_arr = self.rho*np.ones((1, self.N-2))
            self.A = np.diagflat(di)+np.diagflat(rho_arr, 1) + \
                np.diagflat(rho_arr, -1)
            stability_test(self.A)
            for jdx in range(0, self.M):
                gj = self.cal_gj(jdx)
                self.__FDGrid[1:self.N, jdx +
                              1] = self.A.dot(self.__FDGrid[1:self.N, jdx]+gj)

        if method == 'Implicit':
            """
            Implicit method calculation process
            """
            di = (1+2*self.rho)*np.ones((1, self.N-1))
            rho_arr = -self.rho*np.ones((1, self.N-2))
            self.B = np.diagflat(di)+np.diagflat(rho_arr, 1) + \
                np.diagflat(rho_arr, -1)

            for jdx in range(0, self.M):
                gj = self.cal_gj(jdx)
                self.__FDGrid[1:self.N, jdx+1] = np.linalg.inv(
                    self.B).dot(self.__FDGrid[1:self.N, jdx]+gj)

        if method == 'CrankNicolson':
            """
            Crank-Nicolson method calculation process
            """
            diC = 2*(1+self.rho)*np.ones((1, self.N-1))
            rho_arrC = -self.rho*np.ones((1, self.N-2))
            self.C = np.diagflat(diC)+np.diagflat(rho_arrC, 1) + \
                np.diagflat(rho_arrC, -1)
            diD = 2*(1-self.rho)*np.ones((1, self.N-1))
            rho_arrD = self.rho*np.ones((1, self.N-2))
            self.D = np.diagflat(diD)+np.diagflat(rho_arrD, 1) + \
                np.diagflat(rho_arrD, -1)

            for jdx in range(0, self.M):
                gj = self.cal_gj(jdx)
                gj_plus_one = self.cal_gj(jdx+1)
                self.__FDGrid[1:self.N, jdx+1] = np.linalg.inv(self.C).dot(
                    (self.D.dot(self.__FDGrid[1:self.N, jdx])+self.rho*(gj+gj_plus_one)))

        self.FDGrid_full = self.__FDGrid

    def cal_gj(self, idx):
        """
        Inner function in the calculation process
        """
        gj = np.concatenate(([self.__FDGrid[0, idx]], np.zeros(
            (self.N-3)), [self.__FDGrid[self.N, idx]]))
        return gj


def stability_test(A):
    """
    A function to perform stability test for Explicit method
    """
    if np.linalg.norm(A, np.inf) <= 1:
        print("Convergence analysis show that this method is stable")
    else:
        print("Stability of this explict method cannot be granteened")
        warnings.warn("Stability of this explict method cannot be granteened")


def demo_fun_explict(x):
    """
    This is a demo function of specifing boundary conditions
    please indicate Piecewise functions by if else statements
    """
    if 0 <= x <= 0.5:
        return 2*x
    elif x >= 0.5 or x <= 1:
        return 2*(1-x)


# %%
if __name__ == "__main__":
    """
    main function to test code while writing the code
    """
    dx = 0.1
    dt = 0.001
    a = HeatFDGrid(0, dx, 1, dt, 1000*dt)
    a.SetBoundaryConditions(demo_fun_explict)
    a.CalGrid()
    
    plt.subplot(221)
    plt.plot(np.arange(0, 1+dx, dx), a[1])
    plt.axis([0, 1, 0, 1])
    plt.subplot(222)
    plt.plot(np.arange(0, 1+dx, dx), a[:, 10])
    plt.axis([0, 1, 0, 1])
    plt.subplot(223)
    plt.plot(np.arange(0, 1+dx, dx), a[:, 50])
    plt.axis([0, 1, 0, 1])
    plt.subplot(224)
    plt.plot(np.arange(0, 1+dx, dx), a[:, 100])
    plt.axis([0, 1, 0, 1])
    print(a)

# %%
