import numpy as np
import warnings


class BS_FDGrid:
    typecode = 'd'
    """ 
    A Python class to solve Black-Sholes Partial Differenial Equation with Finite Difference method and calculate the price of European and American equity vanilla option.
    This class could implement all three common Schemes of Finite Difference method: Explicit, Implicit, and Crank-Nicolson method
    Pleae find the the details of the three methods in the readme file
    The main part of this class is a numpy matrix that represents the grid, 
    and the number at every discreted time point could be exacted by 1-demension or 2-demonstion index, but only after the instance has completed calculation.
    With 1-D index x, the BS_FDGrid object returns the solution at time x
    With 2-D index specified as list [x,y], the BS_FDGrid object returns the solution at time y and the xth boundary condition
    This class could be directly demonstrated with print() function, it would show the price of the option 
    With Emplicit method, the solution may have unstable solutions, in this case, there will be a warning
    The file contains no error checking 
    
    Sample code to use these codes:
        b = BS_FDGrid(Smax, ds, T, dt)
        #European options
        b.SetBoundaryConditions(S0, K, r, sigma, 'call')
        b.CalGrid()
        print(b)

        c=Am_BS_FDGrid(Smax, ds, T, dt)
        #American options
        c.SetBoundaryConditions(S0, K, r, sigma, 'put')
        c.CalGrid()
        print(c)
    
    Author: Ruinan Lu
    References:
    [1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.
    [2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.
    [3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.
    
    """

    def __init__(self, Smax, ds, T, dt):
        """
        The constructor of the class specifies the following variables to initialise a Finite Difference grid:
        1) Smax: the maximum price the underlying equity could get in one's assumption
        2) ds: step scale for price
        3) T: the maturity of the option, specified as years
        4) dt: time step scale 
        """
        self.Smax = Smax
        self.ds = ds
        self.T = T
        self.dt = dt
        self.N = round(self.T/self.dt)
        self.M = round(self.Smax/self.ds)
        self.dt = self.T/self.N
        self.ds = self.Smax/self.M
        self.__FDGrid = np.zeros((self.M+1, self.N+1))

    def __getitem__(self, idx_list):
        """
        Special method to realise extract solution by index.
        Supports either 1 or 2-D indexing
        """
        if isinstance(idx_list, int):
            x = idx_list
            return self.FDGrid_full[:, x]
        if len(idx_list) == 2:
            x, y = idx_list
            return self.FDGrid_full[x, y]

    def __repr__(self):
        """Special method to output the grid with print() function"""
        return np.str(self.price)

    def SetBoundaryConditions(self, S0, K, r, sigma, catagory='call'):
        """
        Class method to input boundary conditions
        Input the details of the option in this method
        Parameters: SetBoundaryConditions(self, S0, K, r, sigma, catagory='call')
        S0: the price of the option at t0
        K: the strike price of the option
        r: risk-free rate 
        sigma: volatility of the underlying stock
        catagory: 'call' or 'put', 'call' by default
        """
        self.S0 = S0
        self.K = K
        self.r = r
        self.sigma = sigma
        self.catagory = catagory
        self.vetS = np.arange(0, self.Smax+self.ds, self.ds)
        self.veti = np.arange(0, self.M+1)
        self.vetj = np.arange(0, self.N+1)
        if catagory == 'call':
            self.__FDGrid[:, self.N] = np.maximum(self.vetS-self.K, 0)
            self.__FDGrid[0, :] = 0
            self.__FDGrid[self.M, :] = (
                self.Smax-self.K)*np.exp(-self.r*self.dt*(self.N-self.vetj))
        if catagory == 'put':
            self.__FDGrid[:, self.N] = np.maximum(self.K-self.vetS, 0)
            self.__FDGrid[0, :] = self.K * \
                np.exp(-self.r*self.dt*(self.N-self.vetj))
            self.__FDGrid[self.M, :] = 0

        self.FDGrid_ViewBoundary = self.__FDGrid

    def CalGrid(self, method='CrankNicolson'):
        """
        The method to calculate the grid, specify method by 'Explicit' 'Implicit' or 'CrankNicolson', use Crank-Nocolson method by default
        """
        aux = aux = np.zeros(self.M-1)

        if method == 'Explicit':
            self.a = 0.5*self.dt*(self.sigma**2*self.veti-self.r)*self.veti
            self.b = 1-self.dt*(self.sigma**2*(self.veti**2)+self.r)
            self.c = 0.5*self.dt*(self.sigma**2*self.veti+self.r)*self.veti
            coeff = np.diagflat(
                self.b[1:self.M])+np.diagflat(self.a[2:self.M], -1)+np.diagflat(self.c[1:self.M-1], 1)

            for jdx in range(self.N, 0, -1):
                self.__FDGrid[1:self.M, jdx -
                              1] = coeff.dot(self.__FDGrid[1:self.M, jdx])

            warnings.warn(
                "Answers given by explicit method may have a stability issue.")

        if method == 'Implicit':
            self.a = 0.5*(self.r*self.dt*self.veti -
                          self.sigma**2*self.dt*(self.veti**2))
            self.b = 1+self.sigma**2*self.dt*(self.veti**2)+self.r*self.dt
            self.c = -0.5*(self.r*self.dt*self.veti +
                           self.sigma**2*dt*(self.veti**2))
            coeff = np.diagflat(
                self.b[1:self.M])+np.diagflat(self.a[2:self.M], -1)+np.diagflat(self.c[1:self.M-1], 1)

            for jdx in range(self.N, 0, -1):
                aux = self.get_aux(jdx, aux)
                self.__FDGrid[1:self.M, jdx-1] = np.linalg.inv(
                    coeff).dot(self.__FDGrid[1:self.M, jdx]+aux)

        if method == 'CrankNicolson':
            self.a = 0.25*self.dt * \
                (self.sigma**2*(self.veti**2)-self.r*self.veti)
            self.b = -0.5*self.dt*(self.sigma**2*(self.veti**2)+self.r)
            self.c = 0.25*self.dt * \
                (self.sigma**2*(self.veti**2)+self.r*self.veti)

            M1 = np.diagflat(1-self.b[1:self.M])-np.diagflat(
                self.a[2:self.M], -1)-np.diagflat(self.c[1:self.M-1], 1)
            M2 = np.diagflat(1+self.b[1:self.M])+np.diagflat(
                self.a[2:self.M], -1)+np.diagflat(self.c[1:self.M-1], 1)
            for jdx in range(self.N, 0, -1):
                aux = self.get_aux(jdx, aux)
                self.__FDGrid[1:self.M, jdx-1] = np.linalg.inv(
                    M1).dot(M2.dot(self.__FDGrid[1:self.M, jdx]+aux))

        self.price = np.interp(self.S0, self.vetS, self.__FDGrid[:, 0])
        self.FDGrid_full = self.__FDGrid

    def get_aux(self, idx, aux):
        """
        Inner function in the calculation process
        """
        aux[0] = -self.a[1]*self.__FDGrid[0, idx]
        aux[self.M-2] = -self.c[self.M-1]*self.__FDGrid[self.M, idx-1]
        return aux


class Am_BS_FDGrid(BS_FDGrid):
    """
    Inherited object from BS_FDGrid to calculate the price of American options by checking whether one should exercise the option at every descreted time 
    """

    def __init__(self, Smax, ds, T, dt):
        super().__init__(Smax, ds, T, dt)
        self.__FDGrid = self._BS_FDGrid__FDGrid

    def SetBoundaryConditions(self, S0, K, r, sigma, catagory='call'):
        """
        Same method to input boundaty consitions with changes suitable for American options
        No longer need to discount the boundary condition as option could be exercised at any time before maturity
        Record payoff of exercising for further calculation
        Parameters: SetBoundaryConditions(self, S0, K, r, sigma, catagory='call')
        S0: the price of the option at t0
        K: the strike price of the option
        r: risk-free rate 
        sigma: volatility of the underlying stock
        catagory: 'call' or 'put', 'call' by default
        """
        super().SetBoundaryConditions(S0, K, r, sigma, catagory=catagory)
        if catagory == 'call':
            self.__FDGrid[self.M, :] = self.Smax-self.K
            self.payoff = np.maximum(self.vetS[1:self.M]-self.K, 0)
        if catagory == 'put':
            self.__FDGrid[0, :] = self.K
            self.payoff = np.maximum(self.K-self.vetS[1:self.M], 0)

    def CalGrid(self, method='CrankNicolson'):
        """
        The method to calculate the grid, specify method by 'Explicit' 'Implicit' or 'CrankNicolson', use Crank-Nocolson method by default
        """
        aux = aux = np.zeros(self.M-1)

        if method == 'Explicit':
            """
            Explicit method calculation process
            """
            self.a = 0.5*self.dt*(self.sigma**2*self.veti-self.r)*self.veti
            self.b = 1-self.dt*(self.sigma**2*(self.veti**2)+self.r)
            self.c = 0.5*self.dt*(self.sigma**2*self.veti+self.r)*self.veti
            coeff = np.diagflat(
                self.b[1:self.M])+np.diagflat(self.a[2:self.M], -1)+np.diagflat(self.c[1:self.M-1], 1)

            for jdx in range(self.N, 0, -1):
                self.__FDGrid[1:self.M, jdx -
                              1] = coeff.dot(self.__FDGrid[1:self.M, jdx])
                self.__FDGrid[1:self.M, jdx - 1] = np.maximum(
                    self.__FDGrid[1:self.M, jdx - 1], self.payoff)
                warnings.warn(
                    "Answers given by explicit method may have a stability issue.")

        if method == 'Implicit':
            """
            Implicit method calculation process
            """
            self.a = 0.5*(self.r*self.dt*self.veti -
                          self.sigma**2*self.dt*(self.veti**2))
            self.b = 1+self.sigma**2*self.dt*(self.veti**2)+self.r*self.dt
            self.c = -0.5*(self.r*self.dt*self.veti +
                           self.sigma**2*dt*(self.veti**2))
            coeff = np.diagflat(
                self.b[1:self.M])+np.diagflat(self.a[2:self.M], -1)+np.diagflat(self.c[1:self.M-1], 1)

            for jdx in range(self.N, 0, -1):
                aux = self.get_aux(jdx, aux)
                self.__FDGrid[1:self.M, jdx-1] = np.linalg.inv(
                    coeff).dot(self.__FDGrid[1:self.M, jdx]+aux)
                self.__FDGrid[1:self.M, jdx - 1] = np.maximum(
                    self.__FDGrid[1:self.M, jdx - 1], self.payoff)

        if method == 'CrankNicolson':
            """
            Crank-Nicolson method calculation process
            """
            self.a = 0.25*self.dt * \
                (self.sigma**2*(self.veti**2)-self.r*self.veti)
            self.b = -0.5*self.dt*(self.sigma**2*(self.veti**2)+self.r)
            self.c = 0.25*self.dt * \
                (self.sigma**2*(self.veti**2)+self.r*self.veti)

            M1 = np.diagflat(
                1-self.b[1:self.M])-np.diagflat(self.a[2:self.M], -1)-np.diagflat(self.c[1:self.M-1], 1)
            M2 = np.diagflat(
                1+self.b[1:self.M])+np.diagflat(self.a[2:self.M], -1)+np.diagflat(self.c[1:self.M-1], 1)
            M1[0, 0] = M1[0, 0]-self.a[1]
            M1[self.M-2, self.M-2] = M1[self.M-2, self.M-2]-self.c[self.M-1]
            for jdx in range(self.N, 0, -1):
                aux[0] = self.a[1]*self.__FDGrid[0, jdx-1]
                aux[self.M-2] = self.c[self.M-1]*self.__FDGrid[self.M, jdx-1]
                self.__FDGrid[1:self.M, jdx-1] = np.linalg.inv(
                    M1).dot(M2.dot(self.__FDGrid[1:self.M, jdx]+aux))
                self.__FDGrid[1:self.M, jdx - 1] = np.maximum(
                    self.__FDGrid[1:self.M, jdx - 1], self.payoff)

        self.price = np.interp(self.S0, self.vetS, self.__FDGrid[:, 0])
        self.FDGrid_full = self.__FDGrid


if __name__ == "__main__":
    """
    main function to test code while writing the code
    """
    S0 = 50
    K = 50
    r = 0.1
    T = 5/12
    sigma = 0.4
    Smax = 100
    ds = 0.5
    dt = 5/240

    b = BS_FDGrid(Smax, ds, T, dt)
    b.SetBoundaryConditions(S0, K, r, sigma)
    b.CalGrid('Explicit')
    print(b)

    c = Am_BS_FDGrid(Smax, ds, T, dt)
    c.SetBoundaryConditions(S0, K, r, sigma, 'put')
    c.CalGrid()
    print(c)
