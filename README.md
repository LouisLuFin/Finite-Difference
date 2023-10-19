# A MATLAB and Python implementation of Finite Difference method for Heat and Black-Scholes Partial Differential Equation

These codes implement the numerical method of Finite Difference method to solve Heat PDE and Black-Scholes PDE. Specificly, the code for Black-Scholes PDE aims to price vanilla options such as European and American call and put. 

The algorithm is implemented in Python and MATLAB, and the Python code is in Object Oriented discipline and used Numpy to handle matrices. Also, both Python and MATLAB code allow users to write their own function to put into the code to set boundary conditions for the Finite difference grid. The code provided sample user-generated function for setting boundary condition for both Python code and MATLAB code. The Python objects also implemented  special methods in Python classes so as to make it sliceable (`__getitem__()`)and printable(`__repr__()`)

Please find detailed mathematical explanation in the README.pdf file as Github does not support Latex formulas. 

Author: Ruinan Lu

*References:*

*[1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.*

*[2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.*

*[3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.*

