function fx=demo_fun_explict(x)
% This is a demo function of specifing boundary conditions for Heat
% Equation
% Please indicate Piecewise functions by if else statements
% 
%Author: Ruinan Lu
%References:
%[1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.
%[2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.
%[3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.
% 
if x>=0 && x<=0.5
    fx=2*x;
elseif x>=0.5 &&x<=1
    fx=2*(1-x);
else
    error("This is the sample function,don't make farce")
end