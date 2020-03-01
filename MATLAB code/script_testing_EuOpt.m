% MATLAB Script to test EuOptImplicit and EuOptExplicit function
% with MATLAB in-built function blsprice 


% Author: Ruinan Lu
%References:
%[1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.
%[2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.
%[3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.
% 
clc
clear
[call_bls,put_bls]=blsprice(50,50,0.1,5/12,0.4)
p_n=EuOptImplicit(50,50,0.1,5/12,0.4,100,0.5,5/2400,'put')
c_n=EuOptImplicit(50,50,0.1,5/12,0.4,100,0.5,5/2400,'call')
call_fd=EuOptExplicit(50,50,0.1,5/12,0.4,100,0.5,5/240,'call')
put_fd=EuOptExplicit(50,50,0.1,5/12,0.4,100,20.5,5/240,'put')