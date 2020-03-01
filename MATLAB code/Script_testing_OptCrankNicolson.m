% MATLAB Script to test EuOptCrankNicolson function and AmOptCrankNicolson
% with MATLAB in-built function blsprice and binprice


% Author: Ruinan Lu
%References:
%[1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.
%[2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.
%[3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.
% 
clear
clc
[c,p]=blsprice(50,50,0.1,5/12,0.4)
S0=50;
K=50;
r=0.1;
T=5/12;
sigma=0.4;
Smax=100;
ds=0.5;
dt=5/240;


p_n=EuOptCrankNicolson(S0,K,r,T,sigma,Smax,ds,dt,'put')
c_n=EuOptCrankNicolson(S0,K,r,T,sigma,Smax,ds,dt,'call')
c_n1=AmOptCrankNicolson(S0,K,r,T,sigma,Smax,ds,dt,'call')
p_nam=AmOptCrankNicolson(S0,K,r,T,sigma,Smax,ds,dt,'put')
[pr,pn_bin]=binprice(S0,S0,r,T,dt,sigma,0);
pn_bin(1,1)