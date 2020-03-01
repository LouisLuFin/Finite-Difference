function solution=Heat_equation_Crank_Nicolson(xmin,dx,xmax,dt,tmax,f_initial,f_ub,f_lb)
% MATLAB Function to solve heat equation using the Crank-Nocolson Scheme 
% 
% Author: Ruinan Lu
%References:
%[1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.
%[2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.
%[3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.
% 
%% Find the number of nodes on both axises and create a empty grid
N=round((xmax-xmin)/dx);
xmax=xmin+N*dx;
M=round(tmax/dt);
solution=zeros(N+1,M+1);
%% Using rho in the Crank Nicolson equation would reduce computation
rho=dt/dx^2;
vet_x=xmin:dx:xmax;
vet_t=0:dt:tmax;
%% Setting initial and boundary conditions
for i=1:N+1
    solution(i,1)=feval(f_initial,vet_x(i));
end
if f_ub==0
     solution(1,:)=0;
else
for j=1:M+1
    solution(1,j)=feval(f_ub,vet_t(j));
end
end

if f_lb==0
    solution(end,:)=0;
else
    for j=1:M
        solution(end,j)=feval(f_lb,vet_t(j));
    end
end
%% Calculating Parameter matrix C and D
diC=2*(1+rho)*ones(1,N-1);
rho_arrC=-rho*ones(1,N-2);
C=diag(diC)+diag(rho_arrC,1)+diag(rho_arrC,-1);
diD=2*(1-rho)*ones(1,N-1);
rho_arrD=rho*ones(1,N-2);
D=diag(diD)+diag(rho_arrD,1)+diag(rho_arrD,-1);
%% Calculating the numbers in the grid
for j=1:M
    gj=[solution(1,j);zeros(N-3,1);solution(N+1,j)];
    gj_plus_one=[solution(1,j+1);zeros(N-3,1);solution(N+1,j+1)];
    solution(2:N,j+1)=C\(D*solution(2:N,j)+rho*(gj+gj_plus_one));
end
end