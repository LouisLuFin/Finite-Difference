function solution=Heat_equation_Implicit(xmin,dx,xmax,dt,tmax,f_initial,f_ub,f_lb)
% MATLAB Function to solve heat equation using the Implicit Scheme 
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
%% Using rho in the implicit equation would reduce computation
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
%% Calculating Parameter matrix B
di=(1+2*rho)*ones(1,N-1);
rho_arr=-rho*ones(1,N-2);
B=diag(di)+diag(rho_arr,1)+diag(rho_arr,-1);
%% Calculating the numbers in the grid
for j=1:M
    gj=[solution(1,j);zeros(N-3,1);solution(N+1,j)];
    solution(2:N,j+1)=B\(solution(2:N,j)+gj);
end