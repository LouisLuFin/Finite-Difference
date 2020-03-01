function solution=Heat_equation_Explicit(xmin,dx,xmax,dt,tmax,f_initial,f_ub,f_lb)
% MATLAB Function to solve heat equation using the Explicit Scheme with
% stability testing
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

%% Using rho in the explict equation would reduce computation
rho=dt/dx^2;
vet_x=xmin:dx:xmax;
vet_t=0:dt:tmax;

%% stability testing
di=(1-2*rho)*ones(1,N-1);
rho_arr=rho*ones(1,N-2);
A=diag(di)+diag(rho_arr,1)+diag(rho_arr,-1);
norm_A=norm(A,inf);
if norm_A<=1
    disp("Convergence analysis show that this method is stable")
else
    warning("Stability of this explict method cannot be granteened")
end

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
%% Calculating the numbers in the grid
for j=1:M
    gj=[solution(1,j);zeros(N-3,1);solution(N+1,j)];
    solution(2:N,j+1)=A*solution(2:N,j)+gj;
end
%% Below is a calculation method with for loop

% % for j=1:M
% %     for i=2:N
% %         solution(i,j+1)=rho*solution(i-1,j)+(1-2*rho)*solution(i,j)+rho*solution(i+1,j);
% %     end
% % end
 end
    