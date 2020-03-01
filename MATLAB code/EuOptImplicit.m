function price=EuOptImplicit(S0,K,r,T,sigma,Smax,ds,dt,catagory)
% MATLAB Function to price a European option using the Implicit Scheme 
% 
% Author: Ruinan Lu
%References:
%[1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.
%[2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.
%[3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.
% 
if (nargin<9)
    catagory='call';
end
%Set up grid and adjust increments if necessary
M=round(Smax/ds);
ds=Smax/M;
N=round(T/dt);
dt=T/N;
matval=zeros(M+1,N+1);
vetS=linspace(0,Smax,M+1)';
veti=0:M;
vetj=0:N;
%Set up boundary conditions
switch catagory
    case 'call'
        matval(:,N+1)=max(vetS-K,0);
        matval(1,:)=0;
        matval(M+1,:)=(Smax-K)*exp(-r*dt*(N-vetj));
    case 'put'
        matval(:,N+1)=max(K-vetS,0);
        matval(1,:)=K*exp(-r*dt*(N-vetj));
        matval(M+1,:)=0;
    otherwise
        error('catagory must be call or put')
end

%Set up the tridiagonal coefficients matrix
a=0.5*(r*dt*veti-sigma^2*dt*(veti.^2));
b=1+sigma^2*dt*veti.^2+r*dt;
c=-0.5*(r*dt*veti+sigma^2*dt*(veti.^2));
coeff=diag(a(3:M),-1)+diag(b(2:M))+diag(c(2:M-1),1);
[L,U]=lu(coeff);
%Solve the sequence of linear systems
aux=zeros(M-1,1);
for j=N:-1:1
    aux(1)=-a(2)*matval(1,j);
    aux(M-1)=-c(M)*matval(M+1,j);
    matval(2:M,j)=U\(L\(matval(2:M,j+1)+aux));
end
price=interp1(vetS,matval(:,1),S0);