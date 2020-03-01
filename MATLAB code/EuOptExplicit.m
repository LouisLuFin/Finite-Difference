function price=EuOptExplicit(S0,K,r,T,sigma,Smax,ds,dt,catagory)
% MATLAB Function to price a European option using the Explicit Scheme 
% without proper vectoriztion and stability testing
% Author: Ruinan Lu
%References:
%[1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.
%[2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.
%[3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.
% 

if (nargin<9)
    catagory='call';
end
M=round(Smax/ds);
ds=Smax/M;
N=round(T/dt);
dt=T/N;
matval=zeros(M+1,N+1);
vetS=linspace(0,Smax,M+1)';
veti=0:M;
vetj=0:N;
%boundary conditions
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
%coefficients
a=0.5*dt*(sigma^2*veti-r).*veti;
b=1-dt*(sigma^2*veti.^2+r);
c=0.5*dt*(sigma^2*veti+r).*veti;
coeff=diag(a(3:M),-1)+diag(b(2:M))+diag(c(2:M-1),1);
for j=N:-1:1
    matval(
%backward deduction method
% for j=N:-1:1
%     for i=2:M
%         matval(i,j)=a(i)*matval(i-1,j+1)+b(i)*matval(i,j+1)+c(i)*matval(i+1,j+1);
%     end
% end
warning("Answers given by explicit method may have a stability issue.")
price=interp1(vetS,matval(:,1),S0);