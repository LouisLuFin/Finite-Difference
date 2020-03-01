function price=AmOptCrankNicolson_optimizer_method(S0,K,r,T,sigma,Smax,ds,dt,catagory)
% MATLAB Function to price a American option using the Crank-Nicolson Scheme with
% MATLAB build-in optimizer 
% Sometimes Effective, sometimes takes an unacceptable amount of time to
% calculate the result
% 
% Author: Ruinan Lu
%References:
%[1]Brandimarte P. Numerical methods in finance and economics: a MATLAB-based introduction[M]. John Wiley & Sons, 2013.
%[2]Seydel R, Seydel R. Tools for computational finance[M]. Berlin: Springer, 2006.
%[3]Ramalho L. Fluent python: Clear, concise, and effective programming[M]. " O'Reilly Media, Inc.", 2015.
% 


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
intrinsic_value=matval(2:M,N+1);
%Set up the tridiagonal coefficients matrix
a=0.25*dt*(sigma^2*(veti.^2)-r*veti);
b=-0.5*dt*(sigma^2*(veti.^2)+r);
c=0.25*dt*(sigma^2*(veti.^2)+r*veti);
M1=-diag(a(3:M),-1)+diag(1-b(2:M))-diag(c(2:M-1),1);
M2=diag(a(3:M),-1)+diag(1+b(2:M))+diag(c(2:M-1),1);
% set optimizer object options
opts = optimoptions('lsqnonlin');
opts.Display = 'none';
opts.MaxIterations = 50;
opts.FunctionTolerance = 0.01;
opts.OptimalityTolerance = 0.01;
opts.StepTolerance = 0.01;
opts.TypicalX = max(intrinsic_value,1e-4);
opt_cost = nan(1,N);
exitflag = nan(1,N);
aux=zeros(M-1,1);
for j=N:-1:1
    aux(1)=a(2)*(matval(1,j)+matval(1,j+1));
    aux(M-1)=c(M)*(matval(M+1,j)+matval(M+1,j+1));
    costFun=@(F) (M1*F-M2*matval(2:M,j+1)-aux);
    [Fopt,opt_cost(j),~,exitflag(j)] = lsqnonlin(costFun,matval(2:M,j+1),intrinsic_value,[],opts);
    matval(2:M,j)=Fopt;
end
price=interp1(vetS,matval(:,1),S0);

