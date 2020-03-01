dx=0.001;
dt=0.001;
tmax=100*dt;
sol=Heat_equation_Explicit(0,dx,1,dt,tmax,'demo_fun_explict',0,0);
subplot(2,2,1);
plot(0:dx:1,sol(:,2))
axis([0 1 0 1])
subplot(2,2,2);
plot(0:dx:1,sol(:,11))
axis([0 1 0 1])
subplot(2,2,3);
plot(0:dx:1,sol(:,51))
axis([0 1 0 1])
subplot(2,2,4)
plot(0:dx:1,sol(:,101))
axis([0 1 0 1])