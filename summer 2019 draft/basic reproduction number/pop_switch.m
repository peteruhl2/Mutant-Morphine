clear all;

M = 0:200;

%%% q-r parameters
Mh = 100; %2.8534e-3;
rc = 0.16;
rm = 0.52;
qc = 1.23e-6;
qm = 0.25;
n = 8;

eta_r = @(M) (M.^n)./(Mh^n+M.^n);
eta_q = @(M) 1-eta_r(M);

r = rc + (rm - rc)*eta_r(M);
q = qc + (qm - qc)*eta_q(M);

%%% other parameters
lambda = 3690;%3690;
F = 0.1;%0.2;
bl = 1e-9;
bh = 1e-7;
p = 2500; %2500
b = 0.25;%0.005 Vitaly: 0.01 to 0.4
B = 50;%20.5;%10.5;%50
dt = 0.01;
dv = 23;
di = 0.7;
dc = 0.2; %0.63

ep = 3e-5;
mu = 1;%4/24;
eta = 1;
EP = ep./(mu + eta*M);

alp = 6.7e-5;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;

omega_base = 15;
psi = 0.1;
omega = omega_base*exp(-psi.*M);%50

%%% IFE values
Tl = (lambda*(q + dt))./(dt*(q+r+dt));
Th = lambda*r./(dt*(q+r+dt));
C = omega/dc;

%%% Rw and Rm
Rw = ( (1-EP).*(bh*Th + bl*Tl)*p )./( dv*(b*C+di) );
Rm = ( (1-F)*(bh*Th + bl*Tl)*(1+B)*p)./(dv*(di*B+b*C+di));
R0 = max(Rw,Rm);

figure
hold on; box on
plot(M,Rw,'b',"Linewidth",2)
plot(M,Rm,'r',"Linewidth",2)
plot(M,R0,'g--',"Linewidth",2)
xlabel('Morphine concentration (ug/ml)','fontsize',14)
ylabel('Basic reproduction number','fontsize',14)
plot([M(1) M(end)],[1 1],'--')
legend('R_0^w','R_0^m','R_0','Location','Northwest','fontsize',12)
% title('Basic reproduction number')


[Rw(1) Rm(1)]
[Rw(end) Rm(end)]






