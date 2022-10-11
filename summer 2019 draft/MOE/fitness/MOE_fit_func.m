%%% calculates MOE values given one M value
%%% AND F

function E = MOE_fit_func(M,F)

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
% F = 0.2;%0.2;
bl = 1e-9;
bh = 1e-7;
p = 2500; %2500
b = 0.25;%0.005 Vitaly: 0.01 to 0.4
B = 30;%20.5;%10.5;%50
dt = 0.01;
dv = 23;
di = 0.7;
dc = 0.2; %0.63

BL = (1-F)*bl;
BH = (1-F)*bh;

ep = 3e-5;
mu = 1;%4/24;
eta = 1;
EP = ep/(mu + eta*M);

alp = 6.7e-5;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;
ALP = alp/(gamma + xi*M);

omega_base = 15;
psi = 0.1;
omega = omega_base*exp(-psi.*M);%50

Vm = vm_fit_2(M,F);
Th = r*lambda/( (q+BH*Vm+dt)*(r+BL*Vm+dt)-r*q );
Tl = Th*( (q+BH*Vm+dt) /r);
Im = dv*Vm/p;
C = omega/(dc - ALP*Im);

E = [Tl,Th,Vm,Im,C];

end