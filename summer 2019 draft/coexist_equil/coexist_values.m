%%% get coexistence eq values given M

function E = coexist_values(M)

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

%%% get Vw, Vm
x = coexist_sim(M);
Vw = x(1);
Vm = x(2);

%%% other values
%%% Th is pretty messy
dum1 = (q + bh*Vw + BH*Vm + dt)/r;
dum2 = r + bl*Vw + BL*Vm + dt;
Th = lambda/(dum1*dum2 - q);

Tl = (lambda + q*Th)/(r + bl*Vw + BL*Vm + dt);
Iw = dv*Vw/p;
Im = dv*Vm/p;
C = (omega)/(dc - ALP*(Iw+Im));

E = [Tl Th Vw Vm Iw Im C];
end