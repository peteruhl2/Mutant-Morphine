%%% function that calculates Vm MOE value given morphine value
%%% AND F !!!

function Vm = vm_fit(M,F)

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

ep = 3e-5;
mu = 1;%4/24;
eta = 1;

alp = 6.7e-5;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;

omega_base = 15;
psi = 0.1;
omega = omega_base*exp(-psi.*M);%50

ALP = alp/(gamma + xi*M);

BL = (1-F)*bl;
BH = (1-F)*bh;

%%% MOE values 
%%% I did this a long time ago but its right I guess
fun = @(x) BL.*x.*(r.*lambda./((q+BH.*x+dt).*(r+BL.*x+dt)-r.*q)).*...
    ((q+BH.*x+dt)./r) ...
    +BH.*x.*(r.*lambda./((q+BH.*x+dt).*(r+BL.*x+dt)-r.*q))...
    -b./(1+B).*(dv.*x./p).*(omega./(dc-ALP.*(dv.*x./p)))...
    -di.*(dv.*x./p);

options = optimset('Display','off');
%%% CHECK THE SOLVER FOR GUESS AT ANSWER
Vm = fsolve(fun,1e5,options);

end