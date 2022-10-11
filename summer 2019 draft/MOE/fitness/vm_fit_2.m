%%% should get Vm for MOE by solving the model over a v long time
%%% and with F as an input

function x = vm_fit_2(M,fit)

global lambda q r bl bh dt p dv F
global b di B omega dc EP ALP  

F = fit;

Mh = 100; %2.8534e-3;
rc = 0.16;
rm = 0.52;
qc = 1.23e-6;
qm = 0.25;
n = 8;

eta_r = @(M) (M.^n)/(Mh.^n+M.^n);
eta_q = @(M) 1-eta_r(M);

r = rc + (rm-rc)*eta_r(M);
q = qc + (qm - qc)*eta_q(M);

lambda = 3690;%3690;
% F = 0.1;%0.2;
bl = 1e-9;
bh = 1e-7;
p = 2500; %2500
b = 0.25;%0.005 Vitaly: 0.01 to 0.4
B = 30; %30
dt = 0.01;
dv = 23;
di = 0.7;
dc = 0.2; %0.63

ep = 3e-5;
mu = 1;%4/24;
eta = 1;
EP = ep/(mu + eta*M);

alp = 6.7e-5;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;
ALP = alp/(gamma + xi*M);

psi = 0.1;
omega_base = 15;
omega = omega_base*exp(-psi*M);%50

Th0 = 60650; %morphine
%     Th0 = 40980; %control
Tl0 = 1e6-Th0;%60650;
Vw0 = 0;
Vm0 = 200; %%% THIS IS IMPORTANT FOR MOE
Iw0 = 0;
Im0 = 0;
C0 = 15; %10;

y0 = [Tl0 Th0 Vw0 Vm0 Iw0 Im0 C0];

options = odeset('NonNegative',1);

tspan = [0 1e6];

[t y] = ode15s(@mut_model,tspan,y0,options);

%%% get mutant viral load
x = y(end,4);

end