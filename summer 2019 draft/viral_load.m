
M = linspace(0,200,50);
results = zeros(length(M),1);

for i = 1:length(M)
    results(i) = vr_sim(M(i));
end

hold on
plot(M,log10(results),'LineWidth',2)
xlabel('Morphine')
ylabel('log10 Set-point viral load')


function x = vr_sim(M)

global lambda q r bl bh F dt p dv
global b di B omega dc EP ALP 

Mh = 100; %2.8534e-3;
rc = 0.16;
rm = 0.52;
qc = 1.23e-6;
qm = 0.25;
n = 8;

eta_r = @(M) (M^n)/(Mh^n+M^n);
eta_q = @(M) 1-eta_r(M);

r = rc + (rm-rc)*eta_r(M);
q = qc + (qm - qc)*eta_q(M);

lambda = 3690;%3690;
F = 0.03;%0.2;
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
omega_base = 15; %50
omega = omega_base*exp(-psi*M);%50

Th0 = 60650; %morphine
%     Th0 = 40980; %control
Tl0 = 1e6-Th0;%60650;
Vw0 = 200;
Vm0 = 0;
Iw0 = 0;
Im0 = 0;
C0 = 0; %10;

y0 = [Tl0 Th0 Vw0 Vm0 Iw0 Im0 C0];

options = odeset('NonNegative',1);

tspan = [0 1e6];
% tspan = [0 400];

[t y] = ode15s(@mut_model,tspan,y0,options);

x = y(end,3) + y(end,4);
end

