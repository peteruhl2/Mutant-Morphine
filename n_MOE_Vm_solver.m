%finding a numerical value of Vm need for the mutant only equilibrium
clear all; format long

M = 300;

Mh = 100;
rc = 0.16;
rm = 0.52;
qc = 1.23e-6;
qm = 0.25;
n = 7.8731;

eta_r = @(M) (M^n)/(Mh^n+M^n);
eta_q = @(M) 1-eta_r(M);

r = rc + (rm-rc)*eta_r(M);
q = qc + (qm - qc)*eta_q(M);

lambda = 3690;%3690;
F = 0.2;%0.2;
bl = 1e-9;
bh = 1e-7;
p = 2500; %2500
b = 0.25;%0.005 Vitaly: 0.01 to 0.4
B = 25;%10.5;%50
dt = 0.01;
dv = 23;
di = 0.3;
dc = 0.63; %0.63

ep = 3e-5;
mu = 1;%4/24;
eta = 1;

alp = 6.7e-6;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;

omega = 50*exp(-0.05*M);%50

ALP = alp/(gamma + xi*M);

BL = (1-F)*bl;
BH = (1-F)*bh;

% x = linspace(1e5, 3e6);
% x = linspace(0,3e6);

%%%Vm equation
fun = @(x) BL.*x.*(r.*lambda./((q+BH.*x+dt).*(r+BL.*x+dt)-r.*q)).*...
    ((q+BH.*x+dt)./r) ...
    +BH.*x.*(r.*lambda./((q+BH.*x+dt).*(r+BL.*x+dt)-r.*q))...
    -b./(1+B).*(dv.*x./p).*(omega./(dc-ALP.*(dv.*x./p)))...
    -di.*(dv.*x./p);

Vm = linspace(0, 3e6);
y = fun(Vm);


% plot(log10(x),y)
plot(Vm,y)
hold on
% plot([Vm(1) Vm(end)], [0 0], '--')
% plot([log10(x(2)) log10(x(end))], [0 0],'--')
xlabel('V_m')
% ylabel('g(V_m)')
axis([0 3e6 -2000 3000])
% legend('M=0','M=20','M=100')
% title('Solution of V_m equation')

options = optimset('Display','off');

Vm = fsolve(fun,4e6,options)

% format short










