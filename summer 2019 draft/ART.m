%%% file to make ART figures

clear all;

global lambda q r bl bh F dt p dv
global b di B omega dc EP ALP 

%%% morphine pre and post ART
M = 200;
M_ART = 200;

t_treat = 250;

ef1 = 0.9;
ef2 = 0.9;

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
F = 0.1;%0.1;
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
Tl0 = 1e6-Th0;%60650;
Vw0 = 200;
Vm0 = 0;
Iw0 = 0;
Im0 = 0;
C0 = 0; %10;


%%% solve without ART
y0 = [Tl0 Th0 Vw0 Vm0 Iw0 Im0 C0];
options = odeset('NonNegative',1);
tspan = [0 t_treat];

[t, y] = ode15s(@mut_model,tspan,y0,options);

%%% solve with treatment ==================================================

% morphine parameters after treatment begins
r = rc + (rm-rc)*eta_r(M_ART);
q = qc + (qm - qc)*eta_q(M_ART);
EP = ep/(mu + eta*M_ART);
ALP = alp/(gamma + xi*M_ART);
omega = omega_base*exp(-psi*M_ART);%50

% integrase inhibitor
bl = (1 - ef1)*bl;
bh = (1 - ef1)*bh;

% protease inhibitor
p = (1 - ef2)*p;


y0 = y(end,:);
tspan = [t_treat t_treat+50];

[t2, y2] = ode15s(@mut_model,tspan,y0,options);

t = [t; t2];
y = [y; y2];

%%% plot viral load =======================================================

plot(t,log10(y(:,3)+y(:,4)),'Linewidth',2)
hold on
xlabel('Days post infection')
ylabel('log_{10} viral RNA per ml')
% axis([300 400 -1 7])
yline(log10(50),'--')
xline(t_treat)
% legend('M = 0 ug/l Post-ART', 'M = 200 ug/l During ART', 'Detection Level', 'Beginning of ART','fontsize',14)







%%% Functions =============================================================

%%% ode function
function yp = mut_model(t,y)

global lambda q r bl bh F dt p dv
global b di B omega dc EP ALP 

Tl = y(1);
Th = y(2);
Vw = y(3);
Vm = y(4);
Iw = y(5);
Im = y(6);
C = y(7);

yp = ones(7,1);

yp(1) = lambda + q*Th - r*Tl - bl*Vw*Tl - (1-F)*bl*Vm*Tl - dt*Tl;
yp(2) = r*Tl - q*Th - bh*Vw*Th - (1-F)*bh*Vm*Th - dt*Th;
yp(3) = p*Iw - dv*Vw;
yp(4) = p*Im - dv*Vm;
yp(5) = (1-EP)*(bl*Vw*Tl + bh*Vw*Th) - b*Iw*C - di*Iw;
yp(6) = EP*(bl*Vw*Tl + bh*Vw*Th) + (1-F)*bl*Vm*Tl + ...
          (1-F)*bh*Vm*Th - (b/(1+B))*Im*C - di*Im;
yp(7) = omega + ALP*(Iw+Im)*C - dc*C;

end