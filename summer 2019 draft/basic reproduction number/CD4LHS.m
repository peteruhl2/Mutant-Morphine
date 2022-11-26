%%% Latin hypercube sampling for stead state CD4 count
clear all;

% global M

M = 00;

lambda = 3690;
q = 1;
r = 1;
bl = 1e-9;
bh = 1e-7;
p = 2500;
b = 0.25;
B = 30;
dt = 0.01;
dv = 23;
di = 0.7;
dc = 0.2;
ep = 3e-5;
F = 0.1;
omega_base = 15;
mu = 1;
eta = 1;
psi = 0.1;

par = [lambda, q, r, bl, bh, p, b, B, dt, dv, di, dc, ep, F, omega_base, mu, eta, psi];


%%% IC and timespan
Th0 = 60650; %morphine
Tl0 = 1e6-Th0;%60650;
Vw0 = 200;
Vm0 = 0;
Iw0 = 0;
Im0 = 0;
C0 = 0; %10;

y0 = [Tl0 Th0 Vw0 Vm0 Iw0 Im0 C0];

% tspan = [0 1e6];
% [t y] = ode15s(@(t,y) mut_model(t,y,par),tspan,y0);
% 
% Vw = y(end,3)
% Vm = y(end,4)
% 
% return

tspan = [0 1e6];

%%% LHC Sampling here =====================================================
bounds = [1500 10000;
          1e-8 2.8;
          0.005 2.7;
          1e-11 1e-5;
          1e-9 1e-3;
          500 5500;
          0.005 1.8;
          0.1 100;
          0.001 1.2;
          1 50;
          0.01 10;
          0.001 1.6;
          3e-7 3e-3;
          0.001 1.5;
          0.1 40;
          0.01 50;
          0.01 50;
          0.001 1.5];

% %%% better way to get bounds?
% bounds = [0.01*par' 2*par'];

n_sims = 10000;
n_params = length(par);

% take samples
X = lhsdesign(n_sims,n_params);

% apply parameters to columns of X
for i = 1:n_params
%     X(:,i) = par(i)*X(:,i);
    X(:,i) = bounds(i,1) + (bounds(i,2) - bounds(i,1))*X(:,i);
end

%%% Vw, Vm loop ===========================================================
results = zeros(n_sims,1);

tic
parfor i = 1:n_sims
    i
    par = X(i,:);
    
    [t y] = ode15s(@(t,y) mut_model(t,y,par,M),tspan,y0);
    
%     Vw = y(end,3);
%     Vm = y(end,4);
    CD4 = y(end,1) + y(end,2);
    results(i) = CD4;
end
toc

%%% remover outliers ======================================================
lose = isoutlier(results);
keep = ~lose;

results = results(keep);
X = X(keep,:);


%%% get correlations coefficients =========================================
corrs = zeros(n_params,1);
for i = 1:n_params
    temp1 = corrcoef(results(:,1), X(:,i));
%     temp2 = corrcoef(results(:,2), X(:,i));
    
    corrs(i) = temp1(1,2);
%     corrs(i,2) = temp2(1,2);
end

%%% bar plot ==============================================================

figure()
bar(corrs)
xticks([1:18])
xticklabels({'\lambda','q','r','\beta_l','\beta_{h}','p',...
       'b','B','\delta_T','\delta_V','\delta_I','\delta_C',...
       '\epsilon','F','\omega','\mu','\eta','\psi'})
ylim([-1 1])

% xtick font size
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
% yticks(-1:0.2:1)
xtickangle(0)

% legend('V_w','V_m')
title('c)                                                                                                                         ','Fontsize',14)
dim = [.7 .6 .3 .3];
str = "CD4+ PRCCs" + newline + "M = " + M + " ug/ml";
annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',14);










%%% Functions =============================================================

function yp = mut_model(t,y,par,M)
% global M

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

lambda = par(1);
q = par(2);
r = par(3);
bl = par(4);
bh = par(5);
p = par(6);
b = par(7);
B = par(8);
dt = par(9);
dv = par(10);
di = par(11);
dc = par(12);
ep = par(13);
F = par(14);
omega_base = par(15);
mu = par(16);
eta = par(17);
psi = par(18);

% morphine terms
alp = 6.7e-5;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;
ALP = alp/(gamma + xi*M);

EP = ep./(mu + eta*M);
omega = omega_base*exp(-psi.*M);%50

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