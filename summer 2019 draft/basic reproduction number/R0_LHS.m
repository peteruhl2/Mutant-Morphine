%%% Latin hypercube sampling for Rw, Rm
%%% this one actually makes the double bar graph for the figure

clear all;

M = 10;

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

% [x,y] = R0(200,par)
% return

%%% LHC Sampling here =====================================================
% bounds = [1500 10000;
%           1e-8 2.8;
%           0.005 2.7;
%           1e-11 1e-5;
%           1e-9 1e-3;
%           500 5500;
%           0.005 1.8;
%           0.1 100;
%           0.001 1.2;
%           1 50;
%           0.01 10;
%           0.001 1.6;
%           3e-7 3e-3;
%           0.001 0.99;
%           0.1 40;
%           0.01 50;
%           0.01 50;
%           0.001 1.5];

bounds = [0.01*par' 2*par'];

n_sims = 20000;
n_params = length(par);

% take samples
X = lhsdesign(n_sims,n_params);

% apply parameters to columns of X
for i = 1:n_params
%     X(:,i) = par(i)*X(:,i);
    X(:,i) = bounds(i,1) + (bounds(i,2) - bounds(i,1))*X(:,i);
end

%%% Rw, Rm loop ===========================================================
results = zeros(n_sims,2);

for i = 1:n_sims
    i
    par = X(i,:);
    
    [Rw, Rm] = R0(M,par);
    results(i,:) = [Rw, Rm];
end

%%% remover outliers ======================================================
max1 = max(rmoutliers(results(:,1)));
max2 = max(rmoutliers(results(:,2)));
sd1 = std(rmoutliers(results(:,1)));
sd2 = std(rmoutliers(results(:,2)));
real_max = max(max1,max2);
real_sd = max(sd1,sd2);

% keepvals = find(results(:,1) <= real_max + real_sd);
keepvals = find(results(:,1) <= max1);

results = results(keepvals,:);
X = X(keepvals,:);



%%% get correlations coefficients =========================================
corrs = zeros(n_params,2);
for i = 1:n_params
    temp1 = corrcoef(results(:,1), X(:,i));
    temp2 = corrcoef(results(:,2), X(:,i));
    
    corrs(i,1) = temp1(1,2);
    corrs(i,2) = temp2(1,2);
end

%%% bar plot ==============================================================

figure()
bar(corrs)
xticks([1:19])
xticklabels({'\lambda','q','r','\beta_l','\beta_{h}','p',...
       'b','B','\delta_T','\delta_V','\delta_I','\delta_C',...
       '\epsilon','F','\omega','\mu','\eta','\psi'})
ylim([-1 1])

% xtick font size
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
% yticks(-1:0.2:1)
xtickangle(0)
axis([0 19 -1.2 1.2])

legend('R_0^w','R_0^m')
title('b)                                                                                                                         ','Fontsize',14)

%%% Functions =============================================================

function [Rw,Rm] = R0(M,par)
% M = par(end);

% %%% q-r parameters
% Mh = 100; %2.8534e-3;
% rc = 0.16;
% rm = 0.52;
% qc = 1.23e-6;
% qm = 0.25;
% n = 8;
% 
% eta_r = @(M) (M.^n)./(Mh^n+M.^n);
% eta_q = @(M) 1-eta_r(M);
% 
% r = rc + (rm - rc)*eta_r(M);
% q = qc + (qm - qc)*eta_q(M);

lambda = par(1);
% q = (qc + (qm - qc)*eta_q(M))*par(2);
% r = (rc + (rm - rc)*eta_r(M))*par(3);
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
EP = ep./(mu + eta*M);
omega = omega_base*exp(-psi.*M);%50

%%% IFE values
Tl = (lambda*(q + dt))./(dt*(q+r+dt));
Th = lambda*r./(dt*(q+r+dt));
C = omega/dc;

%%% Rw and Rm
Rw = ( (1-EP).*(bh*Th + bl*Tl)*p )./( dv*(b*C+di) );
Rm = ( (1-F)*(bh*Th + bl*Tl)*(1+B)*p)./(dv*(di*B+b*C+di));
% R0 = max(Rw,Rm);

% [Rw Rm]

end


