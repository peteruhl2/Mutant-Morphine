%%% M_thres surface over fitness and escape ratio values
%%% 100x100 takes ~18 hours on dimo computer
%%% got it down quite a bit by using parfor for the inner loop
%%% 150x150 took 14h 44min 
clear all

syms M
assume(M > 0)

npts = 100;

% %%% test
% escape = [0.3026 2.028 4.041];
% fitness = [0.106 0.4045 0.78]; 

%%% in case I want these values again
% escape = linspace(0.2547046,50,npts);
% fitness = linspace(0.000035,0.95332901,npts);

%%% actually want to go from 0 to 50 for B and 
%%% 0 to 1 for F
escape = linspace(0,15,npts);
fitness = linspace(0,1,npts);
results = zeros(length(escape),length(fitness));

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
% B = 30;%20.5;%10.5;%50
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

omega_base = 15;
psi = 0.1;
omega = omega_base*exp(-psi.*M);%50

%%% IFE values
Tl = (lambda*(q + dt))./(dt*(q+r+dt));
Th = lambda*r./(dt*(q+r+dt));
C = omega/dc;

tic
for i = 1:5 %length(escape)
    B = escape(i);
    parfor j = 1:length(fitness)
        
        [i j]
        F = fitness(j);

        %%% Rw and Rm
%         Rw = ( (1-ep)*(bh*Th + bl*Tl)*p )./( dv*(b*C+di) );
%         Rm = ( (1-F)*(bh*Th + bl*Tl)*(1+B)*p)./(dv*(di*B+b*C+di));

        eq = ( (1-EP)*(bh*Th + bl*Tl)*p )./( dv*(b*C+di) ) == ( (1-F)*(bh*Th + bl*Tl)*(1+B)*p)./(dv*(di*B+b*C+di));
%         M_switch = double(subs(vpasolve(eq,M)));
%         results(i,j) = M_switch;
        M_switch = double(subs(solve(eq,M)));
        results(i,j) = NaN;
        
        %%% NEED THIS IF USING solve instead of VPASOLVE 
        %%% Need to find the real solution
        %%% if there is none then M_r = NaN
        for k = 1:length(M_switch)
            if imag(M_switch(k)) == 0
                results(i,j) = M_switch(k);
            end
        end
    end
end
toc

[X Y] = meshgrid(escape,fitness);
contourf(X,Y,(results)')
shading interp
xlabel('Escape ratio(B)')
ylabel('Fitness cost (F)')
zlabel('M_{Thres}')
colorbar
% caxis([0 80])
% axis([0 50 0 1])
title('M_{thresh} in F-B parameter space')

