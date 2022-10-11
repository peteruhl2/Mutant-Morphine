%%% calculate M_thres as function of escape (B)
clear all

syms M

npts = 100;

% fitness = linspace(0.000035,0.95332901,npts);
fitness = linspace(0.0000005,0.9329,npts);
% fitness = 0.2;
results = zeros(1,length(fitness));

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
parfor i = 1:length(fitness)
    i
    F = fitness(i);

    %%% Rw and Rm
    Rw = ( (1-EP)*(bh*Th + bl*Tl)*p )./( dv*(b*C+di) );
    Rm = ( (1-F)*(bh*Th + bl*Tl)*(1+B)*p)./(dv*(di*B+b*C+di));

    eq = Rw == Rm;
    M_switch = double(subs(vpasolve(eq,M)));
    results(i) = M_switch;
%     M_switch = double(subs(solve(eq,M)));
%     results(i) = NaN;
    
    %%% READ dont need anymore
    %%% Need to find the real solution
    %%% if there is none then M_r = NaN
%     for k = 1:length(M_switch)
%         if imag(M_switch(k)) == 0
%             results(i) = M_switch(k);
%         end
%     end
end
toc

hold on
plot(fitness,results)
xlabel('Fitness cost (F)')
ylabel('Threshold morphine concentration (M_{thresh})')
% title('M_{thresh} vs. fitness cost of mutation')