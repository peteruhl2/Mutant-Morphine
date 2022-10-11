clear all

M = 35;
B = [10:100];
% B = 20;

Mh = 2.8534e-3;
rc = 0.16;
rm = 0.52;
qc = 1.23e-6;
qm = 0.25;
n = 7.8731;

eta_r = (M.^n)./(Mh^n + M.^n);
eta_q = 1-eta_r;

r = rc + (rm-rc)*eta_r;
q = qc + (qm - qc)*eta_q;

lambda = 3690;%3690;
% F = 0.2;%0.2;
bl = 1e-9;
bh = 1e-7;
p = 4000; %2500
b = 0.25;%0.005 Vitaly: 0.01 to 0.4
% B = 20.5;%10.5;%50
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

Ep = ep/(mu+eta*M);

Tl = (lambda*(q+dt))/(dt*(q+r+dt));
Th = (lambda*r)/(dt*(q+r+dt));
C = omega/dc;

%%%Rm=1, solved for F as function of B
top = B.*Th*bh*p + B.*Tl*bl*p - B.*di*dv - C*b*dv + Th*bh*p + Tl*bl*p - di*dv;
bottom = p*(B.*Th*bh + B.*Tl*bl + Tl*bh + Tl*bl);
F = top./bottom;

% %%%Rw
% Rw = -(p*(Th*bh*Ep + Tl*bl*Ep - Th*bh - Tl*bl))/...
%     ((C*b+di)*dv);

plot(B,F)
hold on
xlabel('B')
ylabel('F')

% plot([B(1) B(end)],[Rw Rw],'--')