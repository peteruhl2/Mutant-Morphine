clear all; %close all;

M = [0:300];

Mh = 100; %2.8534e-3;
rc = 0.16;
rm = 0.52;
qc = 1.23e-6;
qm = 0.25;
n = 7.8731;
% n = 8;

eta_r = @(M) (M.^n)./(Mh^n+M.^n);
eta_q = @(M) 1-eta_r(M);

% eta_r = (M.^n)./(Mh^n+M.^n);
% eta_q = 1-eta_r;

r = rc + (rm - rc)*eta_r(M);
q = qc + (qm - qc)*eta_q(M);

hold on
plot(M,r,'r')
plot(M,q,'b')
xlabel('Morphine')
legend('r(M)','q(M)')