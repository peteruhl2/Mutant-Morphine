clear all

% q = 1.23e-6; %Morphine
% r = 0.52;

% q = 0.25; %control
% r = 0.16;

M = 0:0.001:0.01;

Mh = 2.8534e-3;
rc = 0.16;
rm = 0.5;
qc = 1.23e-6;
qm = 0.16;
n = 7.8731;

eta_r = @(M)(M.^n)./(Mh^n + M.^n);

eta_q = @(M) 1-eta_r(M);

q = qc + (qm - qc)*eta_q(M);

plot(M,q)
hold on