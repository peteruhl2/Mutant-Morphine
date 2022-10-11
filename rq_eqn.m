clear all

% q = 1.23e-6; %Morphine
% r = 0.52;

% q = 0.25; %control
% r = 0.16;

% Em = 0:0.001:0.1;
Em = linspace(0,200);
r = zeros(length(Em),1);
q = zeros(length(Em),1);

for k = 1:length(Em)

    M = Em(k);
    Mh = 100; %2.8534e-3;
    rc = 0.16;
    rm = 0.5;
    qc = 1.23e-6;
    qm = 0.25;
    n = 7.8731;

    eta_r = (M.^n)./(Mh^n + M.^n);
    eta_q = 1-eta_r;

    r(k) = rc + (rm-rc)*eta_r;
    q(k) = qc + (qm - qc)*eta_q;
end

hold on
plot(Em,r)
plot(Em,q)
xlabel('M')
legend('r','q')






