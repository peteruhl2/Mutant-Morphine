clear all
% close all

Em = 0:300;

Rw = zeros(length(Em),1);
Rm = zeros(length(Em),1);
Rmax = zeros(length(Em),1);

% Q = zeros(length(Em),1);
% R = zeros(length(Em),1);

for k = 1:length(Em)

    M = Em(k);
    
    Mh = 100;
    rc = 0.16;
    rm = 0.52;
    qc = 1.23e-6;
    qm = 0.25;
    n = 7.8731;

    eta_r = (M.^n)./(Mh^n + M.^n);
    eta_q = 1-eta_r;

    r(k) = rc + (rm-rc)*eta_r;
    q(k) = qc + (qm - qc)*eta_q;

    lambda = 3690;%3690;
    F = 0.2;%0.2;
    bl = 1e-9;
    bh = 1e-7;
    p = 2500; %2500
    b = 0.15;%0.005 Vitaly: 0.01 to 0.4
    B = 15;%10.5;%50
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

    Tl = (lambda*(q(k)+dt))/(dt*(q(k)+r(k)+dt));
    Th = (lambda*r(k))/(dt*(q(k)+r(k)+dt));
    C = omega/dc;

    wtop = -p*(Th*bh*Ep + Tl*bl*Ep - Th*bh - Tl*bl);
    wbottom = (C*b + di)*dv;
    
    mtop = -p*(B*F*Th*bh + B*F*Tl*bl - B*Th*bh - B*Tl*bl + F*Th*bh +...
        F*Tl*bl - Th*bh - Tl*bl);
    mbottom = (B*di + C*b +di)*dv;
    
    Rw(k) = wtop./wbottom;
    Rm(k) = mtop./mbottom;
    Rmax(k) = max(Rw(k),Rm(k));
    
%     Q(k) = q;
%     R(k) = r;
end

plot(Em,Rw,'b','Linewidth',2.5)
hold on
plot(Em,Rm,'r','Linewidth',2.5)
plot(Em,Rmax,'g--','Linewidth',2.5)
legend('R_w','R_m','R_0','Location','northwest')
plot([0 Em(end)],[1 1],'--')
xlabel('Morphine Concentration','FontSize',16)
ylabel('Basic Reproduction Number','FontSize',16)
% title('Basic Reproduction Number')

[Rw(1) Rm(1)]
[Rw(end) Rm(end)]
