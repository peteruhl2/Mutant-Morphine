clear all;

global lambda q r bl bh F dt p dv
global b di B omega dc EP ALP 

Morphine = linspace(0,200,500);
results = zeros(length(Morphine),1);
Vw = zeros(length(Morphine),1);
Vm = zeros(length(Morphine),1);
C = zeros(length(Morphine),1);

for i = 1:length(Morphine)
    i
    M = Morphine(i);

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
    %     Th0 = 40980; %control
    Tl0 = 1e6-Th0;%60650;
    Vw0 = 200;
    Vm0 = 0;
    Iw0 = 0;
    Im0 = 0;
    C0 = 0; %10;

    y0 = [Tl0 Th0 Vw0 Vm0 Iw0 Im0 C0];

    options = odeset('NonNegative',1);

    tspan = [0 1e6];
    % tspan = [0 400];

    [t y] = ode15s(@mut_model,tspan,y0,options);

%     plot(t,log10(y(:,3)+y(:,4)))
%     hold on
%     xlabel('Days post infection')
%     ylabel('log_{10} viral RNA per ml')

    % title('Individual viral populations, M = 0')
    
    %%% save total viral load at equilibrium
    results(i) = y(end,3) + y(end,4);
    Vw(i) = y(end,3);
    Vm(i) = y(end,4);
    C(i) = y(end,7);

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% IFE values
%     Tl = (lambda*(q + dt))./(dt*(q+r+dt));
%     Th = lambda*r./(dt*(q+r+dt));
%     C = omega/dc;
% 
%     %%% Rw and Rm
%     Rw = ( (1-EP).*(bh*Th + bl*Tl)*p )./( dv*(b*C+di) );
%     Rm = ( (1-F)*(bh*Th + bl*Tl)*(1+B)*p)./(dv*(di*B+b*C+di));
    
end

figure()
hold on; box on
plot(Morphine,log10(results),'Linewidth',2)
xlabel('Morphine Concentration')
ylabel('log_{10} Steady State Viral Load')

figure()
hold on; box on
plot(Morphine,(C),'Linewidth',2)
xlabel('Morphine Concentration')
ylabel('Steady state CTL count per mL')
title('d)                                                                                                         ')

figure()
hold on; box on
plot(Morphine,Vm./(Vw + Vm)*100,'Linewidth',2)
xlabel('Morphine Concentration')
% ylabel('Mutant % of Viral Load')
ylabel('$\frac{V_m}{V_w + V_m}$','Interpreter','latex','Fontsize',16)
axis([Morphine(1) Morphine(end) 0 120])
title('c)                                                                                                         ')