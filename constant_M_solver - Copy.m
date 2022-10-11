clear all; 

global lambda q r bl bh F dt p dv
global ep mu eta b di B omega alp dc M gamma xi

Em = [0 150];

for M = Em 

    Mh = 2.8534e-3;
    rc = 0.16;
    rm = 0.52;
    qc = 1.23e-6;
    qm = 0.25;
    n = 7.8731;

    eta_r = @(M) (M^n)/(Mh^n+M^n);
    eta_q = @(M) 1-eta_r(M);

    r = rc + (rm-rc)*eta_r(M);
    q = qc + (qm - qc)*eta_q(M);

    % q = 1.23e-6; %Morphine
    % r = 0.52;

    % q = 0.25; %control
    % r = 0.16;

    lambda = 3690;%3690;
    F = 0.2;%0.2;
    bl = 1e-9;
    bh = 1e-7;
    p = 2500; %2500
    b = 0.15;%0.005 Vitaly: 0.01 to 0.4
    B = 20.5;%10.5;%50
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

    Tl0 = 40980;
    Th0 = 1e6-Tl0;%60650;
    Vw0 = 200;
    Vm0 = 00;
    Iw0 = 0;
    Im0 = 0;
    C0 = 0; %10;

    y0 = [Tl0 Th0 Vw0 Vm0 Iw0 Im0 C0];

    options = odeset('NonNegative',1);

    [t y] = ode15s(@cm_test_2,[0 300],y0,options);
    % [t y] = ode45(@cm_test_2,[0 180],y0,options);

    plot(t,log10(y(:,3)+y(:,4)))
    hold on
    xlabel('t (days)')
    ylabel('log10 viral load')
    
end

% legend('Control','Morphine')

% plot(t,log10(y(:,3)),'b')
% hold on
% plot(t,log10(y(:,4)),'r')
% legend('Wild-type', 'Mutant')