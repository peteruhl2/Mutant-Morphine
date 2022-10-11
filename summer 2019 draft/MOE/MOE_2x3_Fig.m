%%% new script to make a big 2x3 figure of all the MOE value stuff
%%% 11/20/2021

%%% g(Vm) plot here =======================================================% M = 10;
EM = [10; 20; 40; 50; 100; 120; 150; 200];
% linS = {'-','--',':','.-'}; % linestyles
linS = {'-'};

for i = 1:length(EM)
    
    M = EM(i);
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
    F = 0.1;%0.2;
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
    ep = ep/(mu + eta*M); % new 11/19/21

    alp = 6.7e-5;%6.7e-6
    gamma = 1; %0.4;%0.4;
    xi = 1;

    omega_base = 15;
    psi = 0.1;
    omega = omega_base*exp(-psi.*M);%50

    ALP = alp/(gamma + xi*M);

    BL = (1-F)*bl;
    BH = (1-F)*bh;

    %%% MOE values 
    %%% I did this a long time ago but its right I guess
    fun = @(x) BL.*x.*(r.*lambda./((q+BH.*x+dt).*(r+BL.*x+dt)-r.*q)).*...
        ((q+BH.*x+dt)./r) ...
        +BH.*x.*(r.*lambda./((q+BH.*x+dt).*(r+BL.*x+dt)-r.*q))...
        -b./(1+B).*(dv.*x./p).*(omega./(dc-ALP.*(dv.*x./p)))...
        -di.*(dv.*x./p);

    % %%% trying the curve again 11/19/21
    % A = @(x) BL.*x.*lambda.*(x.*BH + q.*dt) + BH.*r.*x.*lambda;
    % BB = @(x) (x.*BH + q + dt).*(x.*BL + r + dt) - r*q;
    % C = @(x) dc - ALP*dv*x/p;

    % fun = @(x) A(x)./BB(x) - b*dv*x*omega./((1+B)*p*(C(x))) - di*dv*x/p;

    npts = 1000;
    Vm = linspace(0, 1e6, npts);
    % Vm = linspace(0, 2e5, npts); %%% use this for low M
    y = fun(Vm);

    %%% pick linestyle
    linS_ind = mod(randi([1,length(linS)]),length(EM));

    subplot(2,3,1)
    % plot(log10(x),y)
    plot(log10(Vm),(y),linS{linS_ind},'Linewidth',2)
    hold on
    % plot([(Vm(1)) (Vm(end))], [0 0], '--')
    % plot([log10(x(2)) log10(x(end))], [0 0],'b--')
    xlabel('log_{10} V_m')
    ylabel('g(V_m)')
    yline(0)
    axis([3 6 -500 1500])
%     axis([0 5.5e5 -500 1500])
    % legend('M=0','M=20','M=50','M=100','M=150','M=200')
    % title('Solution of V_m equation')

    options = optimset('Display','off');
    title('a)                           ')
    % Vm = fsolve(fun,4e6,options)
end

% legend('M = 10','M = 20','M = 40','M = 50','M = 100','M = 120','M = 150','M = 200');



%%% original Vm figure here ===============================================
M1 = 0:55;
Vm_results1 = zeros(length(M1),1);

for i=1:length(M1)
    Vm_results1(i) = vm_solver_2(M1(i));
end

M2 = 56:200;
Vm_results2 = zeros(length(M2),1);

for i=1:length(M2)
    Vm_results2(i) = vm_solver_2(M2(i));
end

%%% original 2x2 figure stuff here ========================================
M1 = 0:55;
results1 = zeros(length(M1),5); 

%%% q-r parameters
Mh = 100; %2.8534e-3;
rc = 0.16;
rm = 0.52;
qc = 1.23e-6;
qm = 0.25;
n = 8;

eta_r = @(M) (M.^n)./(Mh^n+M.^n);
eta_q = @(M) 1-eta_r(M);

%%% other parameters
lambda = 3690;%3690;
F = 0.1;%0.2;
bl = 1e-9;
bh = 1e-7;
p = 2500; %2500
b = 0.25;%0.005 Vitaly: 0.01 to 0.4
B = 30;%20.5;%10.5;%50
dt = 0.01;
dv = 23;
di = 0.7;
dc = 0.2; %0.63

BL = (1-F)*bl;
BH = (1-F)*bh;

ep = 3e-5;
mu = 1;%4/24;
eta = 1;

alp = 6.7e-5;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;


omega_base = 15;
psi = 0.1;
% omega = omega_base*exp(-psi.*M);%50

%%% MOE values stable part
for i=1:length(M1)
    r = rc + (rm - rc)*eta_r(M1(i));
    q = qc + (qm - qc)*eta_q(M1(i));
    ALP = alp./(gamma + xi*M1(i));
    omega = omega_base*exp(-psi.*M1(i));
    
    Vm1(i) = vm_solver_2(M1(i));
    Th1(i) = r*lambda/( (q+BH*Vm1(i)+dt)*(r+BL*Vm1(i)+dt)-r*q );
    Tl1(i) = Th1(i)*( (q+BH*Vm1(i)+dt)/r);
    Im1(i) = dv*Vm1(i)/p;
    C1(i) = omega/(dc - ALP*Im1(i));
end

results1 = [Tl1; Th1; Vm1; Im1; C1]';

%%% unstable part
M2 = 56:200;
results2 = zeros(length(M2),1);

for i=1:length(M2)
    r = rc + (rm - rc)*eta_r(M2(i));
    q = qc + (qm - qc)*eta_q(M2(i));
    ALP = alp./(gamma + xi*M2(i));
    omega = omega_base*exp(-psi.*M2(i));
    
    Vm2(i) = vm_solver_2(M2(i));
    Th2(i) = r*lambda/( (q+BH*Vm2(i)+dt)*(r+BL*Vm2(i)+dt)-r*q );
    Tl2(i) = Th2(i)*( (q+BH*Vm2(i)+dt)/r);
    Im2(i) = dv*Vm2(i)/p;
    C2(i) = omega/(dc - ALP*Im2(i));
end

results2 = [Tl2; Th2; Vm2; Im2; C2]';

%%% plot stuff here =======================================================

%%% original Vm here
% figure
subplot(2,3,4)
hold on
% plot(M,log10(results),'*','MarkerSize',2)
plot(M1,log10(Vm_results1),'b','Linewidth',2)
plot(M2,log10(Vm_results2),'b--','Linewidth',2)
xlabel('Morphine')
ylabel('V_m MOE value')
% title('V_m MOE values')
box on
title('b)                           ')


%%% original 2x2 here
subplot(2,3,2)
plot(M1,log10(Tl1),'b','Linewidth',2)
hold on
plot(M2,log10(Tl2),'b--','Linewidth',2)
xlabel('Morphine')
ylabel('log_{10}T_l MOE value')
title('c)                           ')

subplot(2,3,3)
plot(M1,log10(Th1),'b','Linewidth',2)
hold on
plot(M2,log10(Th2),'b--','Linewidth',2)
xlabel('Morphine')
ylabel('log_{10}T_h MOE value')
title('e)                           ')

subplot(2,3,5)
plot(M1,log10(Im1),'b','Linewidth',2)
hold on
plot(M2,log10(Im2),'b--','Linewidth',2)
xlabel('Morphine')
ylabel('log_{10}Im MOE value')
title('d)                           ')

subplot(2,3,6)
plot(M1,C1,'b','Linewidth',2)
hold on
plot(M2,C2,'b--','Linewidth',2)
xlabel('Morphine')
ylabel('C MOE value')
title('f)                           ')