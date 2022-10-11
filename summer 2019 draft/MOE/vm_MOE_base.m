%%% solves the weird Vm equation and shows the zeros,
%%% M is hardcoded

clear all;

% M = 10;
EM = [10; 20; 40; 50; 100; 120; 150; 200];
linS = {'-','--',':','.-'}; % linestyles

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

    % plot(log10(x),y)
    plot(log10(Vm),(y),linS{linS_ind},'Linewidth',2)
    hold on
    % plot([(Vm(1)) (Vm(end))], [0 0], '--')
    % plot([log10(x(2)) log10(x(end))], [0 0],'b--')
    xlabel('V_m')
    ylabel('g(V_m)')
    yline(0)
    axis([3 6 -500 1500])
%     axis([0 5.5e5 -500 1500])
    % legend('M=0','M=20','M=50','M=100','M=150','M=200')
    % title('Solution of V_m equation')

    options = optimset('Display','off');

    % Vm = fsolve(fun,4e6,options)
end

legend('M = 10','M = 20','M = 40','M = 50','M = 100','M = 120','M = 150','M = 200');
