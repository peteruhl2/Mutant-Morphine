clear all

Ef = [0:0.005:1];
% Ef = 0.2;
F = zeros(length(Ef),1);
M = zeros(length(Ef),1);

%%%q and r params, we can just use morphine values
r = 0.52;
q = 1.23e-6;

%most of the paramters
lambda = 3690;%3690;
bl = 1e-9;
bh = 1e-7;
p = 2000; %2500
b = 0.05;%0.005 Vitaly: 0.01 to 0.4
B = 7;%10.5;%50
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

omega = 50;

%%%steady state values
Tl = (lambda*(q+dt))/(dt*(q+r+dt));
Th = (lambda*r)/(dt*(q+r+dt));

wtop = -p*(Th*bh*(ep/(mu+eta*M)) + Tl*bl*(ep/(mu+eta*M)) - Th*bh - Tl*bl);
mtop = -p*(B*F*Th*bh + B*F*Tl*bl - B*Th*bh - B*Tl*bl + F*Th*bh +...
        F*Tl*bl - Th*bh - Tl*bl);

for i = 1:length(Ef)

    F = Ef(i);

    % C = ((omega*exp(-0.05*M))/dc);

    %%%basic reproduction equations
    fun = @(M) (-p*(Th*bh*(ep/(mu+eta*M)) + Tl*bl*(ep/(mu+eta*M)) - Th*bh - Tl*bl))/...
        ((((omega*exp(-0.05*M))/dc)*b + di)*dv) -...
        (-p*(B*F*Th*bh + B*F*Tl*bl - B*Th*bh - B*Tl*bl + F*Th*bh +...
        F*Tl*bl - Th*bh - Tl*bl))...
        /((B*di + ((omega*exp(-0.05*M))/dc)*b +di)*dv);
    
    options = optimset('Display','off');
    M(i) = fsolve(fun,100,options);
end

hold on
plot(Ef,M,'b','LineWidth',2.5);
xlabel('Fitness Cost of Mutation','FontSize',16)
ylabel('Threshold Morphine Concentration','FontSize',16)
% axis([0 1 0 300])
% plot([0.935 1],[0 0],'b','LineWidth',2.5)
set(gca,'box','on')


% M =solve(Rw==Rm,M);
