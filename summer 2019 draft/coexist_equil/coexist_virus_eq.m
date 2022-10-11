%%% function that defines two equations to be solver for Vm and Vw for
%%% coexistence equilibrium
%%% might take M as input, we'll see

function F = coexist_virus_eq(x)

Vw = x(1);
Vm = x(2);

global M

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
F = 0.2;%0.2;
bl = 1e-9;
bh = 1e-7;
p = 2500; %2500
b = 0.25;%0.005 Vitaly: 0.01 to 0.4
B = 30;%20.5;%10.5;%50
dt = 0.01;
dv = 23;
di = 0.3;
dc = 0.63; %0.63

ep = 3e-5;
mu = 1;%4/24;
eta = 1;

EP = ep/(mu + eta*M);

alp = 6.7e-5;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;

omega_base = 50;
psi = 0.1;
omega = omega_base*exp(-psi.*M);%50

ALP = alp/(gamma + xi*M);

BL = (1-F)*bl;
BH = (1-F)*bh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% I c/p'ed all these into the equations 5 and 6 of the model

% Iw = ( dv*Vw/p );
% Im = ( dv*Vm/p );
% C = ( (omega + ALP*( ( dv*Vw/p ) + ( dv*Vm/p ) ))/dc );
% Th = ( lambda/( ( (q+bh*Vw + BH*Vm + dt)/r)*(r + bl*Vw + BL*Vm +dt) - q) );
% Tl = ( lambda + q*( lambda/( ( (q+bh*Vw + BH*Vm + dt)/r)*(r + bl*Vw + BL*Vm +dt) - q) ))/(r + bl*Vw + BL*Vm + dt);

%%% each term gets at least one line
%%% version before substiitutions
% F(1) = (1-EP)*(bl*Vw*Tl + bh*Vw*Th) - b*Iw*C - di*Iw;
F(1) = (1-EP)*(bl*Vw*( lambda + q*( lambda/( ( (q+bh*Vw + BH*Vm + dt)/r)*(r + bl*Vw + BL*Vm +dt) - q) ))/(r + bl*Vw + BL*Vm + dt) ...
    + bh*Vw*( lambda/( ( (q+bh*Vw + BH*Vm + dt)/r)*(r + bl*Vw + BL*Vm +dt) - q) )) ... 
    - b*( dv*Vw/p )*( (omega + ALP*( ( dv*Vw/p ) + ( dv*Vm/p ) ))/dc ) ...
    - di*( dv*Vw/p );

%%% version before substiitutions
% F(2) = Ep*(bl*Vw*Tl + bh*Vw+Th) + BL*Vm*Tl + BH*Vm*Th - (b/(1+B))*Im*C - di*Im;
F(2) = EP*(bl*Vw*( lambda + q*( lambda/( ( (q+bh*Vw + BH*Vm + dt)/r)*(r + bl*Vw + BL*Vm +dt) - q) ))/(r + bl*Vw + BL*Vm + dt) ...
    + bh*Vw+( lambda/( ( (q+bh*Vw + BH*Vm + dt)/r)*(r + bl*Vw + BL*Vm +dt) - q) ))...
    + BL*Vm*( lambda + q*( lambda/( ( (q+bh*Vw + BH*Vm + dt)/r)*(r + bl*Vw + BL*Vm +dt) - q) ))/(r + bl*Vw + BL*Vm + dt) ...
    + BH*Vm*( lambda/( ( (q+bh*Vw + BH*Vm + dt)/r)*(r + bl*Vw + BL*Vm +dt) - q) ) ...
    - (b/(1+B))*( dv*Vm/p )*( (omega + ALP*( ( dv*Vw/p ) + ( dv*Vm/p ) ))/dc ) ...
    - di*( dv*Vm/p );

end