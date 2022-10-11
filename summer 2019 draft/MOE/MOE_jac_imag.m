%%% calculates jacobian eigenvalues given M and MOE values
%%% and gives imaginary parts

function E = MOE_jac_imag(M,Tl,Th,Vm,Im,C)

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
EP = ep/(mu+eta*M);

alp = 6.7e-5;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;
ALP = alp/(gamma + xi*M);

omega_base = 15;
psi = 0.1;
omega = omega_base*exp(-psi.*M);%50

BL = (1-F)*bl;
BH = (1-F)*bh;

%%% no Vw, Iw
Vw = 0;
Iw = 0;

J = zeros(7,7);

J(1,1) = -r - bl*Vw - BL*Vm - dt;
J(1,2) = q;
J(1,3) = -bl*Tl;
J(1,4) = -BL*Tl;
J(1,5) = 0;
J(1,6) = 0;
J(1,7) = 0;

J(2,1) = r;
J(2,2) = -q - bh*Vw - BH*Vm - dt;
J(2,3) = -bh*Th;
J(2,4) = -(BH*Th);
J(2,5) = 0;
J(2,6) = 0;
J(2,7) = 0;

J(3,1) = 0;
J(3,2) = 0;
J(3,3) = -dv;
J(3,4) = 0;
J(3,5) = p;
J(3,6) = 0;
J(3,7) = 0;

J(4,1) = 0;
J(4,2) = 0;
J(4,3) = 0;
J(4,4) = -dv;
J(4,5) = 0;
J(4,6) = p;
J(4,7) = 0;

J(5,1) = (1-EP)*bl*Vw;
J(5,2) = (1-EP)*bh*Vw;
J(5,3) = (1-EP)*(bl*Tl + bh*Th);
J(5,4) = 0;
J(5,5) = -b*C - di;
J(5,6) = 0;
J(5,7) = -b*Iw;

J(6,1) = EP*bl*Vw + BL*Vm;
J(6,2) = EP*bh*Vw + BH*Vm;
J(6,3) = EP*(bl*Tl + bh*Th);
J(6,4) = BL*Tl + BH*Th;
J(6,5) = 0;
J(6,6) = -(b/(1+B))*C - di;
J(6,7) = -(b/(1+B))*Im;

J(7,1) = 0;
J(7,2) = 0;
J(7,3) = 0;
J(7,4) = 0;
J(7,5) = ALP*C;
J(7,6) = ALP*C;
J(7,7) = ALP*(Iw+Im) - dc;

% E = eig(J);

%%% attemp to sort 
E = eig(J);
E = sort((imag(E)));

end