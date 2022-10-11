clear all

J = zeros(7,7);
Em = 0:5:300;

MOE = zeros(7,length(Em));
A = zeros(length(Em),1);

for k=1:length(Em)
    M = Em(k);
    
    Mh = 100;
    rc = 0.16;
    rm = 0.52;
    qc = 1.23e-6;
    qm = 0.25;
    n = 7.8731;

    eta_r = @(M) (M^n)/(Mh^n+M^n);
    eta_q = @(M) 1-eta_r(M);

    r = rc + (rm-rc)*eta_r(M);
    q = qc + (qm - qc)*eta_q(M);

    lambda = 3690;%3690;
    F = 0.2;%0.2;
    bl = 1e-9;
    bh = 1e-7;
    p = 2500; %2500
    b = 0.25;%0.005 Vitaly: 0.01 to 0.4
    B = 25;%10.5;%50
    dt = 0.01;
    dv = 23;
    di = 0.3;
    dc = 0.63; %0.63

    ep = 3e-5;
    mu = 1;%4/24;
    eta = 1;
    EP = ep/(mu+eta*M);

    alp = 6.7e-6;%6.7e-6
    gamma = 1; %0.4;%0.4;
    xi = 1;
    ALP = alp/(gamma+xi*M);

    omega = 50*exp(-0.05*M);%50

    BL = (1-F)*bl;
    BH = (1-F)*bh;

    Vm = n_MOE_Vm_solver_v2(M);
%     Vm = 3.934691879429678e+05;
    Im = dv*Vm/p;
    C = omega/(dc-ALP*Im);
    Th = r*lambda/((q+BH*Vm+dt)*(r+BL*Vm+dt)-r*q);
    Tl = Th*(q+BH*Vm+dt)/r;
    Vw = 0;
    Iw = 0;
        
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

    E = eig(J);

    A(k) = max(real(E));

    MOE(:,k) = E; 

end

% plot(Em,A,'*')
plot(Em,A)
hold on

xlabel('M')
ylabel('Max(Re(\sigma(J)))')
title('Morphine effect on MOE stability')

plot([0 M(end)],[0 0],'--')