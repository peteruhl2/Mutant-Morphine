%%% calculates MOE values for a range of M values

clear all; close all;

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

subplot(2,2,1)
plot(M1,log10(Tl1),'b','Linewidth',2)
hold on
plot(M2,log10(Tl2),'b--','Linewidth',2)
xlabel('Morphine')
ylabel('log_{10}T_l MOE value')

subplot(2,2,2)
plot(M1,log10(Th1),'b','Linewidth',2)
hold on
plot(M2,log10(Th2),'b--','Linewidth',2)
xlabel('Morphine')
ylabel('log_{10}T_h MOE value')

subplot(2,2,3)
plot(M1,log10(Im1),'b','Linewidth',2)
hold on
plot(M2,log10(Im2),'b--','Linewidth',2)
xlabel('Morphine')
ylabel('log_{10}Im MOE value')

subplot(2,2,4)
plot(M1,C1,'b','Linewidth',2)
hold on
plot(M2,C2,'b--','Linewidth',2)
xlabel('Morphine')
ylabel('C MOE value')

% suptitle('Mutant only equilibrium values')

% hold on
% subplot(5,1,1)
% plot(M,Tl)
% xlabel('M')
% ylabel('T_l MOE value')
% 
% subplot(5,1,2)
% plot(M,Th)
% xlabel('M')
% ylabel('T_h MOE value')
% 
% subplot(5,1,3)
% plot(M,log10(Vm))
% xlabel('M')
% ylabel('log_10 Vm MOE value')
% 
% subplot(5,1,4)
% plot(M,Im)
% xlabel('M')
% ylabel('I_m MOE value')
% 
% subplot(5,1,5)
% plot(M,C)
% xlabel('M')
% ylabel('C MOE value')













