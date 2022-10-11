%%% coexistence equilibrium, get Vw and Vm values for a range of M
clear all

npts = 100;

global M 
morphine1 = linspace(0,55,npts);
Vw1 = zeros(1,length(morphine1));
Vm1 = zeros(1,length(morphine1));

tic
for i = 1:length(morphine1)
    M = morphine1(i);
    x = coexist_sim(M);
    Vw1(i) = x(1);
    Vm1(i) = x(2);
end

morphine2 = linspace(55,200,npts);
Vw2 = zeros(1,length(morphine2));
Vm2 = zeros(1,length(morphine2));

for i = 1:length(morphine2)
    M = morphine2(i);
    x = coexist_sim(M);
    Vw2(i) = x(1);
    Vm2(i) = x(2);
end
toc

% figure
hold on
% plot(morphine1,log10(Vw1),'b')
plot(morphine1,log10(Vm1),'r')
plot(morphine2,log10(Vw2),'b')
plot(morphine2,log10(Vm2),'r')
xlabel('Morphine')
ylabel('log10 Viral population densities')
% legend('V_w','V_m')
axis([0 200 -0.5 6])
title('Equilibrium viral loads')