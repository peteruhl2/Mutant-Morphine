%%% shows what the MOE value of Vm is for a range of M

clear all; close all;

M1 = 0:55;
results1 = zeros(length(M1),1);

for i=1:length(M1)
    results1(i) = vm_solver_2(M1(i));
end

M2 = 56:200;
results2 = zeros(length(M2),1);

for i=1:length(M2)
    results2(i) = vm_solver_2(M2(i));
end

% figure
hold on
% plot(M,log10(results),'*','MarkerSize',2)
plot(M1,log10(results1),'b','Linewidth',2)
plot(M2,log10(results2),'b--','Linewidth',2)
xlabel('Morphine')
ylabel('V_m MOE value')
% title('V_m MOE values')
box on