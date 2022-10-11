%%% want to see what happens the the MOE Vm value if M<M_thres (mutant
%%% dominates) and B goes to 0

clear all; %close all;

M = 28;
B = linspace(0,40);
Vmeq = zeros(length(B),1);

for i = 1:length(B)
    Vmeq(i) = vm_escape(M,B(i));
end

hold on
plot(B,(Vmeq))
xlabel('B (escape)')
ylabel('V_m MOE value')

Vmeq(1)