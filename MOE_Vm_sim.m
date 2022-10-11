clear all

% M = 0:0.00001:0.01;
M = 0:300;
Vm = zeros(length(M),1);

for i = 1:length(M)
    Vm(i) = n_MOE_Vm_solver_v2(M(i));
end

plot(M,Vm)
hold on
xlabel('M')
ylabel('Vm')