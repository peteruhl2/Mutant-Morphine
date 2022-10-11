%%% gets eigenvalues for MOE given a range of M
%%% AND a single value of F !!!
% close all; 
clear all;

M = linspace(0,200,100);
F = 0.1;
results = zeros(length(M),7); 
max_eig = zeros(length(M),1);

tic
for i = 1:length(M)
    i
    MOE = MOE_fit_func(M(i),F);
    E = fit_jac_fcn(M(i),F,MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));
    
    results(i,:) = E';
    max_eig(i) = max(real(E));
end
toc

hold on
plot(M,max_eig)
line([M(1) M(end)], [0 0], 'Linestyle','--')
xlabel('Morphine')
ylabel('Real part of maximum eigen value')

% return
% close
hold on
plot(M,results,'b')
line([M(1) M(end)], [0 0], 'Linestyle','--')
axis([M(1) M(end) -0.6 0.2])
xlabel('Morphine')
ylabel('Real part of eigenvalues')