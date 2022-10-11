%%% get eigen values of jacobian for coexist eq
%%% given a range of M
close all;

npts = 5000;
morphine = linspace(0,200,npts);
real_results = zeros(length(morphine),7);
imag_results = zeros(length(morphine),7);
max_eig = zeros(length(morphine),1);

tic
for i = 1:length(morphine)
    i
    M = morphine(i);
    coex = coexist_values(M);
    E = coexist_jac(M,coex(1),coex(2),coex(3),coex(4),coex(5),coex(6),coex(7));
    E_imag = coexist_jac_imag(M,coex(1),coex(2),coex(3),coex(4),coex(5),coex(6),coex(7));
    
    real_results(i,:) = E';
    imag_results(i,:) = E_imag';
    max_eig(i) = max(real(E));
end
toc

% hold on
% plot(morphine,max_eig)
% line([morphine(1) morphine(end)], [0 0], 'Linestyle','--')
% xlabel('Morphine')
% ylabel('Real part of maximum eigen value')
% title('Coex Eigenvalues')

% return
% close
hold on; box on
p = plot(morphine,real_results,'b','LineWidth',1);
line([morphine(1) morphine(end)], [0 0], 'Linestyle','--')
axis([0 200 -0.6 0.2])
xlabel('Morphine')
ylabel('Real part of eigenvalues')
title('Coexist. Eigenvalues')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOE eigenvalues here
%%% gets eigenvalues for MOE given a range of M
% close all;
addpath('C:\Users\COS-DiMoLab\Desktop\mutant morphine project\summer 2019 draft\basic reproduction number')

% npts = 100;
M = linspace(0,200,npts);
real_results_MOE = zeros(length(M),7);
imag_results_MOE = zeros(length(M),7);
max_eig = zeros(length(M),1);

for i = 1:length(M)
    i
    MOE = MOE_func(M(i));
    E = MOE_jac(M(i),MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));
    E_imag = MOE_jac_imag(M(i),MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));
    
    real_results_MOE(i,:) = E';
    imag_results_MOE(i,:) = E_imag';
    max_eig(i) = max(real(E));
end

% hold on
% plot(M,max_eig)
% line([M(1) M(end)], [0 0], 'Linestyle','--')
% xlabel('Morphine')
% ylabel('Real part of maximum eigen value')
% title('MOE Eigenvalues')

% return
% close
% figure()
hold on
p2 = plot(M,real_results_MOE,'r-.','LineWidth',1)
line([M(1) M(end)], [0 0], 'Linestyle','--')
axis([M(1) M(end) -0.55 0.1])
xlabel('Morphine')
ylabel('Real part of eigenvalues')
% title('MOE Eigenvalues')
legend([p(1) p2(1)],{'Coexist.','MOE'},'Fontsize',14)


%%% Plot imagiary parts ===================================================
figure()
hold on; box on
p = plot(morphine,imag_results,'b','LineWidth',1);
line([morphine(1) morphine(end)], [0 0], 'Linestyle','--')
axis([M(1) M(end) -0.6 0.2])
xlabel('Morphine')
ylabel('Imaginary part of eigenvalues')
title('Coexist. Eigenvalues')

hold on
p2 = plot(M,imag_results_MOE,'r-.','LineWidth',1)
line([M(1) M(end)], [0 0], 'Linestyle','--')
axis([0 200 -0.2 0.2])
xlabel('Morphine')
ylabel('Real part of eigenvalues')
% title('MOE Eigenvalues')
legend([p(1) p2(1)],{'Coexist.','MOE'},'Fontsize',14)
