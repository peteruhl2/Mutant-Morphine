%%% gets eigenvalues for MOE given a range of M
% close all;

npts = 500;
M = linspace(0,200,npts);
real_results = zeros(length(M),7);
imag_results = zeros(length(M),7);
max_eig = zeros(length(M),1);

for i = 1:length(M)
    MOE = MOE_func(M(i));
    E = MOE_jac(M(i),MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));
    E_imag = MOE_jac_imag(M(i),MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));
    
    real_results(i,:) = E';
    imag_results(i,:) = E_imag';
    max_eig(i) = max(real(E));
end

% hold on
% plot(M,max_eig)
% line([M(1) M(end)], [0 0], 'Linestyle','--')
% xlabel('Morphine')
% ylabel('Real part of maximum eigen value')
% title('MOE Eigenvalues')

%%% Plot real parts
% return
% close
figure()
hold on
% p2 = plot(M,results,'r-.','LineWidth',0.5);
p2 = plot(M,real_results,'b','LineWidth',0.5);
line([M(1) M(end)], [0 0], 'Linestyle','--')
% axis([0 200 -0.55 0.1])
xlabel('Morphine')
ylabel('Real part of eigenvalues')
box on
axis([0 200 -0.55 0.1])
% title('MOE Eigenvalues')

%%% Plot imaginary parts
figure()
hold on
% p2 = plot(M,results,'r-.','LineWidth',0.5);
p2 = plot(M,results,'b','LineWidth',0.5);
line([M(1) M(end)], [0 0], 'Linestyle','--')
% axis([0 200 -0.55 0.1])
xlabel('Morphine')
ylabel('Imaginary part of eigenvalues')
box on
% title('MOE Eigenvalues')
