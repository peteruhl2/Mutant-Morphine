%%% gets eigenvalues for MOE given a range of M
% close all;

npts = 500;
M = linspace(0,200,npts);
results = zeros(length(M),7); 
max_eig = zeros(length(M),1);

for i = 1:length(M)
    i
    MOE = MOE_func(M(i));
    E = MOE_jac_imag(M(i),MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));
    
    results(i,:) = E';
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
hold on
% p2 = plot(M,results,'r-.','LineWidth',0.5);
p2 = plot(M,results,'b','LineWidth',0.5);
line([M(1) M(end)], [0 0], 'Linestyle','--')
% axis([0 200 -0.55 0.1])
xlabel('Morphine')
ylabel('Imaginary part of eigenvalues')
box on
% title('MOE Eigenvalues')