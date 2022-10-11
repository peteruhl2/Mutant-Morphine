%%% gets eigenvalues for MOE given a range of M
%%% AND a single value of B !!!
% close all;

npts = 100;
M = linspace(0,200,npts);
% M = 55.5;

B = 30;
results = zeros(length(M),7); 
max_eig = zeros(length(M),1);

for i = 1:length(M)
    i
    MOE = MOE_escape_fcn(M(i),B);
    E = escape_jac_fcn(M(i),B,MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));
    
    results(i,:) = E';
    max_eig(i) = max(real(E));
end

% max_eig

hold on
plot(M,max_eig)
line([M(1) M(end)], [0 0], 'Linestyle','--')
xlabel('Morphine')
ylabel('Real part of maximum eigen value')

% return
% close
hold on
plot(M,results)
line([M(1) M(end)], [0 0], 'Linestyle','--')
axis([0 200 -0.6 0.2])
xlabel('Morphine')
ylabel('Real part of eigenvalues')
