%%% gets eigenvalues for MOE given a range of M
%%% AND F !!!
%%% then plots the surface/contour plot
close all; clear all;

npts = 100; %400

M = linspace(0,100,npts);
% B = 10.^(linspace(log10(0.001),log10(30))); 
F = linspace(0,1,npts);
% results = zeros(length(M),7); 
max_eig = zeros(length(M),length(F));

tic
for i = 1:length(M)
    parfor j = 1:length(F)
        [i j]
        MOE = MOE_fit_func(M(i),F(j));
        E = fit_jac_fcn(M(i),F(j),MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));
        
        %%% reproduction numbers here
        x = fit_Rm(M(i),F(j));
        Rw = x(1);
        Rm = x(2);
        
        %%% get max eig
        max_eig(i,j) = max(E);
        
        %%% check if IFE
        if (Rm < 1) && (Rw < 1)
            max_eig(i,j) = -1;
        end
       
        %%% check if coexist
        if (Rw > 1) && (Rw > Rm)
            max_eig(i,j) = 1;
        end
              
    end
end
toc 
        

[X,Y] = meshgrid(M,F);
contourf(X,Y,max_eig')
% surf(X,Y,max_eig)
shading interp
xlabel('Morphine')
ylabel('Fitness cost')
zlabel('Real part of maximum eigen value')
% colorbar
% title('Stability in M-F space')