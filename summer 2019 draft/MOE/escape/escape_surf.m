%%% gets eigenvalues for MOE given a range of M
%%% AND B !!!
%%% then plots the surface/contour plot
close all; clear all;

npts = 400;

M = linspace(0,100,npts);
% B = 10.^(linspace(log10(0.001),log10(30))); 
B = linspace(0,50,npts);
% results = zeros(length(M),7); 
max_eig = zeros(length(M),length(B));

tic
for i = 1:length(M)
    parfor j = 1:length(B)
        [i j]
        MOE = MOE_escape_fcn(M(i),B(j));
        E = escape_jac_fcn(M(i),B(j),MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));

        %%% reproduction numbers here
        x = escape_Rm(M(i),B(j));
        Rw = x(1);
        Rm = x(2);
    
        %%% get max eig
        max_eig(i,j) = max(real(E));
        
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

[X,Y] = meshgrid(M,B);
contourf(X,Y,max_eig')
shading interp
xlabel('Morphine')
ylabel('Escape ratio')
zlabel('Real part of maximum eigen value')
% colorbar
% title('Stability in M-B space')
