%%% uses eigenshuffle to plot the eigenvalues of J for the MOE
%%% sets range of M, calls MOE_func for MOE values,
%%% calls MOE_jac_shuf to get array of J's for eigenshuffle
clear all; close all

npts = 1000;
M = linspace(0,200,npts);
J = zeros(7,7,length(M));

for i = 1:length(M)
    MOE = MOE_func(M(i));
    J(:,:,i) = MOE_jac_shuf(M(i),MOE(1),MOE(2),MOE(3),MOE(4),MOE(5));
end

[Vseq,Dseq] = eigenshuffle(J);

hold on
line([M(1) M(end)], [0 0], 'Linestyle','--')
plot(M,Dseq(1,:),'b'); %good
plot(M,Dseq(2,:),'b'); %good, 1 and 2 are complex conjugates
plot(M,Dseq(3,:),'b'); %weird one
plot(M,Dseq(4,:),'b'); %good
plot(M,Dseq(5,:),'b'); %weird one
plot(M,Dseq(6,:),'b'); %good, large negative
plot(M,Dseq(7,:),'b'); %good, large negative
axis([0 200 -0.8 0.2])
xlabel('Morphine')
ylabel('Real part of eigenvalues')