%%% call the coexistence virus function to solve equations 5 and 6 for Vw
%%% and Vm

clear all; close all;

global M 

M = 56;

%%% solve the two equations for Vw and Vm,
%%% but unfortunately M is global
fun = @coexist_virus_eq;
% x0 = [5e6; 1.5e6];
x0 = [7.726e5; 0];
x = fsolve(fun,x0);

Vw = x(1)
Vm = x(2)




