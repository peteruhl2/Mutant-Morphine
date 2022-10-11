%%% compare MOE values to coexist valuse, in particular for low M
clear all;
addpath('C:\Users\COS-DiMoLab\Desktop\mutant morphine project\summer 2019 draft\coexist_equil');
addpath('C:\Users\COS-DiMoLab\Desktop\mutant morphine project\summer 2019 draft\MOE')

M = 20;

MOE = MOE_func(M);
MOE = [MOE(1),MOE(2),0,MOE(3),0,MOE(4),MOE(5)];

coex = coexist_values(M);