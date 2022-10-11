%%% way %%%%ing better way of getting the F-B-Mthresh contour
%%% uses fsolve instead of gd symbols
%%% just for F curve

clear all; close all

npts = 5000;
morphine = linspace(0,200,npts);
fitness = linspace(0,1,npts);
% escape = linspace(0,50,npts);
escape = 30;
results = zeros(length(escape),length(fitness));
options = optimset('Display','off');

tic
for i = 1:length(escape)
    [i]
    B = escape(i);
    parfor j = 1:length(fitness)
    
        F = fitness(j);
        
        % x = myfun(54.4860,0.1,30)
        %%% solve the equation with F and B value
        [x,fval] = fsolve(@(x) myfun(x,F,B),60,options);
        
        %%% check if positive
        if x > 0
            results(i,j) = x;
        else results(i,j) = NaN;
        end
        
    end
end
toc

hold on
plot(fitness,results,'Linewidth',2)
xlabel('Fitness cost (F)')
ylabel('Threshold morphine concentration (M_{thresh})')

%%% define equation to be solved for M with parameters F and B 
function val = myfun(M,F,B)

        %%% q-r parameters
        Mh = 100; %2.8534e-3;
        rc = 0.16;
        rm = 0.52;
        qc = 1.23e-6;
        qm = 0.25;
        n = 8;

        eta_r = @(M) (M.^n)./(Mh^n+M.^n);
        eta_q = @(M) 1-eta_r(M);

        r = rc + (rm - rc).*eta_r(M);
        q = qc + (qm - qc).*eta_q(M);

        %%% other parameters
        lambda = 3690;%3690;
        % F = 0.2;%0.2;
        bl = 1e-9;
        bh = 1e-7;
        p = 2500; %2500
        b = 0.25;%0.005 Vitaly: 0.01 to 0.4
        % B = 30;%20.5;%10.5;%50
        dt = 0.01;
        dv = 23;
        di = 0.7;
        dc = 0.2; %0.63

        ep = 3e-5;
        mu = 1;%4/24;
        eta = 1;
        EP = ep./(mu + eta*M);

        alp = 6.7e-5;%6.7e-6
        gamma = 1; %0.4;%0.4;
        xi = 1;
        ALP = alp./(gamma + xi*M);

        omega_base = 15;
        psi = 0.1;
        omega = omega_base*exp(-psi.*M);%50

        %%% got this by translating from maple
        val =  lambda * p * (M * eta - ep + mu) * dc * ...
        (((qm + dt) * bl + bh * rc) * Mh ^ (2 * n)...
        + ((qc + dt) * bl + bh * rm) ...
        * M ^ (2 * n) + ((2 * dt + qc + qm) * bl + bh * ...
        (rc + rm)) * Mh ^ n * M ^ n) / dt /...
        (Mh ^ n + M ^ n) / (M * eta + mu) ...
        / ((qm + rc + dt) * Mh ^ n + (qc + rm + dt) ...
        * M ^ n) / dv / (omega_base * exp(-(psi * M)) * b + di ...
        * dc) - (-p * lambda * (((qm + dt) * bl + bh * rc) ...
        * Mh ^ (2 * n) + ((qc + dt) * bl + bh * rm) ...
        * M ^ (2 * n) + ((2 * dt + qc + qm) * bl + bh * ...
        (rc + rm)) * Mh ^ n * M ^ n) * (1 + B) * (-1 + F) * ...
        dc / dv / dt / (Mh ^ n + M ^ n) / ((qm + rc + dt)...
        * Mh ^ n + (qc + rm + dt) * M ^ n) / ...
        (omega_base * exp(-(psi * M)) * b + di * dc * (1 + B)));
end