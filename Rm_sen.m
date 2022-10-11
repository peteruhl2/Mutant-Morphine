%this is gonna be for M=0
%declare symbolic variables
clear all
syms lambda q r bl bh p b B dt dv di dc ep F omega

Tl = (lambda*(q+dt))/(dt*(q+r+dt));
Th = (lambda*r)/(dt*(q+r+dt));
C = omega/dc;

%declare Rm
mtop = -p*(B*F*Th*bh + B*F*Tl*bl - B*Th*bh - B*Tl*bl + F*Th*bh +...
    F*Tl*bl - Th*bh - Tl*bl);
mbottom = (B*di + C*b +di)*dv;
Rm = mtop./mbottom;

%calculate sensitivity indicies 
Slambda = (lambda/Rm)*diff(Rm,lambda);
Sq = (q/Rm)*diff(Rm,q);
Sr = (r/Rm)*diff(Rm,r);
Sbl = (bl/Rm)*diff(Rm,bl);
Sbh = (bl/Rm)*diff(Rm,bh);
Sp = (p/Rm)*diff(Rm,p);
Sb = (b/Rm)*diff(Rm,b);
SB = (B/Rm)*diff(Rm,B);
Sdt = (dt/Rm)*diff(Rm,dt);
Sdv = (dv/Rm)*diff(Rm,dv);
Sdi = (di/Rm)*diff(Rm,di);
Sdc = (dc/Rm)*diff(Rm,dc);
Sep = (ep/Rm)*diff(Rm,ep);
SF = (F/Rm)*diff(Rm,F);
Somega = (omega/Rm)*diff(Rm,omega);

%give values
q = 0.25; %control
r = 0.16; %control
lambda = 3690;%3690;
F = 0.2;%0.2;
bl = 1e-9;
bh = 1e-7;
p = 2500; %2500
b = 0.25;%0.005 Vitaly: 0.01 to 0.4
B = 30;%10.5;%50
dt = 0.01;
dv = 23;
di = 0.3;
dc = 0.63; %0.63

ep = 3e-5;
mu = 1;%4/24;
eta = 1;

alp = 6.7e-6;%6.7e-6
gamma = 1; %0.4;%0.4;
xi = 1;

omega = 50;

%store each senind as double
Slambda = double(subs(Slambda));
Sq = double(subs(Sq));
Sr = double(subs(Sr));
Sbl = double(subs(Sbl));
Sbh = double(subs(Sbh));
Sp = double(subs(Sp));
Sb = double(subs(Sb));
SB = double(subs(SB));
Sdt = ceil(double(subs(Sdt)));
Sdv = double(subs(Sdv));
Sdi = double(subs(Sdi));
Sdc = double(subs(Sdc));
Sep = double(subs(Sep));
SF = double(subs(SF));
Somega = double(subs(Somega));

senind = [Slambda Sq Sr Sbl Sbh Sp Sb SB Sdt Sdv Sdi Sdc Sep SF Somega];

%for label
names = ({'\lambda','q','r','\beta_l','\beta_h','p','b',...
    'B','\delta_T','\delta_V','\delta_I','\delta_C','\epsilon',...
    'F','\omega'});

%make bar graph and label
bar(senind)
set(gca,'xticklabel',names)
axis([0 16 -1.5 1.5])
xlabel('Parameter')
ylabel('Sensitivity')
title('R_0^m Sensitivity indices')