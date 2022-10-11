%this is gonna be for M=0
%declare symbolic variables
clear all
syms lambda q r bl bh p b B dt dv di dc ep F omega

Tl = (lambda*(q+dt))/(dt*(q+r+dt));
Th = (lambda*r)/(dt*(q+r+dt));
C = omega/dc;

%declare Rw
wtop = -p*(Th*bh*ep + Tl*bl*ep - Th*bh - Tl*bl);
wbottom = (C*b + di)*dv;
Rw = wtop./wbottom;

%calculate sensitivity indicies 
Slambda = (lambda/Rw)*diff(Rw,lambda);
Sq = (q/Rw)*diff(Rw,q);
Sr = (r/Rw)*diff(Rw,r);
Sbl = (bl/Rw)*diff(Rw,bl);
Sbh = (bl/Rw)*diff(Rw,bh);
Sp = (p/Rw)*diff(Rw,p);
Sb = (b/Rw)*diff(Rw,b);
SB = (B/Rw)*diff(Rw,B);
Sdt = (dt/Rw)*diff(Rw,dt);
Sdv = (dv/Rw)*diff(Rw,dv);
Sdi = (di/Rw)*diff(Rw,di);
Sdc = (dc/Rw)*diff(Rw,dc);
Sep = (ep/Rw)*diff(Rw,ep);
SF = (F/Rw)*diff(Rw,F);
Somega = (omega/Rw)*diff(Rw,omega);

%give values to variables
q = 0.25; %control
r = 0.16; %control
lambda = 3690;%3690;
F = 0.2;%0.2;
bl = 1e-9;
bh = 1e-7;
p = 4000; %2500
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

%store each Sx as doulbe
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
axis([0 16 -1.2 1.2])
xlabel('Parameter')
ylabel('Sensitivity')
title('R_0^w Sensitivity indices')

