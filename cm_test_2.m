function ydot = cm_test_2(t,y)

global lambda q r bl bh F dt p dv
global ep mu eta b di B omega alp dc M gamma xi

Tl = y(1);
Th = y(2);
Vw = y(3);
Vm = y(4);
Iw = y(5);
Im = y(6);
C = y(7);

ydot = ones(7,1);

ydot(1) = lambda + q*Th - r*Tl - bl*Vw*Tl - (1-F)*bl*Vm*Tl - dt*Tl;
ydot(2) = r*Tl - q*Th - bh*Vw*Th - (1-F)*bh*Vm*Th - dt*Th;
ydot(3) = p*Iw - dv*Vw;
ydot(4) = p*Im - dv*Vm;
ydot(5) = (1-(ep/(mu+eta*M)))*(bl*Vw*Tl + bh*Vw*Th) - b*Iw*C - di*Iw;
ydot(6) = (ep/(mu+eta*M))*(bl*Vw*Tl + bh*Vw*Th) + (1-F)*bl*Vm*Tl + ...
          (1-F)*bh*Vm*Th - (b/(1+B))*Im*C - di*Im;
ydot(7) = omega + (alp/(gamma+xi*M))*(Iw+Im)*C - dc*C;


