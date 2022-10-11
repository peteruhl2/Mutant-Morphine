function yp = mut_model(t,y)

global lambda q r bl bh F dt p dv
global b di B omega dc EP ALP 

Tl = y(1);
Th = y(2);
Vw = y(3);
Vm = y(4);
Iw = y(5);
Im = y(6);
C = y(7);

yp = ones(7,1);

yp(1) = lambda + q*Th - r*Tl - bl*Vw*Tl - (1-F)*bl*Vm*Tl - dt*Tl;
yp(2) = r*Tl - q*Th - bh*Vw*Th - (1-F)*bh*Vm*Th - dt*Th;
yp(3) = p*Iw - dv*Vw;
yp(4) = p*Im - dv*Vm;
yp(5) = (1-EP)*(bl*Vw*Tl + bh*Vw*Th) - b*Iw*C - di*Iw;
yp(6) = EP*(bl*Vw*Tl + bh*Vw*Th) + (1-F)*bl*Vm*Tl + ...
          (1-F)*bh*Vm*Th - (b/(1+B))*Im*C - di*Im;
yp(7) = omega + ALP*(Iw+Im)*C - dc*C;

end