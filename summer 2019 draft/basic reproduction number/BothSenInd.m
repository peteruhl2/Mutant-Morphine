
senind_w = [1,-0.579862941486558,0.603057459146020,0.0159901599015990,0.00984009840098401,1,-0.960375193902936,0,-1,-1,-0.0396248060970640,0.960375193902936,-1.50002250033751e-05,0,-0.960375193902936,7.50011250168753e-06,7.50011250168753e-06,0.0960375193902936];
senind_m = [1,-0.579862941486558,0.603057459146020,0.0159901599015990,0.00984009840098401,1,-0.438779095976966,0.424624931590613,-1,-1,-0.561220904023034,0.438779095976966,0,-0.111111111111111,-0.438779095976966,0,0,0.0438779095976966];

senind = [senind_w; senind_m]';

%for label
names = ({'\lambda','q','r','\beta_l','\beta_h','p','b',...
    'B','\delta_T','\delta_V','\delta_I','\delta_C','\epsilon',...
    'F','\omega','\mu','\eta','\psi'});

%make bar graph and label
figure()
bar(senind)
xticks([1:18])
set(gca,'xticklabel',names)
set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)
axis([0 19 -1.2 1.2])
% xlabel('Parameter')
% ylabel('Sensitivity')
xtickangle(0)
% title('R_0^w Sensitivity indices')
legend('R_0^w','R_0^m')
title('a)                                                                                                                         ','Fontsize',14)
