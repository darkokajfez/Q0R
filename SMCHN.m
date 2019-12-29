%**********************************SMCHN.m************************
%   plotting the normalized Smith chart as figure(nfig)
%   DON'T FORGET TO WRITE hold off WHEN FINIHED WITH YOUR PLOT
%   Copyright (c) D. Kajfez 2010

function U=SMCHN(nfig)

%   compute the detailed Smith chart
phcr=linspace(0,2*pi,100);
cir=exp(1i*phcr);
hor=zeros(10,1);
ver=linspace(-1,1,10);
phiA=linspace(0,2*pi,73);
xA=0.5+0.5*cos(phiA);
yA=0.5*sin(phiA);
phiB=linspace(pi,3*pi/2,19);
xB=1+cos(phiB);
yB=1+sin(phiB);
phiCir=linspace(pi/2,pi,19);
xCir=1+cos(phiCir);
yCir=-1+sin(phiCir);
xcir05=2*cos(phcr)/3+1/3;
ycir05=2*sin(phcr)/3;
xcir2=cos(phcr)/3+2/3;
ycir2=sin(phcr)/3;
xcir02=1/6-5*cos(phcr)/6;
ycir02=5*sin(phcr)/6;
xcir5=5/6-cos(phcr)/6;
ycir5=sin(phcr)/6;
phi05=linspace(0,2*atan(0.5),50);
x05=1-2*sin(phi05);
y05=2-2*cos(phi05);
phi2=linspace(0,2*atan(2),50);
x2=1-0.5*sin(phi2);
y2=0.5-0.5*cos(phi2);
phi02=linspace(0,2*atan(0.2),50);
x02=1-5*sin(phi02);
y02=5-5*cos(phi02);
phi5=linspace(0,2*atan(5),50);
x5=1-0.2*sin(phi5);
y5=0.2-0.2*cos(phi5);
%   compute the position of text notation
alphm=2*atan(0.5);
xm=1-2*sin(alphm)+0.08;  ym=2-2*cos(alphm);
xm2=0.05;   ym2=0.93;
alphm3=2*atan(2);
xm3=1-0.5*sin(alphm3)+0.05; ym3=0.5-0.5*cos(alphm3);
alphm4=2*atan(0.2);
xm4=1-5*sin(alphm4)+0.05; ym4=5-5*cos(alphm4);
alphm5=2*atan(5);
xm5=1-0.2*sin(alphm5)-0.1; ym5=0.2-0.2*cos(alphm5); 
figure(nfig)
%   add the plots to the detailed Smith chart
plot(real(cir),imag(cir),'k',ver,hor,'k',xA,yA,'k',xB,yB,'k',...
    xCir,yCir,'k',xcir05,ycir05,'k',xcir2,ycir2,'k',x05,y05,'k');
hold on;
plot( x05,-y05,'k',x2,y2,'k',x2,-y2,'k',x02,y02,'k',x5,y5,'k',...
    x02,-y02,'k',x5,-y5,'k',xcir02,ycir02,'k',xcir5,ycir5,'k');
csp=[0.6 0.6 0.6];
text(xm,ym,'j0.5','Color',csp);text(xm2,ym2,'j1.0','Color',csp);
text(xm3,ym3,'j2.0','Color',csp);
text(xm,-ym,'-j0.5','Color',csp);text(xm2,-ym2,'-j1.0','Color',csp);
text(xm3,-ym3,'-j2.0','Color',csp);
text(xm4,ym4,'j0.2','Color',csp);    text(xm4,-ym4,'-j0.2','Color',csp);
text(xm5,ym5,'j5.0','Color',csp);   text(xm5,-ym5,'-j5.0','Color',csp);   
text(0,-0.08,0,'1.0','Color',csp);   
text(-0.33,-0.08,'0.5','Color',csp);text(0.33,-0.08,'2.0','Color',csp);
text(-2/3, -0.08,'0.2','Color', csp);text(2/3, -0.08,'5.0','Color', csp);
%text(-1.1,0.95,'(c) 2010 D. Kajfez','Color',csp);
%axis([-1.001  1.001 -1.001 1.001]);axis('square');axis off;
axis([-1.01  1.01 -1.01 1.01]);axis('square');axis off;



