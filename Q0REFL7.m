%********************************* Q0REFL7.m ******************
%   Reflection method of high-Q measurement.
%   Includes the stretching option.
%   Includes  the a posteriori error estimates
%   Includes the line length estimate
%   Needed functions: SMCH.m, SMCHN.m, CLIP.m, kfita.m
%   Dialog prompts, manual pruning
%   Copyright (c) D. Kajfez 2019 

clear;
clear global;
close all;

global   n  x ndeg TL

%ndeg=120;
%ndeg=180;
ndeg=90;

prompt={'Input file name','stretch: y or n?','TL= 0 or 1?','Comment: y or n?'};

dlg_title='Input ';

num_lines=1;

deff={'FD-40.mea','y','1','y'};    
%deff={'v14n1.txt','n','1','n'};
%deff={'v14n0p1.txt','n','1','n'};
%deff={'v14.txt','n','0','n'};
%deff={'v14k50.txt','n','1','n'};
%deff={'v14k80.txt','n','1','n'};
%deff={'v13n0p5.txt','n','0','n'};
%deff={'v13k20.txt','n','1','n'};
%deff={'v13n0p1.txt','n','0','n'};
%deff={'v13.txt','n','0','n'};
%deff={'v14k50n1.txt','n','1','n'};
    












answer=inputdlg(prompt,dlg_title,num_lines,deff);
%   if Cancel is pressed
K=isempty(answer);
if K == 1
    return
end

%   file name
iname=answer{1};

%   Does it need stretching?
need=answer{2};
if need == 'n'
    stre=0;
elseif need == 'y'
    stre=1;
end

%   line length exists?
TL=str2num(answer{3});

%   comment needed or not
cmnt=answer{4};

disp([iname,'  stretch= ',need,' TL= ',answer{3},' comment= ',cmnt])  

%
%   decide whether header exists
%
fid=fopen(iname,'r');
tin=fscanf(fid,'%c');
chr1=tin(1);
if chr1=='!'
    method='header';
else
    method='no header';
end

switch lower(method);
    case('header')
[nsz1 nsz2]=size(tin);
passw=char(35);
for jj=1:500
    word1 = tin(jj);
    if word1 == passw
        jpound=jj;
        break
    end
end
for ii=jpound:jpound+40
            word2= tin(ii);
            if word2 == char(10)
                iline=ii;
                break
            end
end

cleantin=tin(iline+1:nsz2);
nin=str2num(cleantin);  %data matrix

freq=nin(:,1);

otherwise
    nin=str2num(tin);
%    nin=str2double(tin);
    freq=nin(:,1);
end

rawgr=nin(:,2);
rawgi=nin(:,3);
nraw=length(freq);
rawG=rawgr + 1i*rawgi;

refl=abs(rawG);
reflmin=min(refl);
S11=rawG;
fre=freq;

%   relative increment
rawglo=rawG(1:nraw-1);
rawghi=rawG(2:nraw);
distg=rawghi-rawglo;
incrg=abs(distg);
ingmax=max(abs(incrg));
distrel=incrg/ingmax;
delf=fre(2)-fre(1);
frin=fre(1:nraw-1)+delf/2;
subG=rawG;
subgr=real(subG);
subgi=imag(subG);
subref=refl;
hor=[fre(1) fre(nraw)];
ver=[0.5 0.5];

%   plot the entire range of frequency data
%   prune the data by pointing the mouse at the low and high end points
%   in Fig 1, then Fig. 1 will appear again, with pruned portions in red color
figure('renderer','painters','resize','off','visible','off')

figure(1);
axes('fontsize',14)
plot(fre,refl,'bs',frin,distrel,'g-',hor,ver,'--g','linewidth',2);
title(['Q0REFL7.m  ',iname,'  ', num2str(date)],'Fontsize',14);
xlabel('frequency','Fontsize',14);
ylabel('reflection coefficient magnitude','Fontsize',14);
axis([fre(1) fre(nraw) 0.0 1.05]);

if stre==1
    uiwait(msgbox('Begin stretching frequency!'))
    uiwait(msgbox('Choose a low frequency point.','Stretching','modal'))
    [x1 y1]=ginput(1);
    uiwait(msgbox('Choose a high frequency point.','Stretching','modal'))
    [x2 y2]=ginput(1);
    [snp1 sisnp1]=min(abs(fre-x1));
    [snp2 sisnp2]=min(abs(fre-x2));
    fstrelo=fre(1:sisnp1);  fstrehi=fre(sisnp2:nraw); 
    fstremid=fre(sisnp1+1:sisnp2-1);   sreflmid=subref(sisnp1+1:sisnp2-1);    
    fstre1=fre(sisnp1+1);
    nstre=(sisnp2-sisnp1-1);
    fstre2=fre(sisnp2-1);


    % the frequency scale will be stretched
    % Erase the figure and plot plot the streched range for pruning 
	clf
    figure('renderer','painters','resize','off','visible','off')

    figure(1);

    axes('fontsize',14)
    hor=[fstre1,fstre2];
    ver=[0.5 0.5];
    axis([fstre1 fstre2 0.0 1.05]);

    frelo=fstremid(1:nstre-1);    frehi=fstremid(2:nstre);
    frin=0.5*(frelo+frehi);
    distin=distrel(sisnp1+1:sisnp2-2);

    plot(fstremid,sreflmid,'bo',frin,distin,'g-',hor,ver,'--g','linewidth',2);

    xlabel('frequency ');
    ylabel('|S11|');
	sfre=fre(sisnp1+1:sisnp2-1);
	nstre=length(sfre);
	ibeg=sisnp1+1;  iend=sisnp2-1;

    uiwait(msgbox('Begin pruning'))
    uiwait(msgbox('Choose a low frequency point.','Pruning','modal'))
    [x1, y1]=ginput(1);
    uiwait(msgbox('Choose a high frequency point.','Pruning','modal'))
    [x2, y2]=ginput(1);
    [snp1, isnp1]=min(abs(fre-x1));
    [snp2, isnp2]=min(abs(fre-x2));
    grprun1=subgr(sisnp1:isnp1);    giprun1=subgr(sisnp1:isnp1);  reflprun1=subref(sisnp1:isnp1);
    grprun2=subgr(isnp2:sisnp2);    giprun2=subgi(isnp2:sisnp2);    reflprun2=subref(isnp2:sisnp2);
    fprun1=fre(sisnp1:isnp1);  fprun2=fre(isnp2:sisnp2);                               
    fmid=fre(isnp1+1:isnp2-1);  reflmid=refl(isnp1+1:isnp2-1);    
    npruned=isnp2-isnp1-1;
    sf1=fmid(1);    sf2=fmid(npruned);

    %Erase and plot the stretched  range showing pruned parts in red
    clf
figure('renderer','painters','resize','off','visible','off')
    
    figure(1)
    axes('fontsize',14)
    axis([sf1 sf2 0.0 1.05]);

	plot(fprun1,reflprun1,'r.',fmid,reflmid,'bo',fprun2,reflprun2,'r.');
	title(['Q0REFL7.m  ',iname,'  ',' stretched from ',num2str(nraw),...
        ' to ',num2str(nstre) ' points']);
	xlabel('frequency ');
    ylabel('|S11|');
    grid;
else

    %Erase, showing pruned parts in red

    uiwait(msgbox('Begin pruning'))
    uiwait(msgbox('Choose a low frequency point.','Pruning','modal'))
    [x1, y1]=ginput(1);
    uiwait(msgbox('Choose a high frequency point.','Pruning','modal'))
    [x2, y2]=ginput(1);
    [snp1, isnp1]=min(abs(fre-x1));
    [snp2, isnp2]=min(abs(fre-x2));
    grprun1=subgr(1:isnp1);    giprun1=subgi(1:isnp1);  reflprun1=subref(1:isnp1);
    grprun2=subgr(isnp2:nraw);    giprun2=subgi(isnp2:nraw);    reflprun2=subref(isnp2:nraw);  
    gr=subgr((isnp1+1):(isnp2-1));  gi=subgi((isnp1+1):(isnp2-1));
    fprun1=fre(1:isnp1);  fprun2=fre(isnp2:nraw); 
    fmid=freq(isnp1+1:isnp2-1);refl=subref(isnp1+1:isnp2-1);
    f1=fmid(1);
    npruned=fix(isnp2-isnp1-1);
    f2=fmid(npruned);

    clf
    rawGlow=rawG(1:isnp1);
    rawGhi=rawG(isnp2:nraw);


    figure('renderer','painters','resize','off','visible','off')
    figure(1);
    axes('fontsize',14)
    
    
    plot(fprun1,reflprun1,'r.',fmid,refl,'bo',fprun2,reflprun2,'r.','linewidth',1);
    title(['Q0REFL7.m  ',iname,'  ',' pruned from= ',num2str(nraw),...
        '   to= ',num2str(npruned),' points  '],'fontsize',14);
    xlabel('frequency ','fontsize',14);
    ylabel('reflection coefficient magnitude','fontsize',14);
    axis([fre(1), fre(nraw), 0, 1.05]);
    grid;
end

%   pruned S parameters
    PS11=S11(isnp1+1:isnp2-1);
    PS11r=real(PS11);    PS11i=imag(PS11);
    fre=freq(isnp1+1:isnp2-1);

%return
n=length(fre);
%   
%   Pruned data fit
%
%   normalization frequency f00
%
AS11=abs(PS11);
[~, idia]=min(AS11);   
f00=fre(idia);   
f0=f00;     %value at the start of iteration
  ena=ones(n,1);
  eps=ena;
  p=zeros(n,1);
  gam1=PS11;
for ite=1:5
  x=2*(fre-f0)/f0;
  e3=gam1.*x;
  e1=-x;
  e2=-ena;
  F=-gam1;
  x2=abs(x).^2;
  ga2=abs(gam1).^2;
  if ite == 1
      p=1./(x2.*(ena+ga2)+ena);
  else
      p=sig(2)^2./(x2*sig(1)^2+sig(2)^2+x2.*ga2*sig(3)^2);
  end
  E=[e1 e2 e3];
  PE=[(p.*e1) (p.*e2) (p.*e3)];
  C=E'*PE;
  q1=E'*(p.*F);
  D=inv(C);
  g=D*q1;
  eps=g(1)*e1+g(2)*e2+g(3)*e3-F;
  S1sq=eps'*(p.*eps);
  Fsq=F'*(p.*F);
  sumden=C(1,1)*D(1,1)+C(2,2)*D(2,2)+C(3,3)*D(3,3);
  for m=1:3
      sig(m)=sqrt(abs(D(m,m)*S1sq/sumden));
  end
  diam1=2*abs(g(2)*g(3)-g(1))/abs(conj(g(3))-g(3));
  gamc1=(conj(g(3))*g(2)-g(1))/(conj(g(3))-g(3));
  gamd1=g(1)/g(3);
  gam1L2=2*gamc1-gamd1;
  delt=(g(2)-gam1L2)/(gam1L2*g(3)-g(1));
  f0=f0*(1+0.5*real(delt));
end
ft=x;   % names could have been replaced above
a=g;    % names could have been replaced above
%
%   Unreliable input data
%
unrel=0;
DAS=sqrt(abs(S1sq/Fsq));
if DAS > 0.1
    disp('UNRELIABLE INPUT DATA!')
    unrel=1;
end

f011=f0;    %corrected normalization frequency, after iteration
dia=2*abs(a(2)*a(3)-a(1))/abs(conj(a(3))-a(3)); %called d in QFMUM           
ft=2*(fre/f011-1);
gam1com=(a(1)*ft+a(2))./(a(3)*ft+1);    %analytic form of Q-circle
gam1d=a(1)/a(3);                        %detuned reflection coeff.
zs=(1+gam1d)/(1-gam1d);
gam1c=(conj(a(3))*a(2)-a(1))/(conj(a(3))-a(3)); %center of Q-circle
QL1=imag(a(3)); %not accurate!

%   check whether the points are equidistant from the center
equi=abs(gam1-gamc1);
avequi=mean(equi);
sdequi=std(equi);
rough=sdequi/avequi;
if rough > 0.02
    disp(['UNRELIABLE INPUT DATA! ',num2str(rough)])
    forced=1;
    rs=0;
    unrel=2;
end

%   Lossless or lossy decision
%
agamd=abs(gamd1);
if (0.999 < agamd) && (agamd < 1.05)
    ds=2;
    rs=0;
    forced=1;
%    gam1d=(1i*zs-1)/(1i*xs+1);
    disp('AUTOMATIC FORCED LOSSLESS, |Gamma sub d| = 1')

elseif agamd < 0.999   
%   Lossy coupling
    ang1=-atan2(imag(-gam1d),real(-gam1d));
    ang2=atan2(imag(gam1c-gam1d),real(gam1c-gam1d));
    angtot=CLIP(ang1+ang2);
    cob=cos(angtot);
    ds=(1-abs(gam1d)^2)/(1-abs(gam1d)*cob);      
    rs = 2/ds-1;   % disp([' rs= ',num2str(rs)])

elseif agamd > 1.05
    disp('BAD DATA, |Gamma sub d| > 1.05, STOP')
    return
        
end   %end of lossless or lossy decision

z1=(1+gam1)./(1-gam1);
    zstem=(1+gam1d)/(1-gam1d);
    xs=imag(zstem);
    zs=rs+1i*xs;
    zcorr=z1-zs;
if forced == 0
    kapa=dia/(ds-dia);    %total coupling
    
%   find elements of the equivalent circuit

    gam=-2*atan(xs/(1+rs));
    ejgam=exp(1i*gam);
    dejgam=a(2)-gam1d;  
    g0tem=(1+rs)/((kapa)*((1+rs)^2+xs^2));
    r0tem=1/g0tem;
    %disp(['g0tem= ',num2str(g0tem,6)])
    QL1=imag(a(3));
    Q01=QL1*(1+kapa);
end

if forced == 1
    QL1=imag(a(3));
    kapa=dia/(2-dia);
    Q01=QL1*(1+kapa);
    rs=0;
end



%*********************determining the line  length************************

 if   TL==0
    thedeg=linspace(0,ndeg-1,ndeg);
    thedeg=thedeg';
%    disp(['size(thedeg)= ',num2str(size(thedeg))])
    delgam=zeros(ndeg,1);
     xs=imag(zs);
    rs=real(zs);
    ej=1;
 end
 
 if  TL==1
%   find the first minimum of delgam
    [delgam, ndmin] = kfita(gam1);
%    disp(['size(delgam)= ',num2str(size(delgam))])     
      
%     for nt=1:ndmin+3  
%        thedeg=(nt-2);            
%        therad=thedeg*(pi/180);
%        ej=exp(+1i*(2*therad));  %moving toward load
%     end
therad=ndmin*pi/180;    ej=exp(+1i*2*therad);
      
    delgam=delgam(1:ndmin+3);
    thedeg=linspace(-1,ndmin+1,ndmin+3);
    thedeg=thedeg';
%    xs=imag(zs);
%    rs=real(zs);
 end
     

figure('renderer','painters','resize','off','visible','off')
    figure(2)
    plot(thedeg,delgam,'b-')
    xlabel('theta (deg)')
    ylabel('delgam')
    title(['Q0REFL7.m ',iname])
    grid

%   rotate toward load
gam1=gam1*ej;
gam1d=gam1d*ej;
gam1c=gam1c*ej;
gam1com=gam1com*ej;
gam1L2=gam1L2*ej;

%   Find rs, xs, zcorr, kapa
%forced=0;
    ang1=-atan2(imag(-gam1d),real(-gam1d));
    ang2=atan2(imag(gam1c-gam1d),real(gam1c-gam1d));
    angtot=CLIP(ang1+ang2);
    cob=cos(angtot);
    ds=(1-abs(gam1d)^2)/(1-abs(gam1d)*cob);
    if forced == 0
        rs = 2/ds-1;   % disp([' rs= ',num2str(rs)])
    end
    z1=(1+gam1)./(1-gam1);

    zstem=(1+gam1d)/(1-gam1d);
    xs=imag(zstem);
    zs=rs+1i*xs;
    zcorr=z1-zs;
    kapa=dia/(ds-dia);    %total coupling   
    
%   find fL from analytic expression
gam1L=2*gam1c-gam1d;
f1L=(a(2)-gam1L2)/(gam1L2*a(3)-a(1));
f1L=real(f1L);
f1byfn=(1+f1L/2);
fL=f011*f1byfn;
fLG=1e-09*fL;
ftL=2*(fre/fL-1);

%   find the corrected unloaded resonant frequency
gam0raw=(zcorr-1)./(zcorr+1);
ig=imag(gam0raw);
nprod=ig(1)*ig(n);
if nprod < 0
    disp('Pruned range includes f0')
else
    disp('Choose a wider pruned range, because f0 is outside')
    return
end
[i1,mi]=min(abs(ig));
if ig(mi) < 0
    nplus=mi;
    nminus=mi-1;
elseif ig(mi) > 0
    nminus=mi;
    nplus=mi+1;
end
%   interpolate between two points nearest to the horizontal axis
f0corr=fre(nminus)+ig(nminus)*(fre(nplus)...
    -fre(nminus))/(ig(nminus)+abs(ig(nplus)));


%********************************begin a posteriori**********************       
ft0=2*(fre/f0corr-1); 
%   find new Moebius coeficients a0
e1=ft0;
e2=ena;
e3=-gam0raw.*ft0;
E=[e1 e2 e3];
%h=gam0raw;
M=E'*E;
H=E'*gam0raw;
a0=M\H;
gam0com=(a0(1)*ft0+a0(2))./(a0(3)*ft0+1);
gam0c=(conj(a0(3))*a0(2)-a0(1))/(conj(a0(3))-a0(3));
dia0=2*abs(a0(2)*a0(3)-a0(1))/abs(conj(a0(3))-a0(3)); 
gam0d=a0(1)/a0(3);
gam02=2*gam0c-gam0d;
dia1in=[gam0d gam02];
Q0=(2/(1-real(a0(2))))*sqrt(abs(a0(1)*a0(3)));
%disp(['Q0= ',num2str(Q0new,7)])
%g0=(1-a0(2))/(1+a0(2));
g0=(1+rs)/(kapa*((1+rs)^2+xs^2));  %11-7-2019
g0=real(g0);    r0=1/g0;

%   find eQ0
apgam0d=-1;
[fr0,n0]=min(abs(fre-f0corr));
%[frmin,nmin]=min(abs(fre-f00));
if n0 == 1
    disp('pruning start must be lowered')
    return
end
if n0 == n
    disp('pruning end should be increased')
    return
end
ft0sh=ft0(2:n-1);
gam0rawsh=gam0raw(2:n-1);
%g00= (1+rs)/(kapa*((1+rs)^2+xs^2));  %11-7-2019
temp1=(2/g0)./((gam0rawsh+1).^2);
temp2=(a0(1)-a0(2)*a0(3))./((a0(3)*ft0sh+1).^2);
eQ0=-imag(temp1.*temp2);
Q0mean=mean(eQ0);
Q0std=std(eQ0);
gam1f0=gam1com(n0);


figure('renderer','painters','resize','off','visible','off')
SMCHN(3)
hold on;
plot(real(gam0raw),imag(gam0raw),'bo',real(gam0com),imag(gam0com),'r-',...
    real(gam0d),imag(gam0d),'ro',real(gam0c),imag(gam0c),'ro',...
    real(gam02),imag(gam02),'ro',real(dia1in),imag(dia1in),'r-');
    
text(-1.5,0.8,'Q0REFL7','EdgeColor','black','Fontsize',12 )
text(-1.5,0.6,['File name:'],'Fontsize',12)
text(-1.5,0.5,iname,'Fontsize',12)
text(-1.5,-0.9,'Port 0','Fontsize',12)
text(-1.5,-1.0,'UNLOADED RESONATOR','Fontsize',12)
if forced == 1
    text(-1.5,0.3,['FORCED'],'Fontsize',12)
    text(-1.5,0.2,['LOSSLESS'],'Fontsize',12)
end
hold off
axis('manual','on','xy');axis([-1.5 1.001 -1.001 1.001]);axis('equal');
axis off;

%   determine apdQ0 by scalar product
dist0=eQ0-Q0;
scprQ0=dist0'*dist0;
dQ0sp=sqrt(scprQ0/(n-3));  %10-21-2019

apdQ0=Q0std;

disp(['rs= ',num2str(rs,6),' xs= ',num2str(xs,6),...
    ' r0= ',num2str(r0,6),' g0= ',num2str(g0,6)])

%   find fL
xL=(a(2)-gam1L2)/(a(3)*gam1L2-a(1));
fL=f011*real(1+xL/2);
disp(['fL= ',num2str(fL,8),'   f0corr= ',num2str(f0corr,8),' '])

%   use analytic differentiation to find the QLstd
apgam1d=(zs-1)/(zs+1);
ftsh=ft(2:n-1);
gam1sh=gam1(2:n-1);
part1=(a(2)-a(1)/a(3))./(gam1sh-apgam1d).^2;
part2=(a(1)-a(2)*a(3))./(a(3)*ftsh+1).^2;
eQLa=-imag(part1.*part2); 
QLmean=mean(eQLa);
QLstd=std(eQLa);
QL=Q0/(1+kapa);
distL=eQLa-QL;
scprQL=distL'*distL;
dQLsp=sqrt(scprQL/(n-2));

apdQL=QLstd;

%    std value of d
ed=abs((1i*QL*ft+1).*(gam1-gam1d));
ddx=(ed-dia);
scprd=ddx'*ddx;
apdd=std(ed);
dmean=mean(ed);

%    std value of k
if forced==0
apdk=(ds/(ds-dia)^2)*apdd;
disp(['k= ',num2str(kapa,8),' dk= ',num2str(apdk)])
elseif forced==1
apdk=(2/(2-dia)^2)*apdd;
end

disp(['Q0= ',num2str(Q0),' apdQ0= ',num2str(apdQ0)])
disp(['QL= ',num2str(QL),' apdQL= ',num2str(apdQL)])


%   Check the passivity, based on measured data
g1abs=a(1)*conj(a(1));
g2abs=a(2)*conj(a(2));
g3abs=a(3)*conj(a(3));
g13abs=g1abs/g3abs;
check1=0;
if g2abs >  1e-03 && g2abs < 1.001
    check1=1;
end
check2=0;
if g13abs > 1e-03 && g13abs < 1.001
    check2=1;
end
if rs > 0
    check3=1;
else
    check3=0;
end
if QL > 0
    check4=1;
else check4=0;
end
if Q0 > 0
    check5=1;
else 
    check5=0;
end
if kapa > 0
    check6=1;
else check6=0;
end
if check1+check2+check3+check4+check5+check6==6
    disp('eq. ckt. is passive')
else
    disp('eq. ckt. is NOT  passive')
end


oname1='results.txt';


tex01=(['PRUNED from ' num2str(nraw) ' to ' num2str(npruned) ' points']);
%disp(' ')
if nraw > npruned
    disp(tex01)
end

dia1in=[gam1d gam1L2];
figure('renderer','painters','resize','off','visible','off')
SMCH(4);   %50 Ohm Smith chart
hold on;
plot(real(gam1),imag(gam1),'bo',real(gam1com),imag(gam1com),'r-',...
    real(gam1d),imag(gam1d),'ro',real(gam1c),imag(gam1c),'ro',...
    real(gam1L2),imag(gam1L2),'ro',real(dia1in),imag(dia1in),'r-');

plot(real(gam1f0),imag(gam1f0),'r.','markersize',20);

%axis('manual','on','xy');axis([-1.5 1.001 -1.001 1.001]);axis('equal');
%axis off;
text(-1.5,0.8,'Q0REFL7.m','EdgeColor','black','Fontsize',12 )
text(-1.5,0.6,['File name:'],'Fontsize',12)
text(-1.5,0.5,iname,'Fontsize',12)
if forced == 1
    text(-1.5,0.3,['FORCED'],'Fontsize',12)
    text(-1.5,0.2,['LOSSLESS'],'Fontsize',12)
end
text(-1.5,0.0,['Q_L= ',num2str(QL,5)],'Fontsize',12)
text(-1.5,-0.1,['    ',char(43),char(45),' ',num2str(QLstd,4)],'Fontsize',12)
text(-1.5,-.25,['Q_0= ',num2str(Q0,5)],'Fontsize',12)
text(-1.5,-.35,['    ',char(43),char(45),' ',num2str(Q0std,4)],'Fontsize',12)
text(-1.5,-.5,['\kappa=  ',num2str(kapa,8)],'Fontsize',12)
text(-1.5,-.6,['    ',char(43),char(45),' ',...
    num2str(apdk,4)],'Fontsize',12)
text(-1.5,-.75,['r_s= ',num2str(rs,4)],'Fontsize',12)
%text(-1.5,-.9,['f_0= ',num2str(f0corr,7),],'Fontsize',12)
text(-1.5,-.9,['f_0= ',num2str(f0corr,10)],'Fontsize',12)

if unrel==1 
    text(-1.5,-1.,'UNRELIABLE INPUT DATA!','Fontsize',12)
end
if unrel == 2
    text(-1.5,-1.,'UNRELIABLE INPUT DATA!','Fontsize',12)
end
axis('manual','on','xy');axis([-1.5 1.001 -1.001 1.001]);axis('equal');
axis off;

hold  off;
disp(' ')   %so that the next command window starts after a blank line
%*********************writing the comment****************************************
if cmnt == 'y'
comment=input('Enter a comment: ','s');
oname1='results.txt';
fod1=fopen(oname1,'a');
count1=fprintf(fod1,'\n **************************************************************************** ');
count0=fprintf(fod1,'\n Q0REFL7   %5s   yr %4i month %2i day %2i hour %2i min %2i sec %2.0f  \n',iname,clock);
count00=fprintf(fod1,'  comment: %5s \n',comment);


count11=fprintf(fod1,'         QL         dQL       kappa      dkappa          Q0        dQ0 \n');
out06=[QL QLstd kapa apdk Q0 Q0std];
count12=fprintf(fod1,'%+10.4e %+10.4e %+10.4e %+10.4e %+10.4e %+10.4e  \n',out06);

if unrel==1 && forced==0
    count01=fprintf(fod1,' UNRELIABLE DATA \n');
end
if unrel==1 && forced==1
    count01=fprintf(fod1,' UNRELIABLE DATA, FORCED LOSSLESS \n');
end
if unrel==1 && forced==2
    count01=fprintf(fod1,' UNRELIABLE DATA, SPILLOVER \n');
end     
if unrel==0 && forced==1
    count01=fprintf(fod1,'FORCED LOSSLESS \n');
end
if unrel==0 && forced==2
    count01=fprintf(fod1,'SPILLOVER \n');
end
   
if nraw > npruned
    count04=fprintf(fod1,'%s \n',tex01);
end

end














