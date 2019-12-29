function [delgam,  ndmin] = kfita( gama )
%UNTITLED Summary of this function goes here

%   estimate of the transmission line length
%   to be used by Q0REFL7.m
%   computes delgam 
%   assumption: transmission line length is not a function of frequency
%   copyrignt (c) D. Kajfez 2019

global   n  x ndeg TL  

disp('kfita has been called to check the TL length')

%   find Gd2 for theta=0
    e1=x;
    e2=ones(n,1);
    e3=-gama.*x;
    E=[e1 e2 e3];
    M=E'*E;
    h=E'*gama;
    g=M\h;
    Gd2=g(1)/g(3);    



del=zeros(ndeg+2, 1);
zs=zeros(ndeg+2,1);

for nt=1:ndeg+2     
%    thedeg(nt)=(nt-1);
    thedeg(nt)=(nt-2);
    therad=thedeg(nt)*(pi/180);
    ej=exp(+1i*(2*therad));  %moving toward load
    gamrot=gama.*ej;
%   find the new g1 to g3, fast way
    e1=x;
    e2=ones(n,1);
    e3=-gamrot.*x;
    E=[e1 e2 e3];
    M=E'*E;
    h=E'*gamrot;
    g=M\h;
    G1=(g(1)*x+g(2))./(g(3)*x+1);
    Gd=g(1)/g(3);
    
    zs(nt)=(1+Gd)/(1-Gd);   
%   from equivalent circuit    
%   more accurate rs
%    gam1d=g(1)/g(3);                        %detuned reflection coeff.
    gam1c=(conj(g(3))*g(2)-g(1))/(conj(g(3))-g(3)); %center of Q-circle
    ang1=-atan2(imag(-Gd),real(-Gd));
    ang2=atan2(imag(gam1c-Gd),real(gam1c-Gd));
    angtot=CLIP(ang1+ang2);
    cob=cos(angtot);
    ds=(1-abs(Gd)^2)/(1-abs(Gd)*cob);      
    rs = 2/ds-1;    
    xs=imag(zs(nt));
    deno=rs^2+xs^2+1+2*rs;
    gd1r=(rs^2+xs^2-1)/deno;  gd1i=2*xs/deno;
    Gd1=gd1r+1i*gd1i;
%    zs(nt)=rs+1i*xs;
    del(nt)=abs(Gd2*ej-Gd1);
end

        alpha2=atan2(imag(Gd2),real(Gd2)); 
        alphadeg2=alpha2*180/pi;
        G02=gama*exp(-1i*alpha2);   %centered on the real axis
 
%   find the first minimum
imin=zeros(ndeg,1);
im1=0;


for it=1:ndeg-1
    if del(it)>del(it+1) & del(it+2)>del(it+1)
        im1=im1+1;
        imin(im1)=it+1;
%    disp(['imin(',num2str(im1),')= ',num2str(imin(im1))])
    end
end        
        
    xsnew=imag(zs(1));
    rsnew=real(zs(1));


%disp(['size(del)= ',num2str(size(del))])
    delgam=del;
    ndmin=imin(1)-2;
%    disp(['ndmin= ',num2str(ndmin),' rsnew= ',num2str(rsnew),' xsnew= ',...
%    num2str(xsnew),' themin= ',num2str(thedeg(imin(1)))])
disp(['ndeg= ',num2str(ndeg),' theta= ',num2str(ndmin), ' degrees'])
disp('kfita finished')    

end

