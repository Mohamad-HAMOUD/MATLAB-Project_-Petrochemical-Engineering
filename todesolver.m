clear;clc;
wt0(1) = 0.315;%ip0
wt0(2) = 0.303;%np0
wt0(3) = 0.266;
wt0(4) = 0.116;
mwh2=2.016;
mw(1) = 100.206;
mw(2) = 100.206;
mw(3) = 98.1904;
mw(4) = 92.1424;
mw(5) = 51.1112;
mfeed =150;
HHR = 3.5;
temp0 = 490;
T0=temp0+273;
R=8.314;
WHSV=mfeed/100;
%initial molar flow rates
fhc=0;
for k=1:4
   fj0(k) =  (wt0(k)*mfeed)/mw(k);
   fhc = fhc+fj0(k);
end
FIP0 = fj0(1);%10^6 mol/h
FNP0 = fj0(2);
FNA0 = fj0(3);
FAR0 = fj0(4);
FCR0 = 0;
FH20 = HHR*fhc;
FTOT0=FIP0+FNP0+FNA0+FAR0+FCR0+FH20;

P0=18;
%Molar Mass
mwhc = 1/(((wt0(1)+wt0(2))/mw(1))+(wt0(3)/mw(3))+(wt0(4)/mw(4)));
mass=(mfeed/3.6)+ ((FH20*mwh2*1000)/3600);%Kg/s
avmw=(FH20/FTOT0)*mwh2 + (fhc/FTOT0)*mwhc;
rho0=(P0*100*avmw)/(R*T0);%density of entering feed
%equilibrium and rate constants depending on Tref
Tref=input('Please enter a letter of your choice of Tref in degree celcius:\nA: 480\nB: 500\nC: 520\n','s');
switch Tref
    case 'A'
        tref=480;
        c(1) = 4.562;
        c(2) = 0.263;
        c(3) = 29.539;
        c(4) = 0.008;
        kc1(1) = 4.479;
        kc1(2) = 5.217;
        kc1(3) = 150043;
        
    case 'B'
        tref=500;
        c(1) = 8.63;
        c(2) = 0.876;
        c(3) = 55.875;
        c(4) = 0.016;
        kc1(1) = 4.367;
        kc1(2) = 6.548;
        kc1(3) = 339356;
        
    case 'C'
        tref=520;
        c(1) = 15.808;
        c(2) = 2.748;
        c(3) = 102.349;
        c(4) = 0.032;
        kc1(1) = 4.263;
        kc1(2) = 8.124;
        kc1(3) = 736580;
end

rhob=383;
D=4;
W=25000;
rhobmax=(W*4)/(pi*D^3);
%solving ODE
%Reactor 1
[w1,y1]=ode45(@(w,y) tpbr(w,y,mass,rho0,P0,T0,FTOT0,tref,c,kc1,rhob,D),[0 25],[FIP0,FNP0,FNA0,FAR0,FCR0,FH20,T0,P0]);
nrow = size(w1,1);
r11=zeros(nrow,1);
r21=zeros(nrow,1);
r31=zeros(nrow,1);
r41=zeros(nrow,1);
for row = 1 : nrow
    [~, r11(row),r21(row),r31(row),r41(row)] = tpbr(w1(row), y1(row,:),mass,rho0,P0,T0,FTOT0,tref,c,kc1,rhob,D);
end
Tav(1)=(T0+y1(end,7))/2;
Pav(1)=(P0+y1(end,8))/2;
%[w2,y2]=ode45('tpbr',[25 50],[y1(end,1),y1(end,2),y1(end,3),y1(end,4),y1(end,5),y1(end,6),T0,y1(end,8)]);
[w2,y2]=ode45(@(w,y) tpbr(w,y,mass,rho0,P0,T0,FTOT0,tref,c,kc1,rhob,D),[25 50],[y1(end,1),y1(end,2),y1(end,3),y1(end,4),y1(end,5),y1(end,6),T0,y1(end,8)]);
nrow = size(w2,1);
r12=zeros(nrow,1);
r22=zeros(nrow,1);
r32=zeros(nrow,1);
r42=zeros(nrow,1);
for row = 1 : nrow
    [~, r12(row),r22(row),r32(row),r42(row)] = tpbr(w2(row), y2(row,:),mass,rho0,P0,T0,FTOT0,tref,c,kc1,rhob,D);
end
Tav(2)=(y2(1,7)+y2(end,7))/2;
Pav(2)=(y1(end,8)+y2(end,8))/2;
%[w3,y3]=ode45('tpbr',[50 75],[y2(end,1),y2(end,2),y2(end,3),y2(end,4),y2(end,5),y2(end,6),T0,y2(end,8)]);
[w3,y3]=ode45(@(w,y) tpbr(w,y,mass,rho0,P0,T0,FTOT0,tref,c,kc1,rhob,D),[50 75],[y2(end,1),y2(end,2),y2(end,3),y2(end,4),y2(end,5),y2(end,6),T0,y2(end,8)]);
nrow = size(w3,1);
r13=zeros(nrow,1);
r23=zeros(nrow,1);
r33=zeros(nrow,1);
r43=zeros(nrow,1);
for row = 1 : nrow
    [~, r13(row),r23(row),r33(row),r43(row)] = tpbr(w3(row), y3(row,:),mass,rho0,P0,T0,FTOT0,tref,c,kc1,rhob,D);
end
Tav(3)=(y3(1,7)+y3(end,7))/2;
Pav(3)=(y2(end,8)+y3(end,8))/2;
%[w4,y4]=ode45('tpbr',[75 100],[y3(end,1),y3(end,2),y3(end,3),y3(end,4),y3(end,5),y3(end,6),T0,y3(end,8)]);
[w4,y4]=ode45(@(w,y) tpbr(w,y,mass,rho0,P0,T0,FTOT0,tref,c,kc1,rhob,D),[75 100],[y3(end,1),y3(end,2),y3(end,3),y3(end,4),y3(end,5),y3(end,6),T0,y3(end,8)]);
nrow = size(w4,1);
r14=zeros(nrow,1);
r24=zeros(nrow,1);
r34=zeros(nrow,1);
r44=zeros(nrow,1);
for row = 1 : nrow
    [~, r14(row),r24(row),r34(row),r44(row)] = tpbr(w4(row), y4(row,:),mass,rho0,P0,T0,FTOT0,tref,c,kc1,rhob,D);
end
Tav(4)=(y4(1,7)+y4(end,7))/2;
Pav(4)=(y3(end,8)+y4(end,8))/2;
w=[w1;w2;w3;w4];
y=[y1;y2;y3;y4];
r1=[r11;r12;r13;r14];
r2=[r21;r22;r23;r24];
r3=[r31;r32;r33;r34];
r4=[r41;r42;r43;r44];
%yield in weight percent
RNDIP=100*y(:,1)*mw(1)/(mfeed); 
RNDNP=100*y(:,2)*mw(2)/(mfeed); 
RNDNA=100*y(:,3)*mw(3)/(mfeed); 
RNDAR=100*y(:,4)*mw(4)/(mfeed); 
RNDCR=100*y(:,5)*mw(5)/(mfeed); 
RNDH2=100*mwh2*(y(:,6)-FH20)/(mfeed); 
RNDC5P=100-RNDCR-RNDH2;
%RON
RONNP=0;RONIP = 79;RONNA = 68;RONAR = 125;
RON=(RONIP*RNDIP+RONNP*RNDNP+RONNA*RNDNA+RONAR*RNDAR)./RNDC5P;
for n=1:4
coke(n)=4.99*10^6*exp(-8955/Tav(n))*(Pav(n)*0.1)^-0.94*WHSV^-1.28*HHR^-1.33;
end
deltaP=P0-y(end,8);

% %Plotting
% figure(1)
% %subplot(2,2,1) 
% plot(w,y(:,7)-273)
% grid on
% xlabel('weight of catalyst(t)')
% ylabel('Temperature(°C)')
% legend('TEMP')
% axis([0 100 400 500])
% %subplot(2,2,2)
% figure(2)
% plot(w,RNDIP,'b',w,RNDNP,'m',w,RNDNA,'r',w,RNDAR,'c')
% grid on
% xlabel('Weight of catalyst (t)')
% ylabel('Yield (wt%)')
% legend('RNDIP','RNDNP','RNDNA','RNDAR')
% %subplot(2,2,3
% figure(3)
% plot(w,RNDH2,'b',w,RNDCR,'r')
% grid on
% xlabel('Weight of catalyst (t)')
% ylabel('Yield (wt%)')
% legend('RNDH2','RNDCR')
% figure(4)
% plot(w,RNDC5P)
% grid on
% xlabel('Weight of catalyst (t)')
% ylabel('Reformate yield (wt%)')
% figure(5)
% plot(w,RON)
% xlabel('Weight of catalyst(t)')
% grid on
% ylabel('Research Octane Number')
% legend('RON')
% %subplot(2,2,4)
% figure(6)
% plot(w,r1,'b',w,r2,'m',w,r3,'r',w,r4,'c')
% grid on
% xlabel('Weight of catalyst (t)')
% ylabel('Reaction rate (mole/gcat/hr)')
% legend('R1','R2','R3','R4')
% axis([0 100 -0.005 0.03])
% title('Tubular model')
% figure(7)
% plot(w,y(:,8))
% grid on
% xlabel('weight of catalyst(t)')
% ylabel('Total Pressure(bar)')
% title('Pressure drop')
% legend('TPBR')


%for comparison only
%YIELDS IN 4TH REACTOR
% %aromatic yield
% RNDAR4=100*y4(:,4)*mw(4)/(mfeed);
% %plot(w4,RNDAR4)
% %hydrogen yield
% RNDH24=100*mwh2*(y4(:,6)-FH20)/(mfeed); 
% %Cracked products
% RNDCR4=100*y4(:,5)*mw(5)/(mfeed);  
% %Reformate yield
% RNDC5P4=100-RNDCR4-RNDH24;
% %RON 
% s=(size(w,1)-size(w4,1))+1;
% RON4=RON(s:end);
% %coke
% %Coke
% n=[1 2 3 4];
% plot(n,coke,':*')

%Geometry drawing
% R=D/2;
% H=W/(rhob*pi*R^2);
% %[x,y,z]=cylinder(r);
% %z=z*H;
% %surf(x,y,z)
% x=0;
% y=0;
% th=0:pi/100:2*pi;
% a=R*cos(th);
% b=R*sin(th);
% figure
% surf([a; a]+x, [b; b]+y, [ones(1,size(th,2)); zeros(1,size(th,2))]*H, 'FaceColor','b', 'EdgeColor','none')
% % hold on
% % fill3(a+x, b+y, ones(1,size(th,2))*H,'-b')
% % fill3(a+x, b+y, zeros(1,size(th,2)), '-b')
% xlabel('(m)')
% ylabel('(m)')
% zlabel('(m)')
% title('Tubular reactor shape')
