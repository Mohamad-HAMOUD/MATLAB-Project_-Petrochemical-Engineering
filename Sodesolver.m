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
mfeed =150;%ton/h
HHR = 3.5;
temp0 = 490;
T0=temp0+273;
R=8.314;
WHSV=mfeed/100;
%initial molar flow rates
fhc=0;
for k=1:4
   fj0(k) =  ((wt0(k)*mfeed)/mw(k))*1000;
   fhc = fhc+fj0(k);
end
FIP0 = fj0(1);%Kmol/h
FNP0 = fj0(2);
FNA0 = fj0(3);
FAR0 = fj0(4);
FCR0 = 0;
FH20 = HHR*fhc;
FTOT0=FIP0+FNP0+FNA0+FAR0+FCR0+FH20;

P0=8;%bar
%Molar Mass
mwhc = 1/(((wt0(1)+wt0(2))/mw(1))+(wt0(3)/mw(3))+(wt0(4)/mw(4)));
mass=(mfeed/3.6)+ (FH20*mwh2/(3600));%Kg/s
avmw=(FH20/FTOT0)*mwh2 + (fhc/FTOT0)*mwhc;%Kg/Kmol
rho0=(P0*100*avmw)/(R*T0);%density of entering feed Kg/dm3
%Reactor section
%r=L/R;
r=0.7;%chosen value 
rhob=383;
W=25000;%kg
R3=(3*W)/(rhob*pi*(6*r-2*r^3));
D=2*nthroot(R3,3);
L=r*(D/2);

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
% solving ODE
%Reactor 1
[z1,y1]=ode15s(@(z,y) spbr(z,y,mass,rho0,P0,T0,FTOT0,tref,c,kc1,L,D,rhob),[0 2*L],[FIP0,FNP0,FNA0,FAR0,FCR0,FH20,T0,P0]);
w1=rhob*pi*((D^2/4).*z1-(z1-L).^3/3-L^3/3);
nrow = size(z1,1);
r11=zeros(nrow,1);
r21=zeros(nrow,1);
r31=zeros(nrow,1);
r41=zeros(nrow,1);
for row = 1 : nrow
    [~, r11(row),r21(row),r31(row),r41(row)] = spbr(z1(row), y1(row,:),mass,rho0,P0,T0,FTOT0,tref,c,kc1,L,D,rhob);
end
Tav(1)=(T0+y1(end,7))/2;
Pav(1)=(P0+y1(end,8))/2;
%Reactor 2
[z2,y2]=ode15s(@(z,y) spbr(z,y,mass,rho0,P0,T0,FTOT0,tref,c,kc1,L,D,rhob),[0 2*L],[y1(end,1),y1(end,2),y1(end,3),y1(end,4),y1(end,5),y1(end,6),T0,y1(end,8)]);
w2=w1(end) + rhob*pi*((D^2/4).*z2-(z2-L).^3/3-L^3/3);
nrow = size(w2,1);
r12=zeros(nrow,1);
r22=zeros(nrow,1);
r32=zeros(nrow,1);
r42=zeros(nrow,1);
for row = 1 : nrow
    [~, r12(row),r22(row),r32(row),r42(row)] = spbr(w2(row), y2(row,:),mass,rho0,P0,T0,FTOT0,tref,c,kc1,L,D,rhob);
end
Tav(2)=(y2(1,7)+y2(end,7))/2;
Pav(2)=(y1(end,8)+y2(end,8))/2;
%Reactor 3
[z3,y3]=ode15s(@(z,y) spbr(z,y,mass,rho0,P0,T0,FTOT0,tref,c,kc1,L,D,rhob),[0 2*L],[y2(end,1),y2(end,2),y2(end,3),y2(end,4),y2(end,5),y2(end,6),T0,y2(end,8)]);
w3=w2(end) + rhob*pi*((D^2/4).*z3-(z3-L).^3/3-L^3/3);
nrow = size(w3,1);
r13=zeros(nrow,1);
r23=zeros(nrow,1);
r33=zeros(nrow,1);
r43=zeros(nrow,1);
for row = 1 : nrow
    [~, r13(row),r23(row),r33(row),r43(row)] = spbr(w3(row), y3(row,:),mass,rho0,P0,T0,FTOT0,tref,c,kc1,L,D,rhob);
end
Tav(3)=(y3(1,7)+y3(end,7))/2;
Pav(3)=(y2(end,8)+y3(end,8))/2;
%Reactor 4
[z4,y4]=ode15s(@(z,y) spbr(z,y,mass,rho0,P0,T0,FTOT0,tref,c,kc1,L,D,rhob),[0 2*L],[y3(end,1),y3(end,2),y3(end,3),y3(end,4),y3(end,5),y3(end,6),T0,y3(end,8)]);
w4=w3(end) + rhob*pi*((D^2/4).*z4-(z4-L).^3/3-L^3/3);
nrow = size(w4,1);
r14=zeros(nrow,1);
r24=zeros(nrow,1);
r34=zeros(nrow,1);
r44=zeros(nrow,1);
for row = 1 : nrow
    [~, r14(row),r24(row),r34(row),r44(row)] = spbr(w4(row), y4(row,:),mass,rho0,P0,T0,FTOT0,tref,c,kc1,L,D,rhob);
end
Tav(4)=(y4(1,7)+y4(end,7))/2;
Pav(4)=(y3(end,8)+y4(end,8))/2;
%%%%
w=[w1;w2;w3;w4];
y=[y1;y2;y3;y4];
r1=[r11;r12;r13;r14];
r2=[r21;r22;r23;r24];
r3=[r31;r32;r33;r34];
r4=[r41;r42;r43;r44];
%yield in weight percent
RNDIP=100*y(:,1)*mw(1)/(mfeed*1000); 
RNDNP=100*y(:,2)*mw(2)/(mfeed*1000); 
RNDNA=100*y(:,3)*mw(3)/(mfeed*1000); 
RNDAR=100*y(:,4)*mw(4)/(mfeed*1000); 
RNDCR=100*y(:,5)*mw(5)/(mfeed*1000); 
RNDH2=100*mwh2*(y(:,6)-FH20)/(mfeed*1000); 
RNDC5P=100-RNDCR-RNDH2;
%RON
RONNP=0;RONIP = 79;RONNA = 68;RONAR = 125;
RON=(RONIP.*RNDIP+RONNP.*RNDNP+RONNA.*RNDNA+RONAR.*RNDAR)./RNDC5P;
for n=1:4
coke(n)=4.99*10^6*exp(-8955/Tav(n))*(Pav(n)*0.1)^-0.94*(WHSV^-1.28)*(HHR^-1.33);
end
deltaP=P0-y(end,8);

 
%Plotting
% figure(1)
% % % %subplot(2,2,1) 
% plot(w/1000,y(:,7)-273)
% % plot(w/1000,y(:,7)-273,':')%dotted for comparision
% grid on
% xlabel('weight of catalyst(ton)')
% ylabel('Temperature(°C)')
% legend('TEMP')
% axis([0 100 400 500])
% % %subplot(2,2,2)
% figure(2)
% plot(w/1000,RNDIP,'b:',w/1000,RNDNP,'m:',w/1000,RNDNA,'r:',w/1000,RNDAR,'c:')
% grid on
% xlabel('weight of catalyst(t)')
% ylabel('Yield(wt%)')
% legend('RNDIP','RNDNP','RNDNA','RNDAR')
% axis([0 100 0 50])
% % %subplot(2,2,3)
% figure(3)
% plot(w/1000,RNDH2,'b:',w/1000,RNDCR,'r:')
% grid on
% xlabel('Weight of catalyst (t)')
% ylabel('Yield (wt%)')
% legend('RNDH2','RNDCR')
% figure(4)
% plot(w/1000,RNDC5P,'b:')
% grid on
% xlabel('Weight of catalyst (t)')
% ylabel('Reformate yield (wt%)'
% figure(5)
% plot(w/1000,RON,'b:')
% xlabel('weight of catalyst(t)')
% grid on 
% ylabel('Research Octane Number')
% legend('RON')
% axis([0 100 50 100])
% figure(6)
% plot(w/1000,r1,'b',w/1000,r2,'m',w/1000,r3,'r',w/1000,r4,'c')
% grid on
% xlabel('Weight of catalyst (t)')
% ylabel('Reaction rate (mole/gcat/hr)')
% legend('R1','R2','R3','R4')
% axis([0 100 -0.005 0.03])
% title('Spherical model')
% figure(7)
% plot(w/1000,y(:,8),':')
% grid on
% xlabel('weight of catalyst(t)')
% ylabel('Total Pressure(bar)')
% title('Pressure drop')
% axis([0 100 4 14])
% figure(8)
% plot(w/1000,y(:,1),'b:',w/1000,y(:,2),'m:',w/1000,y(:,3),'r:',w/1000,y(:,4),'c:')
% grid on
% xlabel('weight of catalyst(t)')
% ylabel('Molar Flow(Kmol/h)')
%legend('iso-paraffin','n-paraffin','naphthene','aromatics')
% figure(9)
% plot(w/1000,y(:,5),w/1000,y(:,6))
% %plot(w/1000,y(:,6))
% grid on
% xlabel('weight of catalyst(t)')
% ylabel('Molar Flow(Kmol/h)')
% %ylabel(' Hydrogen Molar Flow(Kmol/h)')
% legend('cracked products','hydrogen')

%for comparison only
%YIELDS IN 4TH REACTOR
% %aromatic yield
% RNDAR4=100*y4(:,4)*mw(4)/(mfeed*1000);
% %hydrogen yield
% RNDH24=100*mwh2*(y4(:,6)-FH20)/(mfeed*1000); 
% %Cracked products
% RNDCR4=100*y4(:,5)*mw(5)/(mfeed*1000);  
% %Reformate yield
% RNDC5P4=100-RNDCR4-RNDH24;
% %RON 
% s=(size(w,1)-size(w4,1))+1;
% RON4=RON(s:end);
% %Coke
% n=[1 2 3 4];
% plot(n,coke,':*')

%Geometry Drawing
% R=D/2;%radius
% [x,y,z]=sphere;
% axis equal
% X2=R*x;
% Y2=R*y;
% Z2=R*z;
% mesh(X2,Y2,Z2,'EdgeColor','k')
% axis equal
% hold on
% radius=sqrt(R^2-L^2);
% th=0:pi/100:2*pi;
% a=radius*cos(th);
% b=radius*sin(th);
% H1=0.7*R;
% H2=-0.7*R;
% plot3(a, b, ones(1,size(th,2))*H1,'-r')
% hold on
% plot3(a, b, ones(1,size(th,2))*H2, '-r')
% xlabel('(m)')
% ylabel('(m)')
% zlabel('(m)')
% title('Spherical reactor shape')
% 
