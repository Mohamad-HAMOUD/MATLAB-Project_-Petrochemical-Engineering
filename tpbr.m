function [dydw,R1,R2,R3,R4]=tpbr(w,y,mass,rho0,P0,T0,FTOT0,tref,c,kc1,rhob,D)
dydw=zeros(8,1);
%dydw=[FIP]
%     [FNP] 
%     [FNA]
%     [FAR]
%     [FCR]
%     [FH2]
%     [T]
%     [P] 
R=8.314;
Tre=(1/y(7)-1/(tref+273))/R;
%Heat of reactions
DH=[-6140,55000,197560,-44440];
%Activation Energy
E=[154300,219300,154300,169800];
%Equilibrium Constants
KEQ1=kc1(1)*exp(-DH(1)*Tre);%no unit
KEQ2=kc1(2)*exp(-DH(2)*Tre);%bar 
KEQ3=kc1(3)*exp(-DH(3)*Tre);%bar^3
%Rate constants
k1=c(1)*exp(-E(1)*Tre); %Kmol/Kg(cat).h
k2=c(2)*exp(-E(2)*Tre); 
k3=c(3)*exp(-E(3)*Tre); 
k4=c(4)*exp(-E(4)*Tre);
%Partial Pressures
FTOT=y(1)+y(2)+y(3)+y(4)+y(5)+y(6);
PPIP=y(8)*y(1)/FTOT;%bar 
PPNP=y(8)*y(2)/FTOT; 
PPNA=y(8)*y(3)/FTOT; 
PPAR=y(8)*y(4)/FTOT; 
PPCR=y(8)*y(5)/FTOT; 
PPH2=y(8)*y(6)/FTOT;
%Rate laws
KADSP = 9.176;
KADSN = 19.02;
KADSA = 0;
DENOM=PPH2*((1+KADSP*((PPIP+PPNP)/PPH2)+KADSN*(PPNA)/PPH2))^2;
R1=k1*(PPNP-PPIP/KEQ1)/DENOM; %Kmol/Kg.h
R2=k2*(PPNP-PPNA*PPH2/KEQ2)/DENOM; 
R3=k3*(PPNA-PPAR*PPH2^3/KEQ3)/DENOM; 
R4=(k4*(PPIP*PPH2))/DENOM;
%mole flow
dydw(1)=R1-R4;
dydw(2)=-R1-R2;
dydw(3)=R2-R3;
dydw(4)=R3;
dydw(5)=2*R4;
dydw(6)=R2+3*R3-R4;
%mole fractions
yIP=y(1)/FTOT;
yNP=y(2)/FTOT;
yNA=y(3)/FTOT;
yAR=y(4)/FTOT;
yCR=y(5)/FTOT;
yH2=y(6)/FTOT;
%heat capacities constants
C1=[137810 81020 116600 51920 27617];
C2=[449880 345450 463810 192450 9560];
C3=[1636.9 1553.1 1672 1626.5 2466];
C4=[305300 245900 328940 116800 3700];
C5=[746.85 700.922 781.46 723.6 567];
%Gas phase heat capacity
CpIP=C1(1)+C2(1)*(C3(1)/y(7)/sinh(C3(1)/y(7)))^2+C4(1)*(C5(1)/y(7)/cosh(C5(1)/y(7)))^2;
CpNP=C1(1)+C2(1)*(C3(1)/y(7)/sinh(C3(1)/y(7)))^2+C4(1)*(C5(1)/y(7)/cosh(C5(1)/y(7)))^2;
CpNA=C1(2)+C2(2)*(C3(2)/y(7)/sinh(C3(2)/y(7)))^2+C4(2)*(C5(2)/y(7)/cosh(C5(2)/y(7)))^2;
CpAR=C1(3)+C2(3)*(C3(3)/y(7)/sinh(C3(3)/y(7)))^2+C4(3)*(C5(3)/y(7)/cosh(C5(3)/y(7)))^2;
CpCR=C1(4)+C2(4)*(C3(4)/y(7)/sinh(C3(4)/y(7)))^2+C4(4)*(C5(4)/y(7)/cosh(C5(4)/y(7)))^2;
CpH2=C1(5)+C2(5)*(C3(5)/y(7)/sinh(C3(5)/y(7)))^2+C4(5)*(C5(5)/y(7)/cosh(C5(5)/y(7)))^2;
Cp=(yIP*CpIP+yNP*CpNP+yNA*CpNA+yCR*CpCR+yH2*CpH2+yAR*CpAR)*0.001;%average heat capacity in KJ/Kmole.K
%temperature
dydw(7)=-(DH(1)*R1+DH(2)*R2+DH(3)*R3+DH(4)*R4)/(Cp*FTOT);
%catalyst properties
%u=0.00002;%Pas
Dp=0.002;%2 mm
e=0.4;%void fraction;0.62 solid fraction
% rhoc=2600;%Kg/m3
% rhob=(1-e)*rhoc;%Kg/m3
% rhob=383;
gc=1;
% D=4;%m
% H=2*D;
% rhob=(4*25000)/(pi*H*D^2);
%pressure drop
Ac=pi*(D^2)/4;%m2
beta=(mass/(rho0*gc*Dp*Ac))*((1-e)/(e^3))*(1.75*mass/Ac);%Kg/m2.s2
dydw(8)=-(beta/(Ac*rhob))*(P0/y(8))*(y(7)/T0)*(FTOT/FTOT0)*0.01;%bar
end