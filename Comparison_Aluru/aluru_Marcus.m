clear
clc
close all

%% Physical constants
kb=1.380649*10^-23;
permi= 78.38*8.854188*(10^-12);        %value at T=198.15 and P=0.1 MPa from table 4 of paper
e=1.602177*10^-19;
wallLJ=1;
ionLJ=0;

%% Process parameters

% Uncomment for 1 M
% c0 = 1*1000*6.022e23;  % #/m^3
% N=688;
% psiD=4.5; % need 0.18 C/m^2

%Uncomment for 0.25 M
c0 = 0.25*1000*6.022e23;  % #/m^3
N=688;
psiD=2.838; % need 0.12 C/m^2


electrolyte='NaCl_aluru_2_nonadjusted_LJ_PCCP_radii_Marcus';
T=298.15; % system temperature in K
rho=4/(sqrt(3)*(1.42*sqrt(3)*1e-10)^2);         %particles per unit area of the wall

%% solution parameters 
dx=0.01;                     % this dx is non-dimensional
Linf=(N-1)*dx;

%% iterations
filename=[electrolyte, '_c0_', num2str(c0/(1000*6.022e23)), '_N_', num2str(N), '_dx_', num2str(dx),'.mat'];
% initval=[electrolyte, '_c0_', num2str(c0/(1000*6.022e23)),'_intitialguess.mat'];
% if isfile(filename)
%     return;
% end
nplusf=zeros(1,N);
nminusf=zeros(1,N);
psifinal=zeros(1,N);

%% Electrolyte properties
properties=systemprops(electrolyte,wallLJ,ionLJ,T);
zplus=properties(1);
zminus=properties(2);
epsilonpp=properties(3);
epsilonmm=properties(4);
epsilonpm=properties(5);
epsilonpw=properties(6);
epsilonmw=properties(7);
sigmapp=properties(8);
sigmamm=properties(9);
sigmapm=properties(10);
sigmapw=properties(11);
sigmamw=properties(12);
aplus=properties(13);
aminus=properties(14);
aplus3=properties(15);
aminus3=properties(16);
sigmamm6=sigmamm^6;
sigmapp6=sigmapp^6;
sigmapm6=sigmapm^6;
sigmamm12=sigmamm6^2;
sigmapp12=sigmapp6^2;
sigmapm12=sigmapm6^2;
aplus3c0  =aplus3*c0;
aminus3c0 =aminus3*c0;
lambda_D = sqrt(permi*kb*T/(e^2*(zplus^2*c0*zminus+zminus^2*c0*zplus)));
a=(aplus+aminus)/2;
a3c0=a^3*c0;
four_pi_a3c0 = 4*pi*a3c0;
factorpsi=(e^2*a^2*c0/(kb*T*permi));
if(aplus>aminus)
    aratio=(aplus3c0/aminus3c0)^(1/3);
else
    aratio=(aminus3c0/aplus3c0)^(1/3);
end

%% Solve
xmesh=linspace(0,Linf,N);
guess = @(x) [psiD*exp(-x/lambda_D), -(psiD/lambda_D)*exp(-x/lambda_D)];
solinit=bvpinit(xmesh,guess);
options= bvpset('stats','off','RelTol',0.001,'Stats','on');
sol=bvp5c(@(x,y)odefcn(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,rho,Linf,a),@(psia,psib) aluru_bcfcn(psia,psib,psiD),solinit,options);
y=deval(sol,xmesh);

sol.y=y(1,:);
sol.x=xmesh;

psi=sol.y(1,:);

awall=1.42e-10;
ap=((aplus3c0/c0)^(1/3));
am=((aminus3c0/c0)^(1/3));
apavg=(ap+awall)/2;
amavg=(am+awall)/2;
apm=(ap+am)/2;

ULJplus=-4*pi*rho*epsilonpw.*(apm^2)*((((sigmapw.^6)./(2.*(xmesh).^4))-((sigmapw.^12)./(5.*(xmesh).^10)))...
        +((sigmapw.^6./(2.*((Linf)-(xmesh)).^4))-((sigmapw.^12./(5.*((Linf)-(xmesh)).^10)))))...
        +4*pi*rho*epsilonpw.*(apm^2)*((((sigmapw.^6)./(2.*(Linf/2).^4))-((sigmapw.^12)./(5.*(Linf/2).^10)))...
        + ((sigmapw.^6./(2.*((Linf/2)).^4))-((sigmapw.^12./(5.*((Linf/2)).^10)))));

ULJminus=-4*pi*rho*epsilonmw.*(apm^2)*((((sigmamw.^6)./(2.*(xmesh).^4))-((sigmamw.^12)./(5.*(xmesh).^10)))...
         + ((sigmamw.^6./(2.*((Linf)-(xmesh)).^4))-((sigmamw.^12./(5.*((Linf)-(xmesh)).^10)))))...
         -4*pi*rho*epsilonmw.*(apm^2)*((((sigmamw.^6)./(2.*(Linf/2).^4))-((sigmamw.^12)./(5.*(Linf/2).^10)))...
         + ((sigmamw.^6./(2.*((Linf/2)).^4))-((sigmamw.^12./(5.*((Linf/2)).^10)))));

fofpsi= (1+(zminus*aplus3c0*(exp(-zplus*psi - ULJplus)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
gofpsi= fofpsi + (zplus*aminus3c0*(exp(zminus*psi - ULJminus)-fofpsi)) + (zminus*aplus3c0*fofpsi.*(exp(-zplus*psi - ULJplus)-1));
nplus=fofpsi.*exp(-zplus*psi - ULJplus)./gofpsi;
nminus=exp(zminus*psi - ULJminus)./gofpsi;

xmesh=linspace(0,Linf,N);
guess = @(x) [psiD*exp(-x/lambda_D), -(psiD/lambda_D)*exp(-x/lambda_D)];
solinit=bvpinit(xmesh,guess);
options= bvpset('stats','off','RelTol',0.001,'Stats','on');
sol=bvp5c(@(x,y)odefcn_0(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,rho,Linf,a),@(psia,psib) aluru_bcfcn(psia,psib,psiD),solinit,options);
y=deval(sol,xmesh);
 
sol.y=y(1,:);
sol.x=xmesh;
 
psi_zero=sol.y(1,:);
 
 
ULJminus_zero=zeros(size(ULJminus));
ULJplus_zero=zeros(size(ULJplus));
fofpsi_zero= (1+(zminus*aplus3c0*(exp(-zplus*psi_zero - ULJplus_zero)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
gofpsi_zero= fofpsi_zero + (zplus*aminus3c0*(exp(zminus*psi_zero - ULJminus_zero)-fofpsi_zero)) + (zminus*aplus3c0*fofpsi_zero.*(exp(-zplus*psi_zero - ULJplus_zero)-1));
nplus_zero=fofpsi_zero.*exp(-zplus*psi_zero - ULJplus_zero)./gofpsi_zero;
nminus_zero=exp(zminus*psi_zero - ULJminus_zero)./gofpsi_zero;


set(0,'defaulttextinterpreter','latex');
figure;
plot(sol.x, psi,'DisplayName','with LJ','Color','#7E2F8E');
hold on
plot(sol.x, psi_zero,'DisplayName','no LJ','Color','#7E2F8E','Linestyle','--');
legend boxoff
xlabel('X');
ylabel({'$\psi$'});
f = gcf;
%exportgraphics(f,'Na2SO4_1.png','Resolution',600)

figure;
plot(sol.x, nplus*c0*zminus,'DisplayName','with LJ','Color','#7E2F8E');
hold on
plot(sol.x, nplus_zero*c0*zminus,'DisplayName','no LJ','Color','#7E2F8E','Linestyle','--');
legend boxoff
xlabel('X');
ylabel({'$\c+$'});
f = gcf;
%exportgraphics(f,'Na2SO4_2.png','Resolution',600)

figure;
plot(sol.x, nminus*c0*zplus,'DisplayName','with LJ','Color','#7E2F8E');
hold on
plot(sol.x, nminus_zero*c0*zplus,'DisplayName','no LJ','Color','#7E2F8E','Linestyle','--');
legend boxoff
xlabel('X');
ylabel({'$\c-$'});
f = gcf;
%exportgraphics(f,'Na2SO4_3.png','Resolution',600)
separation=Linf*a
surf_charge_density=permi*y(2,1)*(kb*T/e)/a
final=[sol.x', nplus', nplus_zero', nminus', nminus_zero', psi', psi_zero'];
writematrix(final,"Revised_NaCl_aluru_2_c0_0.25.csv");