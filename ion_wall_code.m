clear
clc
close all

%% Physical constants
kb=1.380649*10^-23;					% boltzmann constant
permi= 78.38*8.854188*(10^-12);			% value at T=298.15 and P=0.1 MPa
e=1.602177*10^-19;					% electronic charge
wallLJ=1;
ionLJ=0;

%% Process parameters
c0 = 5*1000*6.022e23;					% bulk electrolyte concentration in #/m^3
psiD=5;						% wall potential
electrolyte='NaCl';
T=298.15;						% system temperature in K
rho=4/(sqrt(3)*(1.42*sqrt(3)*1e-10)^2);		% particles per unit area of the wall

%% solution parameters 
N=501;							% number of grid points
dx=0.1;						% this dx is non-dimensional
Linf=10;						% Distance between the plates

nplusf=zeros(1,N);					% initialization of solution array
nminusf=zeros(1,N);
psifinal=zeros(1,N);

%% Electrolyte properties
properties=systemprops(electrolyte,wallLJ,ionLJ,T);
zplus=properties(1);
zminus=properties(2);
epsilonpp=properties(3);                        % LJ parameter epsilon for cation-cation interactions
epsilonmm=properties(4);                        % LJ parameter epsilon for anion-anion interactions
epsilonpm=properties(5);                        % LJ parameter epsilon for cation-anion interactions
epsilonpw=properties(6);                        % LJ parameter epsilon for cation-wall interactions
epsilonmw=properties(7);                        % LJ parameter epsilon for anion-wall interactions
sigmapp=properties(8);                          % LJ parameter sigma for cation-cation interactions
sigmamm=properties(9);                          % LJ parameter sigma for anion-anion interactions
sigmapm=properties(10);                         % LJ parameter sigma for cation-anion interactions
sigmapw=properties(11);                         % LJ parameter sigma for cation-wall interactions
sigmamw=properties(12);                         % LJ parameter sigma for anion-wall interactions
aplus=properties(13);                           % cation diameter
aminus=properties(14);                          % anion diameter
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
lambda_D = sqrt(permi*kb*T/(e^2*(zplus^2*c0*zminus+zminus^2*c0*zplus)));		% definition of debye length
a=(aplus+aminus)/2;
a3c0=a^3*c0;
four_pi_a3c0 = 4*pi*a3c0;
factorpsi=(e^2*a^2*c0/(kb*T*permi));
aratio=(aminus3/aplus3)^(1/3);

%% Solve
xmesh=linspace(0,Linf,N);								% definition of the 1-D geometry
guess = @(x) [psiD*exp(-x/lambda_D), -(psiD/lambda_D)*exp(-x/lambda_D)];		% intial guess for the solver
solinit=bvpinit(xmesh,guess);								% guess for bvp solver
options= bvpset('stats','off','RelTol',0.001,'Stats','on');
sol=bvp5c(@(x,y)odefcn(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,rho,Linf,a),@(psia,psib) bcfcn(psia,psib,psiD),solinit,options);
y=deval(sol,xmesh);									% obtaining solution in arrays

sol.y=y(1,:);
sol.x=xmesh;

psi=sol.y(1,:);

awall=1.42e-10;
ap=((aplus3c0/c0)^(1/3));
am=((aminus3c0/c0)^(1/3));
apavg=(ap+awall)/2;
amavg=(am+awall)/2;
apm=(ap+am)/2;

% equations obtained from chemical potential expressions expression

ULJplus=-4*pi*rho*epsilonpw.*(apm^2)*((((sigmapw.^6)./(2.*(xmesh).^4))-((sigmapw.^12)./(5.*(xmesh).^10)))...
    + ((sigmapw.^6./(2.*((Linf)-(xmesh)).^4))-((sigmapw.^12./(5.*((Linf)-(xmesh)).^10)))));
ULJminus=-4*pi*rho*epsilonmw.*(apm^2)*((((sigmamw.^6)./(2.*(xmesh).^4))-((sigmamw.^12)./(5.*(xmesh).^10)))...
    + ((sigmamw.^6./(2.*((Linf)-(xmesh)).^4))-((sigmamw.^12./(5.*((Linf)-(xmesh)).^10)))));

if (am > ap) 		% Chemical Potential Expressions with anion diameter greater than cation diameter
    fofpsi= (1+(zminus*aplus3c0*(exp(-zplus*psi - ULJplus)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zplus*aminus3c0*(exp(zminus*psi - ULJminus)-fofpsi)) + (zminus*aplus3c0*fofpsi.*(exp(-zplus*psi - ULJplus)-1));
    nplus=fofpsi.*exp(-zplus*psi - ULJplus)./gofpsi;
    nminus=exp(zminus*psi - ULJminus)./gofpsi;
else			% Chemical Potential Expressions with cation diameter greater than anion diameter
    fofpsi= (1+(zplus*aminus3c0*(exp(zminus*psi - ULJminus)-1)./(1-zminus*aplus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zminus*aplus3c0*(exp(-zplus*psi - ULJplus)-fofpsi)) + (zplus*aminus3c0*fofpsi.*(exp(zminus*psi - ULJminus)-1));
    nplus=exp(-zplus*psi - ULJplus)./gofpsi;
    nminus=fofpsi.*exp(zminus*psi - ULJminus)./gofpsi;
end

%%%% Calculations without ion-wall LJ interactions (U_LJ set to zero)

xmesh=linspace(0,Linf,N);
guess = @(x) [psiD*exp(-x/lambda_D), -(psiD/lambda_D)*exp(-x/lambda_D)];
solinit=bvpinit(xmesh,guess);
options= bvpset('stats','off','RelTol',0.001);
sol=bvp5c(@(x,y)odefcn_0(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,rho,Linf,a),@(psia,psib) bcfcn(psia,psib,psiD),solinit,options);
y=deval(sol,xmesh);

sol.y=y(1,:);
sol.x=xmesh;

psi_zero=sol.y(1,:);


ULJminus_zero=zeros(size(ULJminus));
ULJplus_zero=zeros(size(ULJplus));

if (am > ap) 		% Chemical Potential Expressions with anion diameter greater than cation diameter
    fofpsi_zero= (1+(zminus*aplus3c0*(exp(-zplus*psi_zero - ULJplus_zero)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
    gofpsi_zero= fofpsi_zero + (zplus*aminus3c0*(exp(zminus*psi_zero - ULJminus_zero)-fofpsi_zero)) + (zminus*aplus3c0*fofpsi_zero.*(exp(-zplus*psi_zero - ULJplus_zero)-1));
    nplus_zero=fofpsi_zero.*exp(-zplus*psi_zero - ULJplus_zero)./gofpsi_zero;
    nminus_zero=exp(zminus*psi_zero - ULJminus_zero)./gofpsi_zero;
else			% Chemical Potential Expressions with cation diameter greater than anion diameter
    fofpsi_zero= (1+(zplus*aminus3c0*(exp(zminus*psi_zero - ULJminus_zero)-1)./(1-zminus*aplus3c0))).^(aratio^3-1);
    gofpsi_zero= fofpsi_zero + (zminus*aplus3c0*(exp(-zplus*psi - ULJplus)-fofpsi)) + (zplus*aminus3c0*fofpsi.*(exp(zminus*psi - ULJminus)-1));
    nplus_zero=exp(-zplus*psi_zero - ULJplus_zero)./gofpsi_zero;
    nminus_zero=fofpsi_zero.*exp(zminus*psi_zero - ULJminus_zero)./gofpsi_zero;
end

fofpsi_zero= (1+(zminus*aplus3c0*(exp(-zplus*psi_zero - ULJplus_zero)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
gofpsi_zero= fofpsi_zero + (zplus*aminus3c0*(exp(zminus*psi_zero - ULJminus_zero)-fofpsi_zero)) + (zminus*aplus3c0*fofpsi_zero.*(exp(-zplus*psi_zero - ULJplus_zero)-1));
nplus_zero=fofpsi_zero.*exp(-zplus*psi_zero - ULJplus_zero)./gofpsi_zero;
nminus_zero=exp(zminus*psi_zero - ULJminus_zero)./gofpsi_zero;

% plotting electric potential profile with and without LJ interaction
figure;
set(0,'defaulttextinterpreter','latex');
plot(sol.x, psi,'b','DisplayName','Wall LJ');
hold on
plot(sol.x, psi_zero,'g','DisplayName','No LJ');
legend boxoff
xlabel('X');
ylabel({'$\psi$'});
box on

% plotting cation concentration profile with and without LJ interaction 
figure;
set(0,'defaulttextinterpreter','latex');
plot(sol.x, nplus,'b','DisplayName','Wall LJ');
hold on
plot(sol.x, nplus_zero,'g','DisplayName','No LJ');
legend boxoff
% xlim([0 10])
xlabel('X');
ylabel({'n$_{+}$'});
box on

 % plotting anion concentration profile with and without LJ interaction
figure;
set(0,'defaulttextinterpreter','latex');
plot(sol.x, nminus,'b','DisplayName','Wall LJ');
hold on
plot(sol.x, nminus_zero,'g','DisplayName','No LJ');
legend boxoff
% xlim([0 10])
xlabel('X');
ylabel({'n$_{-}$'});

box on
%grid on
%set(gca,'GridLineStyle','--')
%ha.FontName = Arial
%ha.FontSize = 11
%ylim([0,10])
f=gcf;
%size = [400 400];
%res = 650;
%set(f,'paperunits','inches','paperposition',[0 0 size/res]);
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'bar','-dpdf','-r0','-bestfit')
exportgraphics(f,'bar.png','Resolution',650,'ContentType','vector')



%%
function properties=systemprops(electrolyte,wallLJ,ionLJ,T)

factoreps=1000/(8.314*T);

if (strcmpi(electrolyte,'MgCl2'))
    % MgCl2 (2:1)
    zplus=2;
    zminus=1;
    
    aplus=0.13e-9;
    aplus3  = aplus^3;
    aminus=0.362e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=3.651900*factoreps;       % value is 1.47398
    epsilonmm=0.076923*factoreps;       % value is 0.03104
    epsilonpm=3.000*factoreps;
      
    sigmapp=(0.116290e-9)/a;         % value is 0.3743
    sigmamm=(0.469906e-9)/a;         % value is 1.5123
    sigmapm=0.300e-9/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        % epsilon of carbon is 0.2330488 kJ/mol
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      % sigma of carbon is 3.400 angstroms
    sigmamw=(sigmamm+sigmaww)/2;

elseif (strcmpi(electrolyte,'NaCl'))
    
    % NaCl (1:1)
    zplus=1;
    zminus=1;
    
    aplus=0.19e-9;
    aplus3  = aplus^3;
    aminus=0.362e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=1.472356*factoreps;
    epsilonmm=0.076923*factoreps;
    epsilonpm=1.438894*factoreps;
        
    sigmapp=(0.221737e-9)/a;
    sigmamm=(0.469906e-9)/a;
    sigmapm=(0.300512e-9)/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        % epsilon of carbon is 0.2330488 kJ/mol
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      % sigma of carbon is 3.400 angstroms
    sigmamw=(sigmamm+sigmaww)/2;
    
    elseif (strcmpi(electrolyte,'KCl'))
    
    % KCl (1:1)
    zplus=1;
    zminus=1;
    
    aplus=0.266e-9;
    aplus3  = aplus^3;
    aminus=0.362e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=1.985740*factoreps;
    epsilonmm=0.076923*factoreps;
    epsilonpm=1.400000*factoreps;
     
    sigmapp=(0.230140e-9)/a;
    sigmamm=(0.469906e-9)/a;
    sigmapm=(0.339700e-9)/a;
    
    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        % epsilon of carbon is 0.2330488 kJ/mol
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      % sigma of carbon is 3.400 angstroms
    sigmamw=(sigmamm+sigmaww)/2;
    
end
    
if (wallLJ==0)
    epsilonpw=0;
    epsilonmw=0;
end

if (ionLJ==0)
    epsilonpp=0;
    epsilonmm=0;
    epsilonpm=0;
end

properties=[zplus,zminus,epsilonpp,epsilonmm,epsilonpm,epsilonpw,epsilonmw,sigmapp,sigmamm,sigmapm,sigmapw,sigmamw,aplus,aminus,aplus3,aminus3];
end

% function to solve the poisson equation with LJ interaction
function differ=odefcn(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,rho,Linf,a)
psi=y(1);
dpsidx=y(2);
awall=1.42e-10;
ap=((aplus3c0/c0)^(1/3));
am=((aminus3c0/c0)^(1/3));
apavg=(ap+awall)/2;
amavg=(am+awall)/2;
apm=(ap+
ULJplus=  -4*pi*rho*epsilonpw*(apm^2)*(((sigmapw^6/(2*(x)^4))-((sigmapw^12/(5*(x)^10))))...
          + ((sigmapw^6/(2*((Linf)-(x))^4))-((sigmapw^12/(5*((Linf)-(x))^10)))));
          
ULJminus= -4*pi*rho*epsilonmw*(apm^2)*(((sigmamw^6/(2*(x)^4))-((sigmamw^12/(5*(x)^10))))...
          + ((sigmamw^6/(2*((Linf)-(x))^4))-((sigmamw^12/(5*((Linf)-(x))^10)))));
if(x==0 || x==Linf)
     ULJplus=inf;
     ULJminus=inf;
end

if (am > ap) 		% Chemical Potential Expressions with anion diameter greater than cation diameter
    fofpsi= (1+(zminus*aplus3c0*(exp(-zplus*psi - ULJplus)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zplus*aminus3c0*(exp(zminus*psi - ULJminus)-fofpsi)) + (zminus*aplus3c0*fofpsi.*(exp(-zplus*psi - ULJplus)-1));
    nplus=fofpsi.*exp(-zplus*psi - ULJplus)./gofpsi;
    nminus=exp(zminus*psi - ULJminus)./gofpsi;
else			% Chemical Potential Expressions with cation diameter greater than anion diameter
    fofpsi= (1+(zplus*aminus3c0*(exp(zminus*psi - ULJminus)-1)./(1-zminus*aplus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zminus*aplus3c0*(exp(-zplus*psi - ULJplus)-fofpsi)) + (zplus*aminus3c0*fofpsi.*(exp(zminus*psi - ULJminus)-1));
    nplus=exp(-zplus*psi - ULJplus)./gofpsi;
    nminus=fofpsi.*exp(zminus*psi - ULJminus)./gofpsi;
end

differ=zeros(2,1);
differ(1)=dpsidx;
differ(2)=(nminus-nplus)./(zplus+zminus); % Poisson Equation
end

% Boundary Conditions
function res=bcfcn(psia,psib,psiD)

res=[ psia(1)-psiD
      psib(1)-psiD 
];
  %    psib(2)+psib(1)/Linf];
end

% function to solve the poisson equation without LJ interaction
function differ=odefcn_0(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,rho,Linf,a)
psi=y(1);
dpsidx=y(2);
awall=1.42e-10;
ap=((aplus3c0/c0)^(1/3));
am=((aminus3c0/c0)^(1/3));
apavg=(ap+awall)/2;
amavg=(am+awall)/2;

if (am > ap) 		% Chemical Potential Expressions with anion diameter greater than cation diameter
    fofpsi= (1+(zminus*aplus3c0*(exp(-zplus*psi)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zplus*aminus3c0*(exp(zminus*psi)-fofpsi)) + (zminus*aplus3c0*fofpsi.*(exp(-zplus*psi)-1));
    nplus=fofpsi.*exp(-zplus*psi)./gofpsi;
    nminus=exp(zminus*psi)./gofpsi;
else			% Chemical Potential Expressions with cation diameter greater than anion diameter
    fofpsi= (1+(zplus*aminus3c0*(exp(zminus*psi)-1)./(1-zminus*aplus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zminus*aplus3c0*(exp(-zplus*psi)-fofpsi)) + (zplus*aminus3c0*fofpsi.*(exp(zminus*psi)-1));
    nplus=exp(-zplus*psi)./gofpsi;
    nminus=fofpsi.*exp(zminus*psi)./gofpsi;
end

differ=zeros(2,1);
differ(1)=dpsidx;
differ(2)=(nminus-nplus)./(zplus+zminus);
end
