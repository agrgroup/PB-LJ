clear
clc
close all

%% Physical constants
kb=1.380649*10^-23;			% boltzmann constant
permi= 78.38*8.854188*(10^-12);	%value at T=298.15 and P=0.1 MPa
e=1.602177*10^-19;			% electronic charge
wallLJ=1;
ionLJ=0;

%% Process parameters
c0 = 1*1000*6.022e23;  		% bulk electrolyte concentration in #/m^3
psiD=5;				% wall potential
electrolyte='hydrated_NaCl';
T=298.15;					% system temperature in K
rho=4/(sqrt(3)*(1.42*sqrt(3)*1e-10)^2);	% particles per unit area of the wall

%% solution parameters 
N=201;			% number of grid points
dx=0.1;		% dx is non-dimensional
Linf=(N-1)*dx;		% distance between the plates

%% initialization of solution arrays
nplusf=zeros(1,N);
nminusf=zeros(1,N);
psifinal=zeros(1,N);

%% Electrolyte properties
properties=systemprops(electrolyte,wallLJ,ionLJ,T);
zplus=properties(1);			% cation valence
zminus=properties(2);			% anion valence
epsilonpp=properties(3);		% LJ parameter epsilon for cation-cation interactions
epsilonmm=properties(4);		% LJ parameter epsilon for anion-anion interactions
epsilonpm=properties(5);		% LJ parameter epsilon for cation-anion interactions
epsilonpw=properties(6);		% LJ parameter epsilon for cation-wall interactions
epsilonmw=properties(7);		% LJ parameter epsilon for anion-wall interactions
sigmapp=properties(8);			% LJ parameter sigma for cation-cation interactions
sigmamm=properties(9);			% LJ parameter sigma for anion-anion interactions
sigmapm=properties(10);		% LJ parameter sigma for cation-anion interactions
sigmapw=properties(11);		% LJ parameter sigma for cation-wall interactions
sigmamw=properties(12);		% LJ parameter sigma for anion-wall interactions
aplus=properties(13);			% cation diameter
aminus=properties(14);			% anion diameter
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
lambda_D = sqrt(permi*kb*T/(e^2*(zplus^2*c0*zminus+zminus^2*c0*zplus))); % debye length
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
xmesh=linspace(0,Linf,N);							% definition of the 1-D geometry
guess = @(x) [psiD*exp(-x/lambda_D), -(psiD/lambda_D)*exp(-x/lambda_D)];	% intial guess for the solver
solinit=bvpinit(xmesh,guess);							% guess for bvp solver
options= bvpset('stats','off','RelTol',0.001,'Stats','on');
sol=bvp5c(@(x,y)odefcn(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,rho,Linf,a),@(psia,psib) bcfcn(psia,psib,psiD),solinit,options);
y=deval(sol,xmesh);								% obtaining solution in arrays

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


ULJplus_lw=-4*pi*rho*epsilonpw.*((((sigmapw.^6)./(2.*(xmesh/a).^4))-((sigmapw.^12)./(5.*(xmesh/a).^10))));
ULJplus_rw=-4*pi*rho*epsilonpw.*(apavg^2)*(((sigmapw.^6./(2.*((Linf/a)-(xmesh/a)).^4))-((sigmapw.^12./(5.*((Linf/a)-(xmesh/a)).^10)))));

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
    fofpsi= (1+(zminus*aplus3c0*(exp(-zplus*psi_zero)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zplus*aminus3c0*(exp(zminus*psi_zero)-fofpsi)) + (zminus*aplus3c0*fofpsi.*(exp(-zplus*psi_zero)-1));
    nplus_zero=fofpsi.*exp(-zplus*psi_zero)./gofpsi;
    nminus_zero=exp(zminus*psi_zero)./gofpsi;
else			% Chemical Potential Expressions with cation diameter greater than anion diameter
    fofpsi= (1+(zplus*aminus3c0*(exp(zminus*psi_zero)-1)./(1-zminus*aplus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zminus*aplus3c0*(exp(-zplus*psi_zero)-fofpsi)) + (zplus*aminus3c0*fofpsi.*(exp(zminus*psi_zero)-1));
    nplus_zero=exp(-zplus*psi_zero)./gofpsi;
    nminus_zero=fofpsi.*exp(zminus*psi_zero)./gofpsi;
end

% saving the double layer properties (electric potential, cation concentration, anion concentration)
final=[sol.x', nplus', nplus_zero', nminus', nminus_zero', psi', psi_zero'];
writematrix(final,"solution.csv");

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
