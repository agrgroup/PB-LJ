clear
clc
close all

%% Physical constants
kb=1.380649*10^-23;                             % boltzmann constant
permi= 78.38*8.854188*(10^-12);                 % value at T=298.15 and P=0.1 MPa
e=1.602177*10^-19;                              % electronic charge

%% Process parameters
c0 = 5*1000*6.022e23;                           % bulk electrolyte concentration in #/m^3          
psiD=5;                                         % wall potential
electrolyte='NaCl';
T=298.15;                                       % system temperature in K
rho=4/(sqrt(3)*(1.42*sqrt(3)*1e-10)^2);         % particles per unit area of the wall
%% solution parameters 
N=201;                                          % number of grid points
dx=0.1;                                         % this dx is non-dimensional
Linf=(N-1)*dx;                                  % Distance between the plates
%% declaration of variables
filename=[electrolyte, '_c0_', num2str(c0/(1000*6.022e23)), '_N_', num2str(N), '_dx_', num2str(dx),'.mat'];

nplusf=zeros(1,N);                              % initialization of solution array
nminusf=zeros(1,N);
psifinal=zeros(1,N);

%% Electrolyte properties
properties=systemprops(electrolyte, T);
zplus=properties(1);
zminus=properties(2);
epsilonpp=properties(3);                        % LJ parameter for cation-cation interactions
epsilonmm=properties(4);                        % LJ parameter for anion-anion interactions
epsilonpm=properties(5);                        % LJ parameter for cation-anion interactions
epsilonpw=properties(6);                        % LJ parameter for cation-wall interactions
epsilonmw=properties(7);                        % LJ parameter for anion-wall interactions
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
lambda_D = sqrt(permi*kb*T/(e^2*(zplus^2*c0*zminus+zminus^2*c0*zplus)));        % definition of debye length
a=(aplus+aminus)/2;
a3c0=a^3*c0;
four_pi_a3c0 = 4*pi*a3c0;
factorpsi=(e^2*a^2*c0/(kb*T*permi));
aratio=(aplus3/aminus3)^(1/3);

%% Solve
xmesh=linspace(0,Linf*a/lambda_D,N);                                        % definition of the 1-D geometry
guess = @(x) [psiD*exp(-x/lambda_D), -(psiD/lambda_D)*exp(-x/lambda_D)];    % intial guess for the solver
solinit=bvpinit(xmesh,guess);                                               % guess for bvp solver
options= bvpset('stats','off','RelTol',0.001);
sol=bvp5c(@(x,y)odefcn(x,y,zplus,zminus,aplus3c0,aminus3c0,aratio),@(psia,psib) bcfcn(psia,psib,psiD),solinit,options);
y=deval(sol,xmesh);                                                         %obtaining solution in arrays

sol.y=y(1,:);
sol.x=xmesh;

figure;
psi=sol.y(1,:);
% plot(sol.x*lambda_D/a, psi);

xlabel('X');
ylabel({'\bf \psi'});

fofpsi= (1+(zminus*aplus3c0*(exp(-zplus*psi)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);           % equations obtained from chemical potential expressions expression
gofpsi= fofpsi + (zplus*aminus3c0*(exp(zminus*psi)-fofpsi)) + (zminus*aplus3c0*fofpsi.*(exp(-zplus*psi)-1));
nplus=fofpsi.*exp(-zplus*psi)./gofpsi;
nminus=exp(zminus*psi)./gofpsi;

psiInit=sol.y(1,1:N);     
xlistp=linspace(0,Linf,10*N);   
xlistn=linspace(0,Linf,10*N);   
xmeshnew=linspace(0,Linf,N);

nplusold=nplus;
nminusold=nminus;
psiold=psiInit;

steptol=1e-6;                           % non-linear least squares solver options
functol=1e-6;
optol=1e-6;
maxfcneval=500000;
maxiter=80;

lb=zeros(1,3*N);                        % concentrations and electric potential cannot be negative 
ub=[(1/(aplus3c0*zminus))*ones(1,N), (1/(aminus3c0*zplus))*ones(1,N), inf*ones(1,N)];       % upper bound constraint on concentrations to prevent the sum of volume fractions from exceeding unity 

options = optimoptions('lsqnonlin','Display','iter','StepTolerance',steptol,'FunctionTolerance',functol,'MaxIterations',maxiter,'MaxFunctionEvaluations',maxfcneval,'optimalityTolerance',optol);
xfinal=lsqnonlin(@(x) residual(x,aplus3c0,aminus3c0,aplus,aminus,a,four_pi_a3c0,factorpsi,aratio,...
                    N,Linf,dx,rho,psiD,zplus,zminus,epsilonpp,epsilonmm,epsilonpm,epsilonpw,epsilonmw,...
                    sigmapw,sigmamw,sigmamm6,sigmapp6,sigmapm6,sigmamm12,sigmapp12,sigmapm12,xlistp,xlistn)...
                    ,x0,lb,ub,options);
                
nplusf=xfinal(1:N);
nminusf=xfinal((N+1):(2*N));
psifinal=xfinal((2*N)+1:(3*N));

solx=sol.x*lambda_D/a;
beta=psifinal(1,1:N);
gamma=nplusf(1,1:N);
delta=nminusf(1,1:N);
figure;
plot(solx, nplusold,"Color",'green');                           % plotting cation concentration profile with and without LJ interaction             
hold on
plot(xmeshnew,nplusf,"Color",'red');
hold off
xlim([0 (N-1)*dx]);
xlabel('X');
ylabel({'\bf n_+'});
legend('No LJ','Ion LJ');

figure;
plot(solx, nminusold,"Color",'green');                          % plotting anion concentration profile with and without LJ interaction
hold on
plot(xmeshnew, nminusf,"Color",'red');

hold off
xlim([0 (N-1)*dx]);
xlabel('X');
ylabel({'\bf n_-'});
legend('No LJ','Ion LJ');

figure;
plot(solx, psiold,"Color",'green');

hold on
plot(xmeshnew, psifinal,"color","red");                         % plotting electric potential profile with and without LJ interaction
hold off
xlim([0 (N-1)*dx]);
xlabel('X');
ylabel({'\bf \Psi'});
legend('No LJ','Ion LJ');
data=[xmeshnew' ,beta' , gamma' , delta', psiold', nplusold', nminusold'];
writematrix(data,'new_ion_ion_both_c0_5.csv');
%% defining properties for different electrolytes
function properties=systemprops(electrolyte,T)

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
    
properties=[zplus,zminus,epsilonpp,epsilonmm,epsilonpm,epsilonpw,epsilonmw,sigmapp,sigmamm,sigmapm,sigmapw,sigmamw,aplus,aminus,aplus3,aminus3];
end


% calculating the residual which is then minimized using a least-squares solver

function y=residual(x,aplus3c0,aminus3c0,aplus,aminus,a,four_pi_a3c0,factorpsi,aratio,...
                    N,Linf,dx,rho,psiD,zplus,zminus,epsilonpp,epsilonmm,epsilonpm,epsilonpw,epsilonmw,...
                    sigmapw,sigmamw,sigmamm6,sigmapp6,sigmapm6,sigmamm12,sigmapp12,sigmapm12,xlistp,xlistn)

nplus=x(1:N);
nminus=x(N+1:(N*2));
psi=x((2*N+1):(3*N));
y=zeros(1,3*N);                                 
aplusn=aplus/a;                                 % non-dimensionalized cation diameter
awall=1.42e-10/a;                               % non-dimensionalized carbon diameter in the wall
apluswn=(aplusn+awall)/2;                       % average cation-wall diameter
aminusn=aminus/a;                               % non-dimensionalized anion diameter
aminuswn=(aminusn+awall)/2;                     % average anion-wall diameter
aavgn=(aplusn+aminusn)/2;
Np=length(xlistp);
Nn=length(xlistn);
xlist=linspace(0,Linf,N);
sigmapp=sigmapp6^(1/6);
sigmamm=sigmamm6^(1/6);
sigmapm=sigmapm6^(1/6);
for i=1:N
    xcurr=(i-1)*dx;     
    
          
        deltax_curr=xcurr-xlistp;                       % terms used in the denominator of the ion-ion LJ interaction equation
        deltax_Linf=(Linf/2) -xlistp;
        deltax_curr_4= deltax_curr.^4;
        deltax_curr_10=deltax_curr.^10;
        deltax_Linf_4= deltax_Linf.^4;
        deltax_Linf_10=deltax_Linf.^10;
        
        lim1=(xlistp>=0 & xlistp<=xcurr-aplusn) | (xlistp>=(xcurr + aplusn) & xlistp<=Linf);            % cation-cation limits to prevent cations from overlapping 
        lim2=(xlistp>=0 & xlistp<=(Linf/2-aplusn)) | (xlistp>=(Linf/2+aplusn) & xlistp<=Linf);
        lim3=(xlistp>=0 & xlistp<=xcurr-aavgn) | (xlistp>=(xcurr + aavgn) & xlistp<=Linf);              % cation-anion limits to prevent cation-anion from overlap
        lim4=(xlistp>=0 & xlistp<=(Linf/2-aavgn)) | (xlistp>=(Linf/2+aavgn) & xlistp<=Linf);
        
        % integrands corresponding to the 4 integral limits
        integrand1=zeros(Np,1);
        integrand2=integrand1;
        integrand3=integrand2;
        integrand4=integrand3;
        
        % using the trapz function to calculate the integrals
        integrand1_val = interp1(xlist,nplus,xlistp).*((sigmapp6./(2*deltax_curr_4)) - (sigmapp12./(5*deltax_curr_10)));
        integrand1(lim1)=integrand1_val(lim1);
        integral1=-four_pi_a3c0*epsilonpp*zminus*trapz(xlistp,integrand1);
        
        integrand2_val =-interp1(xlist,nplus,xlistp).*((sigmapp6./(2*deltax_Linf_4)) - (sigmapp12./(5*deltax_Linf_10)));
        integrand2(lim2)=integrand2_val(lim2);
        integral2=-four_pi_a3c0*epsilonpp*zminus*trapz(xlistp,integrand2);

        integrand3_val = interp1(xlist,nminus,xlistp).*((sigmapm6./(2*deltax_curr_4)) - (sigmapm12./(5*deltax_curr_10)));
        integrand3(lim3)=integrand3_val(lim3);
        integral3=-0.5*four_pi_a3c0*epsilonpm*zplus*trapz(xlistp,integrand3);       
        
        integrand4_val = -interp1(xlist,nminus,xlistp).*((sigmapm6./(2*deltax_Linf_4)) - (sigmapm12./(5*deltax_Linf_10)));
        integrand4(lim4)=integrand4_val(lim4);
        integral4=-0.5*four_pi_a3c0*epsilonpm*zplus*trapz(xlistp,integrand4);          
               
    
    logtermplus=log((aplus3c0*nplus(i)*zminus)/(1-aplus3c0*nplus(i)*zminus-aminus3c0*nminus(i)*zplus));
    
        y(i)= zplus*psi(i) + logtermplus...
        - log(aplus3c0*zminus/(1-aplus3c0*zminus-aminus3c0*zplus))...
        + integral1 + integral2+ integral3 + integral4;                  % includes electrostatic interaction and LJ interaction                              
    
        
end
%function for mu-
for i=1:N
    xcurr=(i-1)*dx;     

        deltax_curr=xcurr-xlistn;                                   % terms used in the denominator of the ion-ion LJ interaction equation
        deltax_Linf=(Linf/2) -xlistn;
        deltax_curr_4= deltax_curr.^4;
        deltax_curr_10=deltax_curr.^10;
        deltax_Linf_4= deltax_Linf.^4;
        deltax_Linf_10=deltax_Linf.^10;
               
        lim1=(xlistn>=0 & xlistn<=xcurr-aminusn) | (xlistn>=(xcurr + aminusn) & xlistn<=Linf);          % anion-anion limits to prevent anion-anion overlap
        lim2=(xlistn>=0 & xlistn<=(Linf/2-aminusn)) | (xlistn>=(Linf/2+aminusn) & xlistn<=Linf);
        lim3=(xlistn>=0 & xlistn<=xcurr-aavgn) | (xlistn>=(xcurr + aavgn) & xlistn<=Linf);              % anion-cation limits to prevent anion-cation overlap
        lim4=(xlistn>=0 & xlistn<=(Linf/2-aavgn)) | (xlistn>=(Linf/2+aavgn) & xlistn<=Linf);

        % integrands corresponding to the 4 integral limits
        integrand1=zeros(Nn,1);
        integrand2=integrand1;
        integrand3=integrand2;
        integrand4=integrand3;
        
        % using the trapz function to calculate the integrals
        integrand1_val = interp1(xlist,nminus,xlistn).*((sigmamm6./(2*deltax_curr_4)) - (sigmamm12./(5*deltax_curr_10)));
        integrand1(lim1)=integrand1_val(lim1);
        integral1=-four_pi_a3c0*epsilonmm*zplus*trapz(xlistn,integrand1);
        
        integrand2_val =-interp1(xlist,nminus,xlistn).*((sigmamm6./(2*deltax_Linf_4)) - (sigmamm12./(5*deltax_Linf_10)));
        integrand2(lim2)=integrand2_val(lim2);
        integral2=-four_pi_a3c0*epsilonmm*zplus*trapz(xlistn,integrand2);
      
        integrand3_val = interp1(xlist,nplus,xlistn).*((sigmapm6./(2*deltax_curr_4)) - (sigmapm12./(5*deltax_curr_10)));
        integrand3(lim3)=integrand3_val(lim3);
        integral3=-0.5*four_pi_a3c0*epsilonpm*zminus*trapz(xlistn,integrand3);       
        
        integrand4_val = -interp1(xlist,nplus,xlistn).*((sigmapm6./(2*deltax_Linf_4)) - (sigmapm12./(5*deltax_Linf_10)));
        integrand4(lim4)=integrand4_val(lim4);
        integral4=-0.5*four_pi_a3c0*epsilonpm*zminus*trapz(xlistn,integrand4);    

    logtermminus=log(aminus3c0*nminus(i)*zplus/(1-aminus3c0*nminus(i)*zplus))-((aminus3c0/aplus3c0)*log((1-aminus3c0*nminus(i)*zplus-aplus3c0*nplus(i)*zminus)/(1-aminus3c0*nminus(i)*zplus)));
    
    y(N+i)= -zminus*psi(i) + logtermminus...
        - log(aminus3c0*zplus/(1-aminus3c0*zplus))-((aminus3c0/aplus3c0)*log((1-aminus3c0*zplus-aplus3c0*zminus)/(1-aminus3c0*zplus)))...
        + integral1 + integral2+ integral3 + integral4;                             %includes electrostatic interaction and LJ interaction          
end
y((N*2)+1)= psi(1)-psiD;            % potential is psiD at the left wall (boundary condition)

%solving the poisson equation
for i=2:(N-1)
    y((N*2)+i)= ((psi(i+1)- 2*psi(i) + psi(i-1))/(dx^2)) + factorpsi*(zplus*zminus*nplus(i) - zminus*zplus*nminus(i));
end
y(3*N)=psi(N)-psiD;                 % potential is psiD at the right wall (boundary condition)

end

% function to solve the poisson equation without LJ interaction
function differ=odefcn(x,y,zplus,zminus,aplus3c0,aminus3c0,aratio)
psi=y(1);
dpsidx=y(2);

fofpsi= (1+(zminus*aplus3c0*(exp(-zplus*psi)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
gofpsi= fofpsi + (zplus*aminus3c0*(exp(zminus*psi)-fofpsi)) + (zminus*aplus3c0*fofpsi.*(exp(-zplus*psi)-1));
nplus=fofpsi.*exp(-zplus*psi)./gofpsi;
nminus=exp(zminus*psi)./gofpsi;

differ=zeros(2,1);
differ(1)=dpsidx;
differ(2)=(nminus-nplus)./(zplus+zminus);           % poisson equation
end

% boundary condition function for function named odefcn
function res=bcfcn(psia,psib,psiD)

res=[ psia(1)-psiD
      psib(1)-psiD 
];
  %    psib(2)+psib(1)/Linf];
end