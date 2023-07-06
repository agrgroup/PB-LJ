
% function to solve the poisson equation with LJ interaction
function differ=odefcn(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,rho,Linf,a)
psi=y(1);
dpsidx=y(2);
ap=((aplus3c0/c0)^(1/3));
am=((aminus3c0/c0)^(1/3));
apm=(ap+am)/2;

ULJplus=  -4*pi*rho*epsilonpw*(apm^2)*(((sigmapw^6/(2*(x)^4))-((sigmapw^12/(5*(x)^10))))...
          + ((sigmapw^6/(2*((Linf)-(x))^4))-((sigmapw^12/(5*((Linf)-(x))^10)))));
          
ULJminus= -4*pi*rho*epsilonmw*(apm^2)*(((sigmamw^6/(2*(x)^4))-((sigmamw^12/(5*(x)^10))))...
          + ((sigmamw^6/(2*((Linf)-(x))^4))-((sigmamw^12/(5*((Linf)-(x))^10)))));
if(x==0 || x==Linf)
     ULJplus=inf;
     ULJminus=inf;
end    
if (am > ap) 					% Chemical Potential Expressions with anion diameter greater than cation diameter
    fofpsi= (1+(zminus*aplus3c0*(exp(-zplus*psi - ULJplus)-1)./(1-zplus*aminus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zplus*aminus3c0*(exp(zminus*psi - ULJminus)-fofpsi)) + (zminus*aplus3c0*fofpsi.*(exp(-zplus*psi - ULJplus)-1));
    nplus=fofpsi.*exp(-zplus*psi - ULJplus)./gofpsi;
    nminus=exp(zminus*psi - ULJminus)./gofpsi;
else						% Chemical Potential Expressions with cation diameter greater than anion diameter
    fofpsi= (1+(zplus*aminus3c0*(exp(zminus*psi - ULJminus)-1)./(1-zminus*aplus3c0))).^(aratio^3-1);
    gofpsi= fofpsi + (zminus*aplus3c0*(exp(-zplus*psi - ULJplus)-fofpsi)) + (zplus*aminus3c0*fofpsi.*(exp(zminus*psi - ULJminus)-1));
    nplus=exp(-zplus*psi - ULJplus)./gofpsi;
    nminus=fofpsi.*exp(zminus*psi - ULJminus)./gofpsi;
end
differ=zeros(2,1);
differ(1)=dpsidx;
differ(2)=(nminus-nplus)./(zplus+zminus);	% Poisson Equation
end

