function differ=odefcn_0(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,rho,Linf,a)
psi=y(1);
dpsidx=y(2);
ap=((aplus3c0/c0)^(1/3));
am=((aminus3c0/c0)^(1/3));

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
