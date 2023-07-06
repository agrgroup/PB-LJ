
% Boundary Conditions
function res=bcfcn(psia,psib,psiD)

res=[ psia(1)-psiD
      psib(1)-psiD 
];

end
