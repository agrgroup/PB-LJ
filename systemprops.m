
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%																								                                                                                         
%% PARAMETER REFERENCES																					                                                                                 
%% ionic diameters reference: Israelachvili, Intermolecular and surface forces, Academic Press, 2011.; Nightingale, J. Phys. Chem., 1959, 63, 1381–1387.					             
%% inter-ionic LJ interaction parameters reference: Zeron et al., J. Chem. Phys., 2019, 151, 134504.; Xantheas and Dany, J. Phys. Chem., 1996, 100, 3989–3995				             
%% LJ parameters for graphene (wall) reference:  Cheng and Steele, J. Chem. Phys., 1990, 92, 3858.												                                         
%%																								                                                                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function properties=systemprops(electrolyte,wallLJ,ionLJ,T)

factoreps=1000/(8.314*T);

if (strcmpi(electrolyte,'MgCl2'))
    % MgCl2 (2:1)
    zplus=2;					% cation valence
    zminus=1;					% anion valence
    
    aplus=0.14e-9;				% cation diameter
    aplus3  = aplus^3;
    aminus=0.360e-9;				% anion diameter
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=3.651900*factoreps;		% LJ parameter epsilon for cation-cation interactions       
    epsilonmm=0.076923*factoreps;		% LJ parameter epsilon for anion-anion interactions       
    epsilonpm=3.000*factoreps;	
      
    sigmapp=(0.116290e-9)/a;			% LJ parameter sigma for cation-cation interactions         
    sigmamm=(0.469906e-9)/a;			% LJ parameter sigma for anion-anion interactions         
    sigmapm=0.300e-9/a;

    %%Lorentz-Bethelot Rules
    epsilonww = 0.2330488*factoreps;		% LJ parameter epsilon for wall-wall interactions
    epsilonpw=sqrt(epsilonpp*epsilonww);	% LJ parameter epsilon for cation-wall interactions        
    epsilonmw=sqrt(epsilonmm*epsilonww);	% LJ parameter epsilon for anion-wall interactions 
    
    sigmaww=0.34e-9/a;				% LJ parameter sigma for wall-wall interactions
    sigmapw=(sigmapp+sigmaww)/2;		% LJ parameter sigma for cation-wall interactions      
    sigmamw=(sigmamm+sigmaww)/2;		% LJ parameter sigma for anion-wall interactions

elseif (strcmpi(electrolyte,'hydrated_MgCl2'))
    % MgCl2 (2:1)
    zplus=2;
    zminus=1;
    
    aplus=0.86e-9;
    aplus3  = aplus^3;
    aminus=0.66e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=3.651900*factoreps;       
    epsilonmm=0.076923*factoreps;       
    epsilonpm=3.000*factoreps;
      
    sigmapp=(0.116290e-9)/a;        
    sigmamm=(0.469906e-9)/a;         
    sigmapm=0.300e-9/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);       
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      
    sigmamw=(sigmamm+sigmaww)/2;    
    
elseif (strcmpi(electrolyte,'NaCl'))
    
    % NaCl (1:1)
    zplus=1;
    zminus=1;
    
    aplus=0.194e-9;
    aplus3  = aplus^3;
    aminus=0.360e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=1.472356*factoreps;
    epsilonmm=0.076923*factoreps;
    epsilonpm=1.438894*factoreps;
        
    sigmapp=(0.221737e-9)/a;
    sigmamm=(0.469906e-9)/a;
    sigmapm=(0.300512e-9)/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      
    sigmamw=(sigmamm+sigmaww)/2;

elseif (strcmpi(electrolyte,'hydrated_NaCl'))
    
    % NaCl (1:1)
    zplus=1;
    zminus=1;
    
    aplus=0.72e-9;
    aplus3  = aplus^3;
    aminus=0.66e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=1.472356*factoreps;
    epsilonmm=0.076923*factoreps;
    epsilonpm=1.438894*factoreps;
        
    sigmapp=(0.221737e-9)/a;
    sigmamm=(0.469906e-9)/a;
    sigmapm=(0.300512e-9)/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      
    sigmamw=(sigmamm+sigmaww)/2;    
    
elseif (strcmpi(electrolyte,'KCl'))
    
    % KCl (1:1)
    zplus=1;
    zminus=1;
    
    aplus=0.282e-9;
    aplus3  = aplus^3;
    aminus=0.360e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=1.985740*factoreps;
    epsilonmm=0.076923*factoreps;
    epsilonpm=1.400000*factoreps;
     
    sigmapp=(0.230140e-9)/a;
    sigmamm=(0.469906e-9)/a;
    sigmapm=(0.339700e-9)/a;
    
    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);       
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;     
    sigmamw=(sigmamm+sigmaww)/2;

elseif (strcmpi(electrolyte,'hydrated_KCl'))
    
    % KCl (1:1)
    zplus=1;
    zminus=1;
    
    aplus=0.532e-9;
    aplus3  = aplus^3;
    aminus=0.532e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=1.985740*factoreps;
    epsilonmm=0.076923*factoreps;
    epsilonpm=1.400000*factoreps;
     
    sigmapp=(0.230140e-9)/a;
    sigmamm=(0.469906e-9)/a;
    sigmapm=(0.339700e-9)/a;
    
    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      
    sigmamw=(sigmamm+sigmaww)/2;
    
elseif (strcmpi(electrolyte,'Na2SO4'))
    
    % Na2SO4 (1:2)
    zplus=1;
    zminus=2;
    
    aplus=0.194e-9;
    aplus3  = aplus^3;
    aminus=0.484e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=1.472356*factoreps;
    epsilonmm=2.8947*factoreps;
    epsilonpm=1.438894*factoreps;
        
    sigmapp=(0.221737e-9)/a;
    sigmamm=(0.508620e-9)/a;
    sigmapm=(0.300512e-9)/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;     
    sigmamw=(sigmamm+sigmaww)/2;
    
elseif (strcmpi(electrolyte,'hydrated_Na2SO4'))
    
    % Na2SO4 (1:2)
    zplus=1;
    zminus=2;
    
    aplus=0.72e-9;
    aplus3  = aplus^3;
    aminus=0.76e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=1.472356*factoreps;
    epsilonmm=2.8947*factoreps;
    epsilonpm=1.438894*factoreps;
        
    sigmapp=(0.221737e-9)/a;
    sigmamm=(0.508620e-9)/a;
    sigmapm=(0.300512e-9)/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);       
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;     
    sigmamw=(sigmamm+sigmaww)/2;

elseif (strcmpi(electrolyte,'MgSO4'))
    zplus=2;
    zminus=2;
    
    aplus=0.14e-9;
    aplus3  = aplus^3;
    aminus=0.484e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=3.6519*factoreps;       
    epsilonmm=2.8947*factoreps;       
    epsilonpm=0.000*factoreps;
      
    sigmapp=(0.116290e-9)/a;         
    sigmamm=(0.508620e-9)/a;
    sigmapm=0.000e-9/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      
    sigmamw=(sigmamm+sigmaww)/2;

elseif (strcmpi(electrolyte,'hydrated_MgSO4'))
    zplus=2;
    zminus=2;
    
    aplus=0.86e-9;
    aplus3  = aplus^3;
    aminus=0.66e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=3.6519*factoreps;       
    epsilonmm=2.8947*factoreps;       
    epsilonpm=0.000*factoreps;
      
    sigmapp=(0.116290e-9)/a;         
    sigmamm=(0.508620e-9)/a;
    sigmapm=0.000e-9/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;     
    sigmamw=(sigmamm+sigmaww)/2;   
    
elseif (strcmpi(electrolyte,'CaSO4'))
    zplus=2;
    zminus=2;
    
    aplus=0.206e-9;
    aplus3  = aplus^3;
    aminus=0.484e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=0.507200*factoreps;       
    epsilonmm=2.8947*factoreps;      
    epsilonpm=0.000*factoreps;
      
    sigmapp=(0.266560e-9)/a;        
    sigmamm=(0.508620e-9)/a;
    sigmapm=0.000e-9/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);        
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      
    sigmamw=(sigmamm+sigmaww)/2;

elseif (strcmpi(electrolyte,'hy_CaSO4'))
    zplus=2;
    zminus=2;
    
    aplus=0.4844e-9;
    aplus3  = aplus^3;
    aminus=0.763e-9;
    aminus3 = aminus^3; 
    a=(aplus+aminus)/2;
    
    epsilonpp=0.507200*factoreps;      
    epsilonmm=2.8947*factoreps;       
    epsilonpm=0.000*factoreps;
      
    sigmapp=(0.266560e-9)/a;        
    sigmamm=(0.508620e-9)/a;
    sigmapm=0.000e-9/a;

    epsilonww = 0.2330488*factoreps;
    epsilonpw=sqrt(epsilonpp*epsilonww);       
    epsilonmw=sqrt(epsilonmm*epsilonww);
    
    sigmaww=0.34e-9/a;
    sigmapw=(sigmapp+sigmaww)/2;      
    sigmamw=(sigmamm+sigmaww)/2;    
   
elseif (strcmpi(electrolyte,'NaCl_aluru_2'))

    zplus=1;
    zminus=1;
    aplus=0.21595e-9;
    aplus3=aplus^3;
    aminus=0.48304e-9;
    aminus3 = aminus^3;
    a=(aplus+aminus)/2;
    epsilonpp=1.47545*factoreps;
    epsilonmm=0.05349*factoreps;
    epsilonpm=0;
    sigmapp=(0.21595e-9)/a;
    sigmamm=(0.48304e-9)/a;
    sigmapm=0;
    epsilonww = 0.2334*factoreps;
    epsilonpw=0.2328*factoreps;
    epsilonmw=0.1781*factoreps;
    sigmaww=0.339e-9/a;
    sigmapw=(0.4596e-9)/a; 
    sigmamw=(0.3814e-9)/a;
    
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
