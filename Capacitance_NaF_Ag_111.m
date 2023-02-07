clear all;
close all;

% UFF LJ with Israelesvili radii
params =[0.47088 0.251093 0.25101 0.30816 0.72 0.70];

% UFF LJ with Marcus radii
%params =[0.47088 0.251093 0.25101 0.30816 0.4712 0.5260];

% Interface FF with Israelesvili radii
%params =[5.30010 0.24249 2.82536 0.299565 0.72 0.70];

%Interface FF with Marcus radii
%params =[5.30010 0.24249 2.82536 0.299565 0.4712 0.5260];

    
T         = 298.15;			% system temperature in K
factoreps=1000/(8.314*T);

aplus     = params(5)*1e-9;		% cation diameter
aminus    = params(6)*1e-9;		% anion diameter
a=(aplus+aminus)/2;

epsilonpw = params(1)*factoreps;	% LJ parameter for cation-wall interactions
sigmapw   = params(2)*1e-9/a;		% LJ parameter sigma for cation-wall interactions
epsilonmw = params(3)*factoreps;	% LJ parameter for anion-wall interactions
sigmamw   = params(4)*1e-9/a;		% LJ parameter sigma for anion-wall interactions

zplus  =  0.85;			% cation valence
zminus =  0.85;			% anion valence

rho=4/(sqrt(3)*(4.0853*1e-10)^2);	%particles per unit area of the wall


%% Physical constants
kb=1.380649*10^-23;			% boltzmann constant
permi= 78.38*8.854188*(10^-12);	% dielectric constant value at T=298.15 and P=0.1 MPa
e=1.602177*10^-19;			% charge of an electron


aplus3  = aplus^3;
aminus3 = aminus^3;


% Experimental Capacitance Data for 1 M NaF on Ag(111) electrode from Valette, J. Electroanal. Chem., 1989, 269, 191–203
psiD_range_point1=[-0.745711624	-0.716940217	-0.688167095	-0.659389689	-0.630613139	-0.601832947	-0.573052969	-0.544269135	-0.515482943	-0.486696323	-0.457904775	-0.429107013	-0.400313537	-0.371522631	-0.342730869	-0.313925758	-0.285118418	-0.256298801	-0.227469756	-0.20125535	-0.176350461	-0.151584255	-0.12689007	-0.085930863	-0.05715667	-0.028399619	0.000362146	0.02915798	0.051444809	0.065896135	0.077717644	0.089545988	0.102687728	0.114511594	0.126338995	0.147322181	0.176124906	0.204905063	0.233671006	0.262422646	0.291167075	0.320954266	0.348662893	0.37741076	0.406162255	0.43490255	0.463650388	0.492400583	0.521152921	0.549910186	0.578669102	0.60480396	0.649241917	0.677992754	0.706747877	0.735510821	0.764269836	0.783873468];
C_expt_point1 = [21.52558126	21.09724235	20.76472314	20.67175318	20.53087337	20.59361042	20.64437	20.91072391	21.3088299	21.73089082	22.42843337	23.47332232	24.27866203	24.9402722	25.64979221	27.10550812	28.68578965	30.95237976	33.74597821	36.63126579	39.43956241	43.09841401	45.46862126	48.27011898	47.99748708	46.7666582	45.79933349	46.73642529	49.44517088	52.89929889	55.59664589	58.67607395	62.04935158	64.87845068	67.9051779	70.86524808	72.18761492	72.24835572	71.51459152	69.98124617	68.04483068	65.89884583	64.56109614	62.81683337	61.27538897	59.10792253	57.36206277	55.7479551	54.25362205	53.03477063	51.90814567	50.34673383	48.17849889	46.6003236	45.26169756	44.36024889	43.23921342	42.20719291];

% Experimental Capacitance Data for 0.04 M NaF on Ag(111) electrode from Valette, J. Electroanal. Chem., 1989, 269, 191–203
psiD_range_point04=[-0.74231629	-0.71305287	-0.683790316	-0.654524078	-0.625259791	-0.595992036	-0.566723198	-0.537452842	-0.508180319	-0.478904112	-0.449621186	-0.420336744	-0.390386194	-0.362205854	-0.334548805	-0.302744095	-0.283393379	-0.258524257	-0.232065677	-0.241901684	-0.211955811	-0.208125391	-0.193126919	-0.168337657	-0.137926908	-0.108633395	-0.080714362	-0.052826637	-0.028946199	-0.00505157	0.028194672	0.050863989	0.069565877	0.081607149	0.096312763	0.104352356	0.115053272	0.125758797	0.148489273	0.179147961	0.208431811	0.237685505	0.266924692	0.296163929	0.325398152	0.354632433	0.383867174	0.413104388	0.442349632	0.471594971	0.500835852	0.530079594	0.559321241	0.588297333	0.617808139	0.647049644	0.676299832	0.705543579	0.734786295	0.762927525	0.782658139];
C_expt_point04 = [21.30075736	21.12160086	20.89350764	20.87339548	20.74317569	20.80870278	20.93540077	21.14773801	21.48241705	22.02507714	22.94699679	23.95455569	25.01615521	26.27658537	27.20107971	28.63081648	30.30717205	32.38791071	35.24129856	33.8274086	37.51327055	38.08182062	40.81612355	43.75274608	45.76075745	47.2803459	46.30595509	43.56419142	39.89736933	37.03159472	36.3964908	39.45728961	43.84301213	47.71402029	51.79629723	55.05882094	58.36215053	61.92566033	68.43889677	71.92720747	72.90131874	72.17308061	70.62597778	69.08164221	67.25431842	65.43024736	63.63210111	61.97359929	60.76840294	59.56856632	58.11705516	56.82703516	55.41875143	54.0390116	52.80567916	51.38941409	50.46327354	49.1736031	47.82567465	46.87838321	46.75612456];

% Experimental Capacitance Data for 0.02 M NaF on Ag(111) electrode from Valette, J. Electroanal. Chem., 1989, 269, 191–203
psiD_range_point02=[-0.742314773	-0.71305287	-0.683790316	-0.654524078	-0.625259791	-0.595992036	-0.566723198	-0.537452842	-0.508180103	-0.478900861	-0.449619886	-0.42033566	-0.390604974	-0.363187807	-0.333594703	-0.311142108	-0.279353671	-0.289166001	-0.253225781	-0.223954169	-0.202861669	-0.207198812	-0.176840428	-0.145905857	-0.11663602	-0.088592144	-0.067510341	-0.050272367	-0.03170657	-0.01154541	0.018744381	0.04271362	0.05739599	0.070761804	0.081459938	0.092160456	0.101525744	0.109559775	0.118938274	0.128306641	0.144883562	0.155038442	0.163063149	0.17452428	0.195289583	0.219073019	0.248320314	0.277559585	0.306799743	0.336030438	0.365262059	0.394274991	0.423734515	0.452976132	0.480644403	0.506965238	0.530051086	0.554765112	0.576652957	0.604703223	0.631267367	0.658880298	0.681983352	0.713516666	0.74275562	0.772009766	0.787954597];
C_expt_point02=[21.38639662	21.12160086	20.89350764	20.87339548	20.74317569	20.80870278	20.93540077	21.14773801	21.49465123	22.20858983	23.02040186	24.01572659	25.18190924	26.68554874	28.61254766	30.16731862	33.06241444	31.64409485	35.39636294	38.54042282	41.09522715	40.32363059	44.18942516	45.92159019	46.10466178	43.83371547	40.71873176	37.56189523	34.26525215	30.91234702	28.59548168	29.94146898	32.71163018	36.25762961	39.4039539	42.68485416	45.68741419	48.63592729	52.38426252	55.56064986	60.14459871	62.6623814	65.08460326	68.3904351	71.24883261	72.83413189	71.74469418	70.20230054	68.71000878	66.68349118	64.70926014	62.90591549	61.28089632	59.87098137	57.99316786	56.80457005	55.21778603	54.4331749	53.1363305	52.43425481	51.6121747	50.19792823	49.46048288	48.68180668	47.12152784	46.41886497	45.35953674];

% Experimental Capacitance Data for 0.01 M NaF on Ag(111) electrode from Valette, J. Electroanal. Chem., 1989, 269, 191–203
psiD_range_point01=[-0.742314773	-0.71305287	-0.683790316	-0.654524078	-0.625259791	-0.595992036	-0.566723198	-0.537452409	-0.508178369	-0.478899344	-0.449618369	-0.420330459	-0.391038387	-0.360313811	-0.332033593	-0.308641658	-0.27953106	-0.258352041	-0.234091067	-0.224125881	-0.199527021	-0.169873129	-0.143093288	-0.114044278	-0.094143828	-0.078236686	-0.062325967	-0.043762271	-0.02253447	0.005341182	0.034600918	0.057250765	0.073276445	0.086651662	0.097353373	0.108061043	0.118767918	0.12813261	0.134826046	0.147403	0.154911719	0.162953101	0.173649248	0.186124644	0.201703591	0.229688884	0.258951886	0.288190449	0.317429418	0.346665989	0.375898449	0.405136842	0.434372381	0.463611751	0.492859104	0.522100526	0.551343092	0.576445186	0.60450721	0.631269811	0.653466428	0.680497409	0.710851148	0.740988232	0.769568602	0.786631826];
C_expt_point01=[21.38639662	21.12160086	20.89350764	20.87339548	20.74317569	20.80870278	20.93540077	21.17220637	21.59252466	22.29422908	23.10604112	24.3093469	25.74754892	27.17861003	29.26341367	31.13960201	33.77671753	36.38798377	39.18332856	40.55047724	42.84912721	44.68668147	44.48758231	42.21765594	39.16677715	35.97992511	32.99493704	29.57973893	26.36235202	22.93908629	22.55194874	24.51371043	28.01827487	32.09510176	35.44329001	39.12791819	42.76768772	45.73660376	48.10398069	52.17604226	55.50890744	58.87236313	61.90654077	65.25202717	68.570356	71.33636166	71.13355241	69.55122291	67.99178881	66.29699225	64.37009583	62.77818289	61.0250992	59.48828862	58.4021522	56.98117794	55.62479238	54.16893726	52.9224574	51.75016826	50.65567993	49.61436078	48.40416694	47.47968115	46.29007166	45.78451619];

% Experimental Capacitance Data for 0.005 M NaF on Ag(111) electrode from Valette, J. Electroanal. Chem., 1989, 269, 191–203
psiD_range_point005=[-0.74231564	-0.71305287	-0.683790316	-0.654524078	-0.625259791	-0.595992036	-0.566723198	-0.537452842	-0.508180319	-0.478902525	-0.449616424	-0.420325556	-0.391026956	-0.361722936	-0.327082955	-0.301401801	-0.280177335	-0.259115722	-0.231343195	-0.201821729	-0.170269043	-0.15158622	-0.131380885	-0.114143535	-0.096906072	-0.079675421	-0.062435763	-0.046525707	-0.023976637	0.005232205	0.035561187	0.062474872	0.082483945	0.098513595	0.110550698	0.12125521	0.130617518	0.139983084	0.150689741	0.160047817	0.170754899	0.181466336	0.192173434	0.205553648	0.220243845	0.239864101	0.271320603	0.30148068	0.330720075	0.359953083	0.389185217	0.418419889	0.447662421	0.476909834	0.506149485	0.535394765	0.564636516	0.594150347	0.623130802	0.652374894	0.681626765	0.709547226	0.7401273	0.7693788	0.785334598];
C_expt_point005=[21.3374599	21.12160086	20.89350764	20.87339548	20.74317569	20.80870278	20.93540077	21.14773801	21.48241705	22.11468045	23.21580674	24.58608497	26.3928296	28.50554272	31.45151983	33.85913952	35.91497589	38.33544588	41.60663904	43.47830138	43.7931694	42.33465664	39.80667856	36.61459594	33.42892171	29.85874469	26.79696579	23.7745955	20.04746375	16.78733306	16.68281516	19.03284861	22.11372114	25.84240799	29.47813248	32.9844475	35.81878757	38.83704814	42.46448153	45.05994924	48.71138202	52.60864024	56.26097019	60.61980851	63.83182749	66.41983375	68.45705689	68.83055226	67.29514958	65.39924642	63.45394047	61.65194919	60.29365083	59.21091279	57.68992896	56.48677164	55.0844011	54.02178923	52.88835216	51.61808162	50.78698145	49.89318555	48.94369619	48.09165982	47.65138107];


color=['rgbk'];
color2=['bk'];
color3=['rg'];
for j=1:3
    
    if (j==1)
        psiD_range = psiD_range_point1;
        C_expt = C_expt_point1;
        c0_mol_L = 0.1;
    elseif (j==2)
        psiD_range = psiD_range_point005;
        C_expt = C_expt_point005;
        c0_mol_L = 0.005;
    elseif (j==3)
        psiD_range = psiD_range_point04;
        C_expt = C_expt_point04;
        c0_mol_L = 0.04;
    end
    
    c0 = c0_mol_L*1000*6.022e23;						% electrolyte concentration
    aplus3c0  =aplus3*c0;
    aminus3c0 =aminus3*c0;

    Nexpt=length(psiD_range);
    psiD_range=psiD_range./((kb*T)/e);					% range of electrode potential
    psiD_range=psiD_range(6:Nexpt-6);
    C_expt=C_expt(6:Nexpt-6);
    
    wall_on_both_sides = 1;
    
    
    %% solution parameters
    N=500;
    dx=0.1;                     						% this dx is non-dimensional
    Linf=(N-1)*dx;
    
    aplus3  = aplus^3;
    aminus3 = aminus^3;
    
    lambda_D = sqrt(permi*kb*T/(e^2*(zplus^2*c0*zminus+zminus^2*c0*zplus)));	% definition of debye length
    
    
    if (aminus > aplus)
        aratio=(aminus3/aplus3)^(1/3);
    else
        aratio=(aplus3/aminus3)^(1/3);
    end
    
    Npsi=length(psiD_range);
    
    C_vdw_hy=zeros(1,Npsi);
    for i = 1:Npsi
        xmesh=linspace(0,Linf,N);						% definition of the 1-D geometry
        
        % computing electrode surface charge with pot + 0.01
        psiD_fwd=psiD_range(i)+0.01;
        guess = @(x) [psiD_fwd*exp(-x), -psiD_fwd*exp(-x)];			% intial guess for the solver
        solinit=bvpinit(xmesh,guess);						% guess for bvp solver
        options= bvpset('stats','off','RelTol',0.001);
        sol=bvp5c(@(x,y)odefcn(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,Linf,a,rho,wall_on_both_sides),@(psia,psib) bcfcn(psia,psib,psiD_fwd,wall_on_both_sides),solinit,options);
        y=deval(sol,xmesh);							% obtaining solution in arrays
        sol.y=y(1,:);
        sol.x=xmesh;
        dpsi=y(2,:);
        q_fwd=-dpsi(1);
        
        % computing electrode surface charge with pot - 0.01
        psiD_bkd=psiD_range(i)-0.01;
        guess = @(x) [psiD_bkd*exp(-x), -psiD_bkd*exp(-x)];			% intial guess for the solver
        solinit=bvpinit(xmesh,guess);						% guess for bvp solver
        options= bvpset('stats','off','RelTol',0.001);
        sol=bvp5c(@(x,y)odefcn(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,Linf,a,rho,wall_on_both_sides),@(psia,psib) bcfcn(psia,psib,psiD_bkd,wall_on_both_sides),solinit,options);
        y=deval(sol,xmesh);							% obtaining solution in arrays
        sol.y=y(1,:);
        sol.x=xmesh;
        dpsi=y(2,:);
        q_bkd=-dpsi(1);
        
        C_conv = ((permi/lambda_D)*100);

        C_vdw_hy(i) = (q_fwd-q_bkd)*C_conv/0.02;				% Capacitance calculated as charge stored in electrode per unit change in potential

    end
    
    
    devn = (C_vdw_hy-C_expt)./C_expt;
    
    
    sse = sum(devn.^2)
    rmsd = sqrt(sum(devn.^2)/length(devn))
    
    psiD_range = psiD_range*((kb*T)/e);
    p((j-1)*2+1)=plot(psiD_range,C_vdw_hy,[color(j),'-'],'linewidth',2);
    hold on;
    p((j-1)*2+2)=plot(psiD_range,C_expt,[color(j),'x'],'linewidth',2);
    ylabel('Capacitance (\mu F/cm^2)','fontsize',14);
    xlabel('Potential vs. PZC (V)','fontsize',14)
    set(gca,'fontsize',14);
    ylim([0 100]);   
end

legend(p,{['Fit (', num2str(0.1) ,' M)'],['Expt (', num2str(0.1) ,' M)'],['Fit (', num2str(0.04) ,' M)'],['Expt (', num2str(0.04) ,' M)'],['Fit (', num2str(0.005) ,' M)'],['Expt (', num2str(0.005) ,' M)'] });
rmsd


function differ=odefcn(x,y,c0,zplus,zminus,aplus3c0,aminus3c0,aratio,epsilonpw,epsilonmw,sigmapw,sigmamw,Linf,a,rho,wall_on_both_sides)
psi=y(1);
dpsidx=y(2);

ap=((aplus3c0/c0)^(1/3));
am=((aminus3c0/c0)^(1/3));
apm = (ap+am)/2;

ULJplus_wall = @(wallpos) -4*pi*rho*epsilonpw*(apm^2)*((sigmapw^6/(2*(wallpos-x)^4))-((sigmapw^12/(5*(wallpos-x)^10))));

ULJminus_wall = @(wallpos) -4*pi*rho*epsilonmw*(apm^2)*((sigmamw^6/(2*((wallpos)-(x))^4))-((sigmamw^12/(5*((wallpos)-(x))^10)))) ;

spacing = (3.20e-10/(sqrt(3)*apm));

ULJplus = ULJplus_wall(0); % + ULJplus_wall(-spacing) + ULJplus_wall(-2*spacing);
if (wall_on_both_sides==1)
    ULJplus=ULJplus  + ULJplus_wall(Linf); % + ULJplus_wall(Linf+spacing) + ULJplus_wall(Linf+2*spacing);
end

ULJminus= ULJminus_wall(0); % + ULJminus_wall(-spacing) + ULJminus_wall(-2*spacing);
if (wall_on_both_sides==1)
  ULJminus=ULJminus    + ULJminus_wall(Linf); % + ULJminus_wall(Linf+spacing) + ULJminus_wall(Linf+2*spacing);
end

if(x==0)
     ULJplus=inf;
     ULJminus=inf;
end    

if (x==Linf && wall_on_both_sides==1)
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
differ(2)=(nminus-nplus)./(zplus+zminus);	 % Poisson Equation

end

% boundary condition function for function named odefcn
function res=bcfcn(psia,psib,psiD,wall_on_both_sides)

if (wall_on_both_sides==0)
    res=[ psia(1)-psiD
          psib(1)
    ];
else
      res=[ psia(1)-psiD   %(4.28+0.244+psiD)
            psib(1)        %(4.28+0.244)
    ];  
end
end

