clear
clc
close all

xmesh=[4.8:0.0001:10.0];		% array of inter-molecular (sulphate-sulphate) distances

d_SO=1.49;				% Equilibrium S-O distance in sulphate         		
d_OO=2.43316;				% Equilibrium O-O distance in sulphate 			
a=d_OO/sqrt(2);
sigma_SS=3.55000;			% LJ parameter sigma for S-S interactions
sigma_OO=3.65000;			% LJ parameter sigma for O-O interactions
sigma_SO=(sigma_SS+sigma_OO)/2;	% LJ parameter sigma for S-O interactions (Lorentz-Berthelot Combining Rules)

eps_SS=1.046700;			% LJ parameter epsilon for S-S interactions
eps_OO=0.837400;			% LJ parameter epsilon for O-O interactions
eps_SO=sqrt(eps_SS*eps_OO);		% LJ parameter epsilon for S-O interactions (Lorentz-Berthelot Combining Rules)

y1_hand=zeros(1,length(xmesh));
y2=zeros(1,length(xmesh));

% computing sulphate-sulphate interactions as the sum of atom-atom interactions
for i=1:length(xmesh)
    d_SS=xmesh(i);
    d_SO1=d_SS-(a/2);
    d_SO2=d_SS+(a/2);
    dist_OO1=d_SS-a;
    dist_OO2=d_SS+a;
    dist_OO3=d_SS;
    U_SS=-4*eps_SS*((((sigma_SS)^6)/((d_SS^6))) - (((sigma_SS)^12)/((d_SS^12))));
    U_SO1=-4*eps_SO*((((sigma_SO)^6)/((d_SO1^6))) - (((sigma_SO)^12)/((d_SO1^12))));
    U_SO2=-4*eps_SO*((((sigma_SO)^6)/((d_SO2^6))) - (((sigma_SO)^12)/((d_SO2^12))));
    U_OO1=-4*eps_OO*((((sigma_OO)^6)/((dist_OO1^6))) - (((sigma_OO)^12)/((dist_OO1^12))));
    U_OO2=-4*eps_OO*((((sigma_OO)^6)/((dist_OO2^6))) - (((sigma_OO)^12)/((dist_OO2^12))));
    U_OO3=-4*eps_OO*((((sigma_OO)^6)/((dist_OO3^6))) - (((sigma_OO)^12)/((dist_OO3^12))));
    U_total= U_SS + (2*U_SO1) + (2*U_SO2)+U_OO1+U_OO2+(2*U_OO3); 
    y1_hand(i)=U_total;
end

% sigma-epsilon fitting routine
k=1;
ep=1;
for i=2:length(y1_hand)
    if(y1_hand(i)*y1_hand(i-1) ==0 || y1_hand(i)*y1_hand(i-1)<0)
        k=i;
    end

    if((y1_hand(i) < y1_hand(i-1)))
        ep=i;
    end
end

sig=xmesh(k);			% sigma computed as the distance where the potential cuts the x-axis
eps=y1_hand(ep)-y1_hand(k);	% epsilon computed as the well depth

for i=1:length(xmesh)
    d_SS=xmesh(i);
    y2(i)=4*eps*((((sig)^6)/((d_SS^6))) - (((sig)^12)/((d_SS^12))));
end

figure;
set(0,'defaulttextinterpreter','latex')
plot(xmesh,y1_hand,'red','DisplayName','Datapoints-hand');
hold on
plot(xmesh,y2,"Color",'red',"LineStyle",'--','DisplayName','Fit-hand');
legend boxoff
xlabel("Distance from wall (in \AA)");
ylabel("Interaction energy (in kJ/mol)");

box on
grid on
