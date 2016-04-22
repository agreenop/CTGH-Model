%This program will perform a 3-D simulation on the CTGH and give the
%overall heat transfer, give a vertical temperature and pressure
%distribution along the sub-bundles.
clc;clear;
load('THEEM_3D_Input.mat');
i=1; %Allows CTGH_Geom to start
[tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t,R_curv,loops,spacers,section,bundles,L_tube_avg]=CTGH_geom(tube_material,D_out,t,ST,SL,entry,i);
n=bundles*section; %Total number of vertical sections of CTGH
alpha_press=0.5; %Manifold pressure recovery factor
gamma_press=-0.25; %Manifold pressure recovery factor increment
C_f=1.0; %Coefficient of turning loss
D_manifold=0.280; %Inner diameter of manifold [m]
L_manifold=5.9087; %Length of manifold [m]
%% Values obtained from 0-D model results
Re_l_avg=1660.1;
R_c_avg=L_tube_avg/(loops*pi*2);
De_l=Re_l*sqrt(D_in/(2*R_c_avg));
%End of 0-D results
%%
f_c=(64/Re_l_avg)/(1-(1-(11.6/De_l)^0.45)^2.2); %Average pipe friction factor for CTGH
F_c=tubes_vol*(pi/4)*D_in^2; %Port (outlet leaving manifold) cross sectional area [m^2]
F_m=(pi/4)*D_manifold^2; %Manifold cross sectional area [m^2]
zeta=zeros(n+1,1);
w=zeros(n+1,1);
Q=zeros(n+1,1);
R=zeros(n+1,1);
J=zeros(n+1,1);
B=zeros(n+1,1);
u_c=zeros(n+1,1);
w(1)=1;
for i=1:n
zeta(i)=(1+C_f+f_c*L_tube_avg/D_in);
Q(i)=2/(3*zeta(i))*(alpha_press+2*gamma_press*log(w(i)))*(F_c*n/F_m)^2;
R(i)=-(f_c*L_tube_avg/(4*D_in*zeta(i)))*(F_c*n/F_m)^2;
B(i)=(R(i)+sqrt(Q(i)^3+R(i)^2))^(1/3)+(R(i)-sqrt(Q(i)^3+R(i)^2))^(1/3);
J(i)=(R(i)+sqrt(Q(i)^3+R(i)^2))^(1/3)-(R(i)-sqrt(Q(i)^3+R(i)^2))^(1/3);
w(i+1)
u_c(i+1)
end
save('THEEM_3D_Output.mat');