%This program will perform a 3-D simulation on the CTGH and give the
%overall heat transfer, give a vertical temperature and pressure
%distribution along the sub-bundles.
clc;clear;
load('THEEM_Input_3D.mat');
%% Liquid Mass Flow Rate Distribution
i=1; %Allows CTGH_Geom to start
[L,R_curv,H,tubes_vol,N_T,N_L,tubes,D_in,section,L_tube_avg,vol_cells_gap,slice_total,slice_holder,R_co,vol_wid]=CTGH_geom(THEEM_model,i);
n=bundles*section; %Total number of vertical sections of CTGH
alpha_press=0.5; %Manifold pressure recovery factor
gamma_press=-0.25; %Manifold pressure recovery factor increment
C_f=1.0; %Coefficient of turning loss
D_manifold=0.280; %Inner diameter of manifold [m]
L_manifold=5.9087; %Length of manifold [m]
%Values obtained using 0-D calculations
Flow_area=tubes*pi/4*D_in^2; %Sum total of flow area inside all tubes
v_l=m_l/(rho_l*Flow_area); %Average velocity of salt
Re_l_avg=rho_l*v_l*D_in/mu_l; %Salt Reynolds number
R_c_avg=L_tube_avg/(loops*pi*2);
De_l_avg=Re_l_avg*sqrt(D_in/(2*R_c_avg));
%End of 0-D calculations
f_c=(64/Re_l_avg)/(1-(1-(11.6/De_l_avg)^0.45)^2.2); %Average pipe friction factor for CTGH
F_c=tubes_vol*(pi/4)*D_in^2; %Port (outlet leaving manifold) cross sectional area [m^2]
F_m=(pi/4)*D_manifold^2; %Manifold cross sectional area [m^2]
zeta=zeros(n+1,1);
w=zeros(n+1,1);
Q=zeros(n+1,1);
R=zeros(n+1,1);
J=zeros(n+1,1);
B=zeros(n+1,1);
u_c=zeros(n+1,1);
m_l_manifold=zeros(n+1,1);
m_l_2_D=zeros(n+1,1);
X=0:L_manifold/n:L_manifold; %Evenly distribute X along manifold starting at x=0
x=X/L_manifold; %Nondimensionalize X in terms of the length of the manifold
%Starting values of vectors:
w(1)=1;
u_c(1)=0;
m_l_vol(1)=0;
m_l_manifold(1)=w(1)*m_l/entry;
for i=1:n
zeta(i)=(1+C_f+f_c*L_tube_avg/D_in);
Q(i)=2/(3*zeta(i))*(alpha_press+2*gamma_press*log(w(i)))*(F_c*n/F_m)^2;
R(i)=-(f_c*L_tube_avg/(4*D_in*zeta(i)))*(F_c*n/F_m)^2;
B(i)=sign(R(i))*abs((R(i)+sqrt(Q(i)^3+R(i)^2)))^(1/3)+sign(R(i))*abs((R(i)-sqrt(Q(i)^3+R(i)^2)))^(1/3); %Cubic root was giving the complex roots. Using sign(x)*abs(x)^(1/3) forces Matlab to find the real cubic root.
J(i)=sign(R(i))*abs((R(i)+sqrt(Q(i)^3+R(i)^2)))^(1/3)-sign(R(i))*abs((R(i)-sqrt(Q(i)^3+R(i)^2)))^(1/3);
w(i+1)=exp(-B(i)*x(i+1)/2)*(sin(sqrt(3)*J(i)*(1-x(i+1))/2)/sin(sqrt(3)*J(i)/2));
u_c(i+1)=(F_m/(2*n*F_c))*exp(-B(i)*x(i+1)/2)*(B(i)*sin(sqrt(3)*J(i)*(1-x(i+1))/2)+sqrt(3)*J(i)*cos(sqrt(3)*J(i)*(1-x(i+1))/2))/(sin(sqrt(3)*J(i)/2));
m_l_manifold(i+1)=w(i+1)*m_l/entry;
m_l_2_D(i+1)=u_c(i+1)*m_l/entry*F_c/F_m;
end
%% Gas Mass Flow Rate Distribution
m_g_2_D=repmat(m_g/n,n+1,1);
%% Run 2-D simulation and store matrix outputs
T_l_out_store=cell(n+1,1);
T_g_out_store=cell(n+1,1);
P_l_out_store=cell(n+1,1);
P_g_out_store=cell(n+1,1);
T_l_store=cell(n+1,1);
T_g_store=cell(n+1,1);
P_l_store=cell(n+1,1);
P_g_store=cell(n+1,1);
Q_store=cell(n+1,1);
epsilon_store=zeros(n+1,1);
U_store=zeros(n+1,1);
Area_store=zeros(n+1,1);
T_l_out_store{1}=0;
T_g_out_store{1}=0;
P_l_out_store{1}=0;
T_l_store{1}=0;
T_g_store{1}=0;
P_l_store{1}=0;
P_g_store{1}=0;
Q_store{1}=0;
for i=2:n+1
    [T_l_out,T_g_out,P_l_out,T_l,T_g,P_l,P_g,Q,epsilon,U_avg,A_total]=CTGH_2D(THEEM_model,m_l_2_D(i),m_g_2_D(i));
    T_l_out_store{i}=T_l_out;
    T_g_out_store{i}=T_g_out;
    P_l_out_store{i}=P_l_out;
    P_g_out_store{i}=P_g(size(P_g,1),:);
    T_l_store{i}=T_l;
    T_g_store{i}=T_g;
    P_l_store{i}=P_l;
    P_g_store{i}=P_g;
    Q_store{i}=Q;
    U_store(i)=U_avg;
    Area_store(i)=A_total;
    epsilon_store(i)=epsilon;
end
save('3-D Model/THEEM_Output_3D.mat');

