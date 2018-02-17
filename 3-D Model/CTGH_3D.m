%This program will perform a 3-D simulation on the CTGH and give the
%overall heat transfer, give a vertical temperature and pressure
%distribution along the sub-bundles.
function CTGH_3D(THEEM_model)
load('THEEM_Input_3D.mat');
%% Bundle Geometry
i=1; %Allows CTGH_Geom to start
[~,~,H,tubes_vol,~,~,tubes,D_in,section,L_tube_avg,vol_cells_gap,slice_total,slice_holder,R_co,vol_wid]=CTGH_geom(THEEM_model,i);
[~,~,mu_l,~,rho_l,~,rho_g,mu_g,~,~,~]=Material_prop(liquid,gas,tube_material,T_l_in,T_g_in,P_g_in);
n=bundles*section; %Total number of vertical sections of CTGH
C_f=1.0; %Coefficient of turning loss
D_manifold=0.280; %Inner diameter of manifold [m]
disk_thick=.003; %Thickness of plate separating bundles
H_sub=D_out*ST*((layer_num+1)/2); %Height of sub-bundle, excluding spacer disk [m]
H_bank=(H_sub+disk_thick)*bundles; %Height of entire tube bank, including spacer disk [m]
L_manifold=H_bank; %Length of manifold [m]
R_c_avg=L_tube_avg/(loops*pi*2);
if L_manifold<=30
    alpha_press=0.5; %Manifold pressure recovery factor at first volume
    gamma_press=0.146; %Manifold pressure recovery factor increment
else
    alpha_press=0.6; %Manifold pressure recovery factor at first volume
    gamma_press=0.15; %Manifold pressure recovery factor increment
end
%% Liquid Mass Flow Rate Distribution
%Values obtained using 0-D calculations
Flow_area_l=tubes*pi/4*D_in^2; %Sum total of flow area inside all tubes
v_l=m_l/(rho_l*Flow_area_l); %Average velocity of salt
Re_l_avg=rho_l*v_l*D_in/mu_l; %Salt Reynolds number
De_l_avg=Re_l_avg*sqrt(D_in/(2*R_c_avg));
%End of 0-D calculations
f_c_l=(64/Re_l_avg)/(1-(1-(11.6/De_l_avg)^0.45)^2.2); %Average pipe friction factor for CTGH
F_c_l=tubes_vol*(pi/4)*D_in^2; %Port (outlet leaving manifold) cross sectional area [m^2]
F_m_l=(pi/4)*D_manifold^2; %Manifold cross sectional area [m^2]
zeta_l=zeros(n+1,1);
k_l=zeros(n+1,1);
w_l=zeros(n+1,1);
Q_l=zeros(n+1,1);
R_l=zeros(n+1,1);
J_l=zeros(n+1,1);
B_l=zeros(n+1,1);
u_c_l=zeros(n+1,1);
m_l_manifold=zeros(n+1,1);
m_l_2_D=zeros(n+1,1);
X=0:L_manifold/n:L_manifold; %Evenly distribute X along manifold starting at x=0
x=X/L_manifold; %Nondimensionalize X in terms of the length of the manifold
%Starting values of vectors:
w_l(1)=1;
u_c_l(1)=0;
m_l_manifold(1)=w_l(1)*m_l/entry;
for i=1:n
zeta_l(i)=(1+C_f+f_c_l*L_tube_avg/D_in);
k_l(i)=alpha_press+2*gamma_press*log(w_l(i));
Q_l(i)=2/(3*zeta_l(i))*k_l(i)*(F_c_l*n/F_m_l)^2;
R_l(i)=-(f_c_l*L_tube_avg/(4*D_in*zeta_l(i)))*(F_c_l*n/F_m_l)^2;
B_l(i)=sign(R_l(i)+sqrt(Q_l(i)^3+R_l(i)^2))*abs((R_l(i)+sqrt(Q_l(i)^3+R_l(i)^2)))^(1/3)+sign(R_l(i)-sqrt(Q_l(i)^3+R_l(i)^2))*abs((R_l(i)-sqrt(Q_l(i)^3+R_l(i)^2)))^(1/3); %Cubic root was giving the complex roots. Using sign(x)*abs(x)^(1/3) forces Matlab to find the real cubic root.
J_l(i)=sign(R_l(i)+sqrt(Q_l(i)^3+R_l(i)^2))*abs((R_l(i)+sqrt(Q_l(i)^3+R_l(i)^2)))^(1/3)-sign(R_l(i)-sqrt(Q_l(i)^3+R_l(i)^2))*abs((R_l(i)-sqrt(Q_l(i)^3+R_l(i)^2)))^(1/3);
w_l(i+1)=exp(-B_l(i)*x(i+1)/2)*(sin(sqrt(3)/2*J_l(i)*(1-x(i+1)))/sin(sqrt(3)*J_l(i)/2));
u_c_l(i+1)=(F_m_l/(2*n*F_c_l))*exp(-B_l(i)*x(i+1)/2)*(B_l(i)*sin(sqrt(3)*J_l(i)*(1-x(i+1))/2)+sqrt(3)*J_l(i)*cos(sqrt(3)*J_l(i)*(1-x(i+1))/2))/(sin(sqrt(3)*J_l(i)/2));
m_l_manifold(i+1)=w_l(i+1)*m_l/entry;
m_l_2_D(i+1)=u_c_l(i+1)*m_l/entry*F_c_l/F_m_l;
end
%% Gas Mass Flow Rate Distribution
D_center_bund=2*R_ci;
F_c_g=H*pi*D_center_bund; %Port (outlet leaving manifold) cross sectional area [m^2]
F_m_g=(pi/4)*D_center_bund^2; %Manifold cross sectional area [m^2]
[~,Eu] = StaggeredPressureDrop(ST,SL,u_max_app,rho_g,N_L,Re_g);
f_c_g=Eu;
zeta_g=zeros(n+1,1);
k_g=zeros(n+1,1);
w_g=zeros(n+1,1);
Q_g=zeros(n+1,1);
R_g=zeros(n+1,1);
J_g=zeros(n+1,1);
B_g=zeros(n+1,1);
u_c_g=zeros(n+1,1);
m_g_center_bund=zeros(n+1,1);
m_g_2_D=zeros(n+1,1);
Y=0:H_bank/n:H_bank; %Evenly distribute Y along bundle starting at y=0
y=Y/H_bank; %Nondimensionalize Y in terms of the height of the bundle
%Starting values of vectors:
w_g(1)=1;
u_c_g(1)=0;
m_g_center_bund(1)=w_g(1)*m_g;
for i=1:n
zeta_g(i)=(1+C_f+f_c_g*H_bank/D_in);
k_g(i)=alpha_press+2*gamma_press*log(w_g(i));
Q_g(i)=2/(3*zeta_g(i))*k_g(i)*(F_c_g*n/F_m_g)^2;
R_g(i)=-(f_c_g*L_tube_avg/(4*D_in*zeta_g(i)))*(F_c_g*n/F_m_g)^2;
B_g(i)=sign(R_g(i)+sqrt(Q_g(i)^3+R_g(i)^2))*abs((R_g(i)+sqrt(Q_g(i)^3+R_g(i)^2)))^(1/3)+sign(R_g(i)-sqrt(Q_g(i)^3+R_g(i)^2))*abs((R_g(i)-sqrt(Q_g(i)^3+R_g(i)^2)))^(1/3); %Cubic root was giving the complex roots. Using sign(x)*abs(x)^(1/3) forces Matlab to find the real cubic root.
J_g(i)=sign(R_g(i)+sqrt(Q_g(i)^3+R_g(i)^2))*abs((R_g(i)+sqrt(Q_g(i)^3+R_g(i)^2)))^(1/3)-sign(R_g(i)-sqrt(Q_g(i)^3+R_g(i)^2))*abs((R_g(i)-sqrt(Q_g(i)^3+R_g(i)^2)))^(1/3);
w_g(i+1)=exp(-B_g(i)*y(i+1)/2)*(sin(sqrt(3)/2*J_g(i)*(1-y(i+1)))/sin(sqrt(3)*J_g(i)/2));
u_c_g(i+1)=(F_m_g/(2*n*F_c_g))*exp(-B_g(i)*y(i+1)/2)*(B_g(i)*sin(sqrt(3)*J_g(i)*(1-y(i+1))/2)+sqrt(3)*J_g(i)*cos(sqrt(3)*J_g(i)*(1-y(i+1))/2))/(sin(sqrt(3)*J_g(i)/2));
m_g_center_bund(i+1)=w_g(i+1)*m_g;
m_g_2_D(i+1)=u_c_g(i+1)*m_g*F_c_g/F_m_g;
end
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