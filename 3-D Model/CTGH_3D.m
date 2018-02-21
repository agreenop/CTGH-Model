%This program will perform a 3-D simulation on the CTGH and give the
%overall heat transfer, give a vertical temperature and pressure
%distribution along the sub-bundles.
function CTGH_3D
load('THEEM_Input_3D.mat');
%% Bundle Geometry
i=1; %Allows CTGH_Geom to start
[~,~,H,tubes_vol,N_T,~,~,D_in,section,L_tube_avg]=CTGH_geom(THEEM_model,i);
CTGH_0D(THEEM_model,i)
load('3-D Model/THEEM_Output_temp_0D.mat');
n=bundles*section; %Total number of vertical sections of CTGH
C_f=1.0; %Coefficient of turning loss
D_manifold=0.280; %Inner diameter of manifold [m]
t_manifold=0.003; %Manifold thickness [m]
L_manifold=H_bank; %Length of manifold [m]
if L_manifold/D_manifold<=30
    alpha_press=0.5; %Manifold pressure recovery factor at first volume
    gamma_press=0.146; %Manifold pressure recovery factor increment
else
    alpha_press=0.6; %Manifold pressure recovery factor at first volume
    gamma_press=0.15; %Manifold pressure recovery factor increment
end
%% Liquid Mass Flow Rate Distribution
f_c_l=f_l; %Average friction factor for tubes in tube bank calculated by 0-D code
F_c_l=tubes_vol*(pi/4)*D_in^2; %Port (outlet leaving manifold) cross sectional area [m^2]
F_m_l=(pi/4)*D_manifold^2; %Manifold cross sectional area [m^2]
W_0_l=m_l/(entry*rho_l*F_m_l); %Manifold inlet velocity
Re_l_man=zeros(n+1,1);
f_m_l=zeros(n+1,1);
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
delta_p=zeros(n+1,1);
P_l_inlet=zeros(n+1,1);
X=0:L_manifold/n:L_manifold; %Evenly distribute X along manifold starting at x=0
x=X/L_manifold; %Nondimensionalize X in terms of the length of the manifold
%Starting values of vectors:
w_l(1)=1;
u_c_l(1)=0;
m_l_manifold(1)=w_l(1)*m_l/entry;
delta_p(1)=0;
P_l_inlet(1)=P_l_in;
for i=1:n
Re_l_man(i)=rho_l*W_0_l*w_l(i)*D_manifold/mu_l; %Reynolds number of liquid in manifold
if Re_l_man(i)<=2200
    f_m_l(i)=64/Re_l_man(i); %Friction factor in manifold if laminar
elseif Re_l_man(i)>2200 && Re_l_man(i)<=10^5
    f_m_l(i)=0.790*log(Re_l_man(i)-1.64)^(-2); %Friction factor if transitional
else
    f_m_l(i)=0.0032+0.221/(Re_l_man(i)^0.237); %Friction factor if turbulent
end
zeta_l(i)=(1+C_f+f_c_l*L_tube_avg/D_in);
k_l(i)=alpha_press+2*gamma_press*log(w_l(i));
Q_l(i)=2/(3*zeta_l(i))*k_l(i)*(F_c_l*n/F_m_l)^2;
R_l(i)=-(f_m_l(i)*L_manifold/(4*D_manifold*zeta_l(i)))*(F_c_l*n/F_m_l)^2;
B_l(i)=sign(R_l(i)+sqrt(Q_l(i)^3+R_l(i)^2))*abs((R_l(i)+sqrt(Q_l(i)^3+R_l(i)^2)))^(1/3)+sign(R_l(i)-sqrt(Q_l(i)^3+R_l(i)^2))*abs((R_l(i)-sqrt(Q_l(i)^3+R_l(i)^2)))^(1/3); %Cubic root was giving the complex roots. Using sign(x)*abs(x)^(1/3) forces Matlab to find the real cubic root.
J_l(i)=sign(R_l(i)+sqrt(Q_l(i)^3+R_l(i)^2))*abs((R_l(i)+sqrt(Q_l(i)^3+R_l(i)^2)))^(1/3)-sign(R_l(i)-sqrt(Q_l(i)^3+R_l(i)^2))*abs((R_l(i)-sqrt(Q_l(i)^3+R_l(i)^2)))^(1/3);
w_l(i+1)=exp(-B_l(i)*x(i+1)/2)*(sin(sqrt(3)/2*J_l(i)*(1-x(i+1)))/sin(sqrt(3)*J_l(i)/2));
u_c_l(i+1)=(F_m_l/(2*n*F_c_l))*exp(-B_l(i)*x(i+1)/2)*(B_l(i)*sin(sqrt(3)*J_l(i)*(1-x(i+1))/2)+sqrt(3)*J_l(i)*cos(sqrt(3)*J_l(i)*(1-x(i+1))/2))/(sin(sqrt(3)*J_l(i)/2));
m_l_manifold(i+1)=w_l(i+1)*m_l/entry;
m_l_2_D(i+1)=u_c_l(i+1)*m_l*F_c_l/F_m_l;
delta_p(i+1)=L_manifold*f_m_l(i)/(4*D_manifold*B_l(i)*sin(sqrt(3)/2*J_l(i))^2)*(exp(-B_l(i)*x(i+1))-1)-L_manifold*f_m_l(i)/(4*4*D_manifold*(B_l(i)^2+3*J_l(i)^2)*sin(sqrt(3)/2*J_l(i))^2)*(B_l(i)*exp(-B_l(i)*x(i+1))*cos(sqrt(3)*J_l(i)*(1-x(i+1)))-sqrt(3)*J_l(i)*exp(-B_l(i)*x(i+1))*sin(sqrt(3)*J_l(i)*(1-x(i+1)))+B_l(i)*cos(sqrt(3)*J_l(i))+sqrt(3)*J_l(i)*sin(sqrt(3)*J_l(i)))-k_l(i)*exp(-B_l(i)*x(i+1))*sin(sqrt(3)*J_l(i)*(1-x(i+1))/2)^2/sin(sqrt(3)*J_l(i)/2)^2;
P_l_inlet(i+1)=P_l_in+delta_p(i+1)*rho_l*W_0_l^2*10^-5;
end
%% Gas Mass Flow Rate Distribution
D_center_bund=2*R_ci;
if H_bank/D_center_bund<=30
    alpha_press_g=0.5; %Manifold pressure recovery factor at first volume for gas
    gamma_press_g=0.146; %Manifold pressure recovery factor increment for gas
else
    alpha_press_g=0.6; %Manifold pressure recovery factor at first volume for gas
    gamma_press_g=0.15; %Manifold pressure recovery factor increment for gas
end
if mod(N_T,2)==1 %If number of tubes per volume is even
    F_c_g=(H-N_T*D_out/2)*(pi*D_center_bund-entry*(D_manifold+2*t_manifold)); %Port (outlet leaving manifold) cross sectional area [m^2]
else % If odd
    F_c_g=(H-(N_T+1)*D_out/2)*pi*D_center_bund; %Port (outlet leaving manifold) cross sectional area [m^2]
end
F_m_g=(pi/4)*D_center_bund^2; %Manifold cross sectional area [m^2]
W_0_g=m_g/(rho_g*F_m_g); %Gas manifold inlet velocity
Re_g_man=zeros(n+1,1);
u_c_g=zeros(n+1,1);
m_g_center_bund=zeros(n+1,1);
m_g_2_D=zeros(n+1,1);
Eu_man_g=zeros(n+1,1);
P_g_inlet=zeros(n+1,1);
w_g=zeros(n+1,1);
Y=0:H_bank/n:H_bank; %Evenly distribute Y along bundle starting at y=0
y=Y/H_bank; %Nondimensionalize Y in terms of the height of the bundle
for i=1:n+1
    w_g(i)=1-y(i);
    Re_g_man(i)=rho_g*W_0_g*w_g(i)*D_manifold/mu_g; %Reynolds number of liquid in manifold
    u_c_g(i)=F_m_g/(n*F_c_g);
    m_g_center_bund(i)=w_g(i)*m_g;
    m_g_2_D(i)=u_c_g(i)*m_g*F_c_g/F_m_g;
    if Re_g_man(i)<=2200
        Eu_man_g(i)=(alpha_press_g-16*H_bank/(D_center_bund*Re_g_man(1)))*(1-(1-y(i))^2)-2*gamma_press_g*((1-y(i))^2*log(1-y(i))+1/2*y(i)*(2-y(i)));
    elseif Re_g_man(i)>2200 && Re_g_man(i)<=10^5
        Eu_man_g(i)=alpha_press_g*(1-(1-y(i))^2)-0.058*H_bank/(D_center_bund*Re_g_man(1)^0.25)*(1-(1-y(i))^2.75)-2*gamma_press_g*((1-y(i))^2*log(1-y(i))+1/2*y(i)*(2-y(i)));
    else
        Eu_man_g(i)=alpha_press_g*(1-(1-y(i))^2)-0.0032*H_bank/(D_center_bund*6)*(1-(1-y(i))^3)-0.04*H_bank/(D_center_bund*Re_g_man(1)^0.237)*(1-(1-y(i))^2.763)-2*gamma_press_g*((1-y(i))^2*log(1-y(i))+1/2*y(i)*(2-y(i)));
    end
    P_g_inlet(i)=rho_g*W_0_g^2*Eu_man_g(i)*10^-5+P_g_in; %Pressure in bar
end
delete('3-D Model/THEEM_Output_temp_0D.mat');
%% Run 2-D simulation and store matrix outputs
T_l_out_store=cell(n,1);
T_g_out_store=cell(n,1);
P_l_out_store=cell(n,1);
P_g_out_store=cell(n,1);
T_l_store=cell(n,1);
T_g_store=cell(n,1);
P_l_store=cell(n,1);
P_g_store=cell(n,1);
Q_store=cell(n,1);
epsilon_store=zeros(n,1);
U_store=zeros(n,1);
Area_store=zeros(n,1);
for i=1:n
    h=waitbar((i-1)/n,sprintf('Progress = %2.2f%%',(i-1)/n*100));
    [T_l_out,T_g_out,P_l_out,T_l,T_g,P_l,P_g,Q,epsilon,U_avg,A_total]=CTGH_2D(THEEM_model,m_l_2_D(i+1),m_g_2_D(n+1-i),P_l_inlet(i+1),P_g_inlet(n+1-i));
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
    delete(h)
end

save('3-D Model/THEEM_Output_3D.mat');