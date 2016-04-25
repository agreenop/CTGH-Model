%This program will perform a 3-D simulation on the CTGH and give the
%overall heat transfer, give a vertical temperature and pressure
%distribution along the sub-bundles.
clc;clear;
load('THEEM_3D_Input.mat');
%% Liquid Mass Flow Rate Distribution
i=14; %Allows CTGH_Geom to start
[tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t,R_curv,loops,spacers,section,bundles,L_tube_avg]=CTGH_geom(tube_material,D_out,t,ST,SL,entry,i);
n=bundles*section; %Total number of vertical sections of CTGH
alpha_press=0.5; %Manifold pressure recovery factor
gamma_press=-0.25; %Manifold pressure recovery factor increment
C_f=1.0; %Coefficient of turning loss
D_manifold=0.280; %Inner diameter of manifold [m]
L_manifold=5.9087; %Length of manifold [m]
%Values obtained from 0-D model results
Re_l_avg=1660.1;
R_c_avg=L_tube_avg/(loops*pi*2);
De_l_avg=Re_l_avg*sqrt(D_in/(2*R_c_avg));
%End of 0-D results
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
w(1)=1;
u_c(1)=0;
m_l_vol(1)=0;
m_l_manifold(1)=w(1)*m_l;
for i=1:n
zeta(i)=(1+C_f+f_c*L_tube_avg/D_in);
Q(i)=2/(3*zeta(i))*(alpha_press+2*gamma_press*log(w(i)))*(F_c*n/F_m)^2;
R(i)=-(f_c*L_tube_avg/(4*D_in*zeta(i)))*(F_c*n/F_m)^2;
B(i)=sign(R(i))*abs((R(i)+sqrt(Q(i)^3+R(i)^2)))^(1/3)+sign(R(i))*abs((R(i)-sqrt(Q(i)^3+R(i)^2)))^(1/3); %Cubic root was giving the complex roots. Using sign(x)*abs(x)^(1/3) forces Matlab to find the real cubic root.
J(i)=sign(R(i))*abs((R(i)+sqrt(Q(i)^3+R(i)^2)))^(1/3)-sign(R(i))*abs((R(i)-sqrt(Q(i)^3+R(i)^2)))^(1/3);
w(i+1)=exp(-B(i)*x(i+1)/2)*(sin(sqrt(3)*J(i)*(1-x(i+1))/2)/sin(sqrt(3)*J(i)/2));
u_c(i+1)=(F_m/(2*n*F_c))*exp(-B(i)*x(i+1)/2)*(B(i)*sin(sqrt(3)*J(i)*(1-x(i+1))/2)+sqrt(3)*J(i)*cos(sqrt(3)*J(i)*(1-x(i+1))/2))/(sin(sqrt(3)*J(i)/2));
m_l_manifold(i+1)=w(i+1)*m_l;
m_l_2_D(i+1)=u_c(i+1)*m_l*F_c/F_m;
end
%% Gas Mass Flow Rate Distribution
m_g_2_D=repmat(m_g/n,n+1,1);
%% Run 2-D simulation and store matrix outputs
T_l_out_store=cell(n+1,1);
T_g_out_store=cell(n+1,1);
P_l_out_store=cell(n+1,1);
T_l_store=cell(n+1,1);
T_g_store=cell(n+1,1);
P_l_store=cell(n+1,1);
P_g_store=cell(n+1,1);
Q_store=cell(n+1,1);
epsilon_store=zeros(n+1,1);
T_l_out_store{1}=0;
T_g_out_store{1}=0;
P_l_out_store{1}=0;
T_l_store{1}=0;
T_g_store{1}=0;
P_l_store{1}=0;
P_g_store{1}=0;
Q_store{1}=0;
for i=2:n+1
    [T_l_out,T_g_out,P_l_out,T_l,T_g,P_l,P_g,Q,epsilon]=CTGH_2D_calculations(m_l_2_D(i),m_g_2_D(i));
    T_l_out_store{i}=T_l_out;
    T_g_out_store{i}=T_g_out;
    P_l_out_store{i}=P_l_out;
    T_l_store{i}=T_l;
    T_g_store{i}=T_g;
    P_l_store{i}=P_l;
    P_g_store{i}=P_g;
    Q_store{1}=Q;
    epsilon_store(i)=epsilon;
end
save('THEEM_3D_Output.mat');

%% Plot liquid manifold flow distribution
% figure(1)
% plot(X(2:n+1),m_l_vol(2:n+1)) %Individual port flow rates distribution over manifolds
% hold on
% plot([0,X(n+1)],[m_l/n,m_l/n]) %Average port flow rate
% hold off
% figure(2)
% plot(X,m_l_manifold) %Manifold flow rate
%% Plot liquid manifold flow distribution in 3-D
% y=R_curv:D_manifold/10:R_curv+D_manifold;
% y3=-D_manifold/2:D_manifold/10:D_manifold/2;
% X1=repmat(X,11,1);
% Y1=repmat(y',1,n+1);
% Y2=repmat(-1*y',1,n+1);
% Y3=repmat(y3',1,n+1);
% Z=repmat((m_l_manifold./4)',11,1);
% figure(3)
% hold on
% surf(Y1,X1,Z,'linestyle','none');
% surf(Y2,X1,Z,'linestyle','none');
% surf(Y3,X1,Z,'linestyle','none');
% hold off
% view(2)
% axis([-2,2,0,10])