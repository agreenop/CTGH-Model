clc;clear;
load('THEEM_3D_Output.mat');
%% Plot liquid manifold flow distribution
figure(1)
plot(X(2:n+1),entry*m_l_2_D(2:n+1)) %Individual port flow rates distribution over manifolds
hold on
plot([0,X(n+1)],[m_l/n,m_l/n]) %Average port flow rate
hold off
title('Port Mass Flow Rate Distribution Along Manifold')
xlabel('Distance from Manifold Inlet (m)')
ylabel('Port Mass Flow Rates (kg/s)')
legend('Calculated Mass Flow Rates','Average Mass Flow Rate')
figure(2)
plot(X,m_l_manifold) %Manifold flow rate
title('Manifold Mass Flow Rate Distribution')
xlabel('Distance from Manifold Inlet (m)')
ylabel('Manifold Mass Flow Rates (kg/s)')
%% Plot Temperature Distribution Outlets
T_l_out_avg=zeros(size(T_l_out_store,1)-1,1);
T_g_out_avg=zeros(size(T_l_out_store,1)-1,1);
P_l_out_avg=zeros(size(T_l_out_store,1)-1,1);
P_g_out_avg=zeros(size(T_l_out_store,1)-1,1);
Q_sum=zeros(size(T_l_out_store,1)-1,1);

for x=2:size(T_l_out_store,1)
    T_l_out_avg(x-1)=mean(T_l_out_store{x});
    T_g_out_avg(x-1)=mean(T_g_out_store{x});
    P_l_out_avg(x-1)=mean(P_l_out_store{x});
    P_g_out_avg(x-1)=mean(P_g_out_store{x});
    Q_sum(x-1)=sum(nansum(Q_store{x-1}));
end
figure(3)
hold on
plot(X(2:size(X,2)),T_l_out_avg) %Manifold flow rate
plot([X(2),X(size(X,2))],[T_l_in,T_l_in])
legend('Outlet','Inlet','location','best')
hold off
title([liquid ' Outlet Temperature Distribution'])
xlabel('Distance from Manifold Inlet (m)')
ylim([500,720])
ylabel([liquid 'Temperature (\circC)'])
figure(4)
hold on
plot(X(2:size(X,2)),T_g_out_avg) %Manifold flow rate
plot([X(2),X(size(X,2))],[T_g_in,T_g_in])
hold off
legend('Outlet','Inlet','location','best')
title([gas ' Outlet Temperature Distribution'])
xlabel('Distance from Liquid Manifold Inlet (m)')
ylabel([gas 'Temperature (\circC)'])
Area_tot=sum(Area_store);
U_mean=mean(U_store);
T_l_out_mean=mean(T_l_out_avg);
T_g_out_mean=mean(T_g_out_avg);
P_l_out_mean=mean(P_l_out_avg);
P_g_out_mean=mean(P_g_out_avg);
T_l_avg=(T_l_out_mean+T_l_in)/2;
T_g_avg=(T_g_out_mean+T_g_in)/2;
P_g_avg=(P_g_out_mean+P_g_in)/2;
[~,Cp_g,~,~,~] = Air_prop(T_g_avg,P_g_avg);
[~,Cp_l,~,~,~,~] = Flibe_prop(T_l_avg);
C_min=min(m_g*Cp_g,m_l*Cp_l);
Q_max=C_min*(T_l_in-T_g_in);
Q_total=sum(Q_sum);
e1=Q_total/Q_max;
LMTD_tot=((T_l_in-T_g_out_mean)-(T_l_out_mean-T_g_in))/log((T_l_in-T_g_out_mean)/(T_l_out_mean-T_g_in));
F_factor=Q_total/(U_mean*Area_tot*LMTD_tot);
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