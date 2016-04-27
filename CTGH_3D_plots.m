clc;clear;
load('THEEM_3D_Output.mat');
%% Plot liquid manifold flow distribution
figure(1)
plot(X(2:n+1),m_l_2_D(2:n+1)) %Individual port flow rates distribution over manifolds
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
for x=2:size(T_l_out_store,1)
    T_l_out_avg(x-1)=mean(T_l_out_store{x});
    T_g_out_avg(x-1)=mean(T_g_out_store{x});
end
figure(3)
plot(X(2:size(X,2)),T_l_out_avg) %Manifold flow rate
title([liquid ' Outlet Temperature Distribution'])
xlabel('Distance from Manifold Inlet (m)')
ylabel([liquid 'Temperature (\circC)'])
figure(4)
plot(X(2:size(X,2)),T_g_out_avg) %Manifold flow rate
title([gas ' Outlet Temperature Distribution'])
xlabel('Distance from Manifold Inlet (m)')
ylabel([liquid 'Temperature (\circC)'])
Area_tot=3869.5; %m^2
U_avg=
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