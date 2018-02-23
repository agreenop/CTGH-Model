clc;clear;close all;
var_empty=[]; %Initializes while loop
while isempty(var_empty) %Will not proceed until user selects a real output file or quits the script
    [FileName,PathName] = uigetfile('*.mat'); %User selects an output file
    if FileName~=0 %A .mat file needs to be selected
        vars = whos('-file',[PathName,FileName]); %Stores variables from .mat file without loading into workspace
        mat_size=numel(vars); %Number of variables in mat file
        load([PathName,FileName],'THEEM_model'); %Loads which THEEM Model the file is using
        if mat_size>50 && strcmp(THEEM_model,'3D')%Check to see if mat file has the correct number of variablesand is a 3D output file. Input file currently has 25 variables, so 50 was chosen to be safe
            load([PathName,FileName])
            var_empty=1;
        else %If incorrect output file, will give error and restart program
            h=errordlg({'This is not a 3D output file.','Please select another file.'},'3D Output File Error');
            uiwait(h);
        end
    else
        clc;clear;
        return; %If user does not select a file, script terminates
    end
end
X_liq=H_bank-Y; %Positions for the gas manifolds
%% Plot port flow distribution
figure(1)
hold on
plot(X_liq(2:n+1),m_l_2_D(2:n+1)) %Individual port flow rates distribution over manifolds
plot([0,X_liq(1)],[m_l/n,m_l/n]) %Average port flow rate
hold off
axis([Y(1),Y(numel(Y)),0,m_l_2_D(n+1)]);
axis 'auto y'
title([liquid ' Flow Rate Distribution to Bundle'])
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Mass Flow Rate (kg/s)')
legend('3D Model','2D Model','Location','northeast')
figure(2)
plot(Y(2:n+1),m_g_2_D(2:n+1)) %Individual port flow rates distribution over manifolds
axis([Y(1),Y(numel(Y)),0,m_g_2_D(n+1)]);
axis 'auto y'
title([gas,' Flow Rate Distribution to Bundle'])
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Mass Flow Rate (kg/s)')
%% Plot port flow distribution
figure(3)
hold on
grid on
manifold_plot=plot(X_liq,entry*m_l_manifold,Y,m_g_center_bund);
hold off
manifold_plot(1).LineStyle='-';
manifold_plot(1).Color='k';
manifold_plot(2).LineStyle='-.';
manifold_plot(2).Color='b';
axis([Y(1),Y(numel(Y)),0,entry*m_l_manifold(1)]);
axis 'auto y'
title('Manifold Flow Rate Vertical Distribution')
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Mass Flow Rate (kg/s)')
legend([liquid ' Mass Flow'],[gas ' Mass Flow'],'location','north')
%% Plot Pressure Distributions in Manifold
figure(4)
plot(X_liq,P_l_inlet) %Liquid Pressure in Manifold
title([liquid ' Pressure in Manifold'])
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Pressure (bar)')
axis([Y(1),Y(numel(Y)),0,P_l_inlet(n+1)]);
axis 'auto y'
figure(5)
plot(Y,P_g_inlet) %Gas Pressure in Gas Manifold
title([gas ' Pressure in Center of Bundle'])
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Pressure (bar)')
axis([Y(1),Y(numel(Y)),0,P_g_inlet(1)]);
axis 'auto y'
%% Plot Liquid Outlet Temperature & Heat Transfer Distribution
T_l_out_avg=zeros(size(T_l_out_store,1),1);
P_l_out_avg=zeros(size(T_l_out_store,1),1);
Q_sum_store=zeros(size(T_l_out_store,1),1);
for i1=1:size(T_l_out_store,1)
    T_l_out_avg(i1)=mean(T_l_out_store{i1});
    P_l_out_avg(i1)=mean(P_l_out_store{i1});
    Q_sum_store(i1)=sum(nansum(Q_store{i1}));
end
figure(6)
plot(X_liq(1:n),T_l_out_avg) %Plot liquid outlet temperature
title([liquid ' Outlet Temperature Distribution'])
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Temperature (\circC)')
axis([Y(1),Y(numel(Y)),0,T_l_out_avg(1)]);
axis 'auto y'
delta_P_l_avg=P_l_inlet(2:n+1)-P_l_out_avg; %Calculate liquid pressure drop across each cross section
figure(7)
plot(X_liq(1:n),delta_P_l_avg) %Plot liquid pressure drop
title([liquid ' Pressure Drop Distribution'])
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Pressure Drop (bar)')
axis([Y(1),Y(numel(Y)),0,delta_P_l_avg(1)]);
axis 'auto y'
figure(8)
plot(X_liq(1:n),Q_sum_store) %Plot heat transfer distribution
title('Verttical Heat Transfer Distribution across Bundle')
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Heat Transfer (W)')
axis([Y(1),Y(numel(Y)),0,Q_sum_store(1)]);
axis 'auto y'
%% Plot Gas Outlet Temperature & Heat Transfer
T_g_out_avg=zeros(size(T_l_out_store,1),1);
T_g_out_max=zeros(size(T_l_out_store,1),1);
T_g_out_min=zeros(size(T_l_out_store,1),1);
T_g_out_std=zeros(size(T_l_out_store,1),1);
P_g_out_avg=zeros(size(T_l_out_store,1),1);
P_g_out_max=zeros(size(T_l_out_store,1),1);
P_g_out_min=zeros(size(T_l_out_store,1),1);
P_g_out_std=zeros(size(T_l_out_store,1),1);
delta_P_g_avg=zeros(size(T_l_out_store,1),1);
delta_P_g_max=zeros(size(T_l_out_store,1),1);
delta_P_g_min=zeros(size(T_l_out_store,1),1);
for j1=1:size(T_g_out_store,1)
    T_g_out_avg(j1)=mean(T_g_out_store{j1});
    T_g_out_max(j1)=max(T_g_out_store{j1});
    T_g_out_min(j1)=min(T_g_out_store{j1});
    T_g_out_std(j1)=std(T_g_out_store{j1});
    P_g_out_avg(j1)=mean(P_g_out_store{j1});
    P_g_out_max(j1)=max(P_g_out_store{j1});
    P_g_out_min(j1)=min(P_g_out_store{j1});
    P_g_out_std(j1)=std(P_g_out_store{j1});
end
for k1=1:size(P_g_out_store,1)
    delta_P_g_avg(k1)=P_g_inlet(size(P_g_out_store,1)+1-k1)-P_g_out_avg(k1);
    delta_P_g_max(k1)=P_g_inlet(size(P_g_out_store,1)+1-k1)-P_g_out_min(k1);
    delta_P_g_min(k1)=P_g_inlet(size(P_g_out_store,1)+1-k1)-P_g_out_max(k1);
end
figure(9)
hold on
plot(Y(1:n),T_g_out_avg,'-k') %Plot gas outlet temperature distribution
plot(Y(1:n),T_g_out_min,'-.b') %Minimum outlet temperature at each cross section
plot(Y(1:n),T_g_out_max,'--r') %Minimum outlet temperature at each cross section
title([gas ' Outlet Temperature Distribution'])
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Temperature (\circC)')
axis([Y(1),Y(numel(Y)),0,T_g_out_max(1)]);
axis 'auto y'
legend('Mean Outlet Temperature','Minimum Temperature','Maximum Temperature','Location','Best')
figure(10)
hold on
plot(Y(1:n),delta_P_g_avg,'-k') %Plot gas pressure drop distribution
plot(Y(1:n),delta_P_g_min,'-.b') %Minimum pressure drop at each cross section
plot(Y(1:n),delta_P_g_max,'--r') %Maximum pressure drop at each cross section
title([gas ' Pressure Drop Distribution'])
xlabel('Vertical Position in Tube Bundle (m)')
ylabel('Pressure Drop (bar)')
axis([Y(1),Y(numel(Y)),0,delta_P_g_max(1)]);
axis 'auto y'
legend('Mean Pressure Drop','Minimum Pressure Drop','Maximum Pressure Drop','Location','Best')
Area_tot=sum(Area_store);
U_mean=mean(U_store);
%% Mean Outlet Conditions & Effectiveness
Cp_l_out=zeros(size(T_l_out_avg,1),1);
Cp_g_out=zeros(size(T_l_out_avg,1),1);
for l1=1:size(T_g_out_avg,1)
    [Cp_l_out(l1),Cp_g_out(l1)]=Material_prop(liquid,gas,tube_material,T_l_out_avg(l1),T_g_out_avg(l1),P_g_out_avg(l1));
end
T_l_out_mean_tot=sum((m_l_2_D(2:n+1).*Cp_l_out.*T_l_out_avg))/sum((m_l_2_D(2:n+1).*Cp_l_out));
T_g_out_mean_tot=sum((flipud(m_g_2_D(1:n)).*Cp_g_out.*T_g_out_avg))/sum((flipud(m_g_2_D(1:n)).*Cp_g_out));
P_l_out_mean_tot=sum((m_l_2_D(2:n+1).*P_l_out_avg))/m_l;
P_g_out_mean_tot=sum((flipud(m_g_2_D(1:n)).*P_g_out_avg))/m_g;
T_l_avg=(T_l_out_mean_tot+T_l_in)/2;
T_g_avg=(T_g_out_mean_tot+T_g_in)/2;
P_g_avg=(P_g_out_mean_tot+P_g_in)/2;
delta_P_l_total=P_l_in-P_l_out_mean_tot;
delta_P_g_total=P_g_in-P_g_out_mean_tot;
[Cp_l,Cp_g]=Material_prop(liquid,gas,tube_material,T_l_avg,T_g_avg,P_g_in);
C_min=min(m_g*Cp_g,m_l*Cp_l);
Q_max=C_min*(T_l_in-T_g_in);
Q_total=sum(Q_sum_store);
epsilon_3D=Q_total/Q_max;
LMTD_tot=((T_l_in-T_g_out_mean_tot)-(T_l_out_mean_tot-T_g_in))/log((T_l_in-T_g_out_mean_tot)/(T_l_out_mean_tot-T_g_in));
F_factor=Q_total/(U_mean*Area_tot*LMTD_tot);
fprintf('The effectiveness of this heat exchanger Is %4.4f.\n',epsilon_3D)
fprintf('The %s outlet temperature is %1.1f%cC.\n',liquid,T_l_out_mean_tot,char(176))
fprintf('The %s outlet temperature is %1.1f%cC.\n',gas,T_g_out_mean_tot,char(176))
fprintf('The %s pressure drop is %1.2f bar.\n',liquid,delta_P_l_total)
fprintf('The %s pressure drop is %1.4f bar.\n',gas,delta_P_g_total)
fprintf('The overall heat transfer is %1.3e W.\n',Q_total)