%This function will calculate the properties of the 'hot' coolant and
%"cool" gas.  The first time it will assume that the inlet condtion
%properties remain constant throughout the system.  The second time it will
%use the average gas temperature.
function [UA,Cp_l,Cp_g,mu_l,rho_l,u_max_app,rho_g]=heat_properties(inlet_prop,gas,liquid,tube_material,D_out,t,ST,SL,T_l_in,T_g_in,P_l_in,P_g_in,T_g,T_l,P_g,P_l,m_g_vol,i,j,i1,j1)
if inlet_prop==1 %First time, properties will be calculated at inlet temperatures and pressures
    T_l_avg=T_l_in;
    T_g_avg=T_g_in;
    P_g_avg=P_g_in;
    P_l_avg=P_l_in;
else %Every other time, properties will be calculated at average temperatures and pressures for each volume.
    T_l_avg=(T_l(i1,j1)+T_l(i,j))/2;
    T_g_avg=(T_g(i,j)+T_g(i+1,j))/2;
    P_g_avg=(P_g(i,j)+P_g(i+1,j))/2;
    P_l_avg=(P_l(i1,j1)+P_l(i,j))/2;
end
[tubes_vol,N_T,N_L,tubes,D_in,L,H,k_t,rho_t,Cp_t]=CTGH_geom(tube_material,D_out,t);
switch liquid %Liquid properties depending on type of liquid
    case 'Fluoride Salt'
        [mu_l,Cp_l,k_l,rho_l,nu_l,Pr_l] = Flibe_prop(T_l_avg);
    case 'Water'
        [mu_l,Cp_l,k_l,rho_l,nu_l,Pr_l] = Water_prop(T_l_avg);
end
switch gas %Gas properties depending on type of gas
    case 'Air'
        [rho_g,Cp_g,mu_g,k_g,Pr_g] = Air_prop(T_g_avg,P_g_avg);
end
%This next part finds UA.
Nu_l=3.66; %Nusselt number for fully developed laminar flow in a pipe
h_l=k_l*Nu_l/D_in; %Liquid heat transfer coefficient 
R_l=1/(tubes_vol*pi*D_in*L*h_l); %Liquid thermal resistance for pipes
R_t=log(D_out/D_in)/(2*pi*tubes_vol*k_t*L); %Metal thermal resistance for pipes
%Air flow across cross flow tubes.  See Khan paper.
u_app_g=m_g_vol/(rho_g*H*L); %Approach velocity of gas
u_max_app=max((ST/(ST-1))*u_app_g,(ST/(2*(sqrt(SL^2+(ST/2)^2)-1)))*u_app_g); %Max velocity of gas/ velocity of gas between tubes
Re_g=D_out*u_max_app*rho_g/(mu_g); %Reynolds number for gas based on max velocity
C2=(0.588+0.004*ST)*(0.858+0.04*ST-0.008*ST^2)^(1/SL);
Nu_df_g=C2*Re_g^(1/2)*Pr_g^(1/3)+0.001*Re_g; 
C1=(1.21+1.64*N_L^1.44)/(1.87+N_L^1.44);
Nu_g=C1*Nu_df_g; %Nusselt number for gas
%End of Khan Paper equations
h_g=k_g*Nu_g/D_out; %Gas heat transfer coefficient 
R_g=1/(tubes_vol*pi*D_out*L*h_g); %Gas thermal resistance
UA=1/(R_l+R_t+R_g); %Total UA for volume based on thermal resistances