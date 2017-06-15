%This function will calculate the properties of the 'hot' coolant and
%"cool" gas.  The first time it will assume that the inlet condtion
%properties remain constant throughout the system.  The second time it will
%use the average gas temperature.
function [UA,Cp_l,Cp_g,rho_l,u_max_app,rho_g,Re_g,h_g,Area,Re_l,f_l,De_l,h_l,Nu_l,Nu_g]=heat_properties(inlet_prop,T_g,T_l,P_g,m_g_vol,i,j,i1,j1,m_l_t,T_s_out_matrix,THEEM_model)
if strcmp(THEEM_model, '3D')
    load('THEEM_Input_3D.mat');
else
    load('THEEM_Input_2D.mat');
end
if inlet_prop==1 %First time, properties will be calculated at inlet temperatures and pressures
    T_l_avg=T_l_in;
    T_g_avg=T_g_in;
    P_g_avg=P_g_in;
    [~,~,~,~,~,~,~,~,~,Pr_s,~]=Material_prop(liquid,gas,tube_material,T_l_avg,T_l_avg,P_g_avg); %Estimates gas Prandtl number at tube surface
else %Every other time, properties will be calculated at average temperatures and pressures for each volume.
    T_l_avg=(T_l(i1,j1)+T_l(i,j))/2;
    T_g_avg=(T_g(i,j)+T_g(i+1,j))/2;
    P_g_avg=(P_g(i,j)+P_g(i+1,j))/2;
    T_s_out=T_s_out_matrix(i,j);
    [~,~,~,~,~,~,~,~,~,Pr_s,~]=Material_prop(liquid,gas,tube_material,T_l_avg,T_s_out,P_g_avg); %Estimates gas Prandtl number at tube surface
end
[L,R_curv,H,tubes_vol,~,N_L,~,D_in]=CTGH_geom(THEEM_model,i);
[Cp_l,Cp_g,mu_l,k_l,rho_l,Pr_l,rho_g,mu_g,k_g,Pr_g,k_t]=Material_prop(liquid,gas,tube_material,T_l_avg,T_g_avg,P_g_avg);
% [~,~,~,~,~,~,~,~,~,Pr_s,~]=Material_prop(liquid,gas,tube_material,T_l_avg,T_l_avg,P_g_avg); %Estimates gas Prandtl number at tube surface
%This next part finds UA.
Re_l=4*m_l_t/(pi*D_in*mu_l); %Reynolds number for liquid
De_l=Re_l*sqrt(D_in/(2*R_curv)); %Dean number for liquid through a curved pipe
Re_c=2300*(1+12*sqrt(D_in/(2*R_curv))); %Critical reynolds number for a cruved pipe
if Re_l<=Re_c %Laminar flow
    Nu_l=((3.657+4.343/(1+957/(De_l^2*Pr_l))^2)^3+1.158*(De_l/(1+.477/Pr_l))^(3/2))^(1/3); %Nusselt number for fully developed laminar flow in a curved pipe with uniform wall temp (Manlapaz/Churchill)
    if De_l<=11.6
        f_l=64/Re_l; %Friction Factor for laminar flow in pipe
    elseif De_l>11.6 && De_l<=2000
         f_l=(64/Re_l)/(1-(1-(11.6/De_l)^0.45)^2.2);
    elseif De_l>2000
        f_l=7.0144*sqrt(De_l)/Re_l;
    end
elseif strcmp(liquid,'Sodium')==1 %Turbulent and liquid metal
   f_l=(0.790*log(Re_l)-1.64)^(-2); %Friction factor for turbulent flow for smooth straight pipe
   Nu_l=5.0+0.025*(Re_l*Pr_l)^0.8; %Nusselt number for straight pipe with turbulent liquid metal
elseif Re_l>Re_c %Transition Zone Flow/ Turbulent Flow
   f_l=.336*(D_in/(2*R_curv))^0.1*Re_l^-0.2; %Friction factor for turbulent flow for smooth curved pipe
   Nu_l=((f_l/8)*(Re_l-1000)*Pr_l)/(1+12.7*(f_l/8)^0.5*(Pr_l^(2/3)-1)); %Nusselt number for straight pipe with turbulence
end    
h_l=k_l*Nu_l/D_in; %Liquid heat transfer coefficient 
R_l=1/(tubes_vol*pi*D_in*L*h_l); %Liquid thermal resistance for pipes
R_t=log(D_out/D_in)/(2*pi*tubes_vol*k_t*L); %Metal thermal resistance for pipes
%Air flow across cross flow tubes. 
u_app_g=m_g_vol/(rho_g*H*L); %Approach velocity of gas
u_max_app=max((ST/(ST-1))*u_app_g,(ST/(2*(sqrt(SL^2+(ST/2)^2)-1)))*u_app_g); %Max velocity of gas/ velocity of gas between tubes
Re_g=D_out*u_max_app*rho_g/(mu_g); %Reynolds number for gas based on max velocity
N_L_list=[1,2,3,4,5,7,10,13,16,20];
C2_list=[0.64,0.76,0.84,0.89,0.92,0.95,0.97,0.98,0.99,1];
tube_count=2*N_L*i;
    if tube_count<20
        C2=interp1(N_L_list,C2_list,tube_count);
    else
        C2=1;
    end
if Re_g<2*10^5
    if ST/SL <2
        m=0.60;
        C1=0.35*(ST/SL)^(1/5);
    else
        m=0.60;
        C1=0.022;
    end
else
    m=0.84;
    C1=0.022;
end
if inlet_prop==1
    Nu_g=C2*C1*Re_g^m*Pr_g^0.36; %Tube bank Nusselt correlation, 0.7<Pr<2000, 1000<Re_g<2*10^6 (Incropera Eq. 7.56)
else
    Nu_g=C2*C1*Re_g^m*Pr_g^(0.36)*(Pr_g/Pr_s)^(1/4);
end    
h_g=k_g*Nu_g/D_out; %Gas heat transfer coefficient 
R_g=1/(tubes_vol*pi*D_out*L*h_g); %Gas thermal resistance
UA=1/(R_l+R_t+R_g); %Total UA for volume based on thermal resistances
Area=tubes_vol*pi*D_out*L; %Outer surface area of tubes in volume