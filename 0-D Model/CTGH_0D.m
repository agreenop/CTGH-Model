%This program will create a 0-D model of the Coiled Gas Tube Heater.  This
%is based off the 0-D Excel Spreadsheet.  February 23,2016
function CTGH_0D(THEEM_model,i)
if strcmp(THEEM_model,'0D') %Runs for 0-D model only.
    load('THEEM_Input_0D.mat');
else %Runs for optimization code.
    fname1=sprintf('Optimization_Files/Inputs/Input%d.mat',i);
    load(fname1);
end
T_g_avg=(T_g_in+T_g_out)/2; %Average gas outlet temp. [degC]
T_l_avg=(T_l_in+T_l_out)/2; %Average liquid outlet temp. [degC]
LMTD=((T_l_in-T_g_out)-(T_l_out-T_g_in))/log((T_l_in-T_g_out)/(T_l_out-T_g_in)); %Log Mean Temp. Difference
%% Thermodynamic Properties
[Cp_l,Cp_g,mu_l,k_l,rho_l,Pr_l,rho_g,mu_g,k_g,Pr_g,k_t]=Material_prop(liquid,gas,tube_material,T_l_avg,T_g_avg,P_g_in);
Vol_flow_g=m_g/rho_g; %Volumetric flow rate of gas [m^3/s]
Vol_flow_l=m_l/rho_l; %Volumetric flow rate of liquid [m^3/s]
%% CTAH Geometry and Bundle Height
D_in=D_out-2*t; 
disk_thick=.003; %Thickness of plate separating bundles
D_curve_inner=2*R_ci; %Average tube bundle inside diameter [m]
tube_col=entry*(tube_layer*2)*loops; %Number of tube columns, inlcuding heating rods & accounting for staggered arrangement
tubes_manifold=layer_num*(tube_layer-heat_rod); %Number of tubes per manifold per sub-bundle
tubes=entry*tubes_manifold*bundles; %Total number of tubes in CTGH
bank_depth=tube_col*SL*D_out+spacers*spacer_width; %Depth of tube bank based on tubes, tube spacing, and spacer gaps [m]
D_curve_outer=D_curve_inner+2*bank_depth; %Average tube bundle outside diameter [m]
D_curv_avg=(D_curve_outer+D_curve_inner)/2; %Diameter of the middle of the tube bundle [m]
H=D_out*ST*((layer_num+1)/2); %Height of sub-bundle, excluding spacer disk [m]
H_bank=(H+disk_thick)*bundles; %Height of entire tube bank, including spacer disk [m]
Area_avg=pi*D_curv_avg*(H_bank-disk_thick*bundles); %Average/Middle of bundle cross sectional area
L_tube=loops*pi*D_curv_avg; %Average length of each tube in bundle
Area_surf=pi*D_out*L_tube*tubes; %Surface area of tubes, based on outside diameter
%% Analysis of Gas Cross Flow
v_g_avg=Vol_flow_g/Area_avg; %Average velocity of gas based on bundle flow area [m/s]
v_g_max=v_g_avg*max((ST/(ST-1)),(ST/(2*(sqrt(SL^2+(ST/2)^2)-1)))); %Maximum gas velocity between tubes [m/s]
Re_g=rho_g*v_g_max*D_out/mu_g; %Gas Reynolds number
N_L_list=[1,2,3,4,5,7,10,13,16,20];
C2_list=[0.64,0.76,0.84,0.89,0.92,0.95,0.97,0.98,0.99,1];
    if tube_col<20
        C2=interp1(N_L_list,C2_list,tube_col);
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
Nu_g=C2*C1*Re_g^m*Pr_g^0.36; %Tube bank Nusselt correlation, 0.7<Pr<2000, 1000<Re_g<2*10^6 (Incropera Eq. 7.56)
h_g=k_g*Nu_g/D_out; %Gas heat transfer coefficient
%% Analysis of Salt Flow in Tubes
Nu_l=3.66; %Nusselt number for liquid assuming laminar flow
h_l=k_l*Nu_l/D_in; %Gas heat transfer coefficient
%% Overall Heat Transfer Coefficient U
R_g=1/h_g; %Gas thermal resistance based on outer diameter
R_t=log(D_out/D_in)*D_out/(2*k_t); %Metal thermal resistance for pipes based on outer diameter
R_l=(D_out/D_in)*1/h_l; %Liquid thermal resistance for pipes based on outer diameter
U=1/(R_g+R_t+R_l); %Overall heat transfer coefficient [W/m^2*K] based on outer diameter
%% Surface Area and Tube Requirements
Q_tot=m_g*Cp_g*(T_g_out-T_g_in)/10^6; %Gas thermal power [MW]
A_ideal=Q_tot*10^6/(U*LMTD); %Ideal (F=1) surface area based on outer diameter[m^2]
F=A_ideal/Area_surf; %Reguired F factor needed to obtain heat transfer, Q_tot
L_ideal=A_ideal/(pi*D_out); %Ideal total tube length
tubes_ideal=round(L_ideal/L_tube,0); %Estimated ideal number of tubes, assuming average tube length
C_min=min(m_g*Cp_g,m_l*Cp_l);
Q_max=C_min*(T_l_in-T_g_in)/10^6;
epsilon=Q_tot/Q_max;
%% Pressure drop across CTGH
% Gas Pressure Drop
f=0.4; %Gas friction factor from Incropera Fig. 7.14
chi=1.00; %Correction factor from Incropera Fig. 7.14
deltaP_tube_g=1/2*f*chi*rho_g*v_g_max^2; %Pressure drop across tube column [Pa]
deltaP_g=deltaP_tube_g*tube_col*10^-5; %Pressure drop across entire tube bundle [bar]

%Salt Pressure Drop
Flow_area=tubes*pi/4*D_in^2; %Sum total of flow area inside all tubes
v_l=m_l/(rho_l*Flow_area); %Average velocity of salt
Re_l=rho_l*v_l*D_in/mu_l; %Salt Reynolds number
f_l=64/Re_l; %Salt friction factor (assuming laminar flow)
deltaP_l=1/2*f_l*rho_l*L_tube*v_l^2/D_in*10^-5; %Salt pressure drop across bundle [bar]
%% Save variables to output files
if strcmp(THEEM_model,'0D') %Runs for 0-D model only.
    save('0-D Model/THEEM_Output_0D.mat');
else %Runs for Optimization Code only
    fname2=sprintf('Optimization Program/Optimization_Files/Outputs/Output%d.mat',i);
    save(fname2,'tubes','D_curve_outer','H_bank','Area_surf','v_g_max','Re_g','U','A_ideal','F','deltaP_g','deltaP_l','bank_depth');
    range_output=sprintf('I%d:T%d',i+1,i+1);
    B=[tubes,D_curve_outer,H_bank,Area_surf,v_g_max,Re_g,U,A_ideal,F,deltaP_g,deltaP_l,bank_depth];
    xlswrite('Optimization Program/Optimization_Files/Optimization_Results.xlsx',B,range_output);
end