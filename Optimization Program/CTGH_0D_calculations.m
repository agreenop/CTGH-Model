%This program will create a 0-D model of the Coiled Gas Tube Heater.  This
%is based off the 0-D Excel Spreadsheet.  February 23,2016
function CTGH_0D_calculations(i)
fname1=sprintf('Optimization_Files/Inputs/Input%d.mat',i);
load(fname1);
T_g_avg=(T_g_in+T_g_out)/2; %Average gas outlet temp. [degC]
T_l_avg=(T_l_in+T_l_out)/2; %Average liquid outlet temp. [degC]
LMTD=((T_l_in-T_g_out)-(T_l_out-T_g_in))/log((T_l_in-T_g_out)/(T_l_out-T_g_in)); %Log Mean Temp. Difference
%% Thermodynamic Properties
%Average air properties
[rho_g,cp_g,mu_g,k_g,Pr_g] = Air_prop(T_g_avg,P_g_in); %Avg. Heater Air Properties
Vol_flow_g=m_g/rho_g; %Volumetric flow rate of gas [m^3/s]

%Average flibe properties
[mu_l,cp_l,k_l,rho_l,nu_l,Pr_l] = Flibe_prop(T_l_avg);
Vol_flow_l=m_l/rho_l; %Volumetric flow rate of liquid [m^3/s]

%Metal (316 SS) properties
k_t=13.40; %316 SS thermal conductivity [W/m*K]
rho_t=8238; %316 SS density [kg/m^3]
Cp_t=468; %316 SS specific heat [J/kg*K]
%% CTAH Geometry and Bundle Height
D_in=D_out-2*t; 
bundles=36; %Number of sub-bundles
disk_thick=.003; %Thickness of plate separating bundles
% row_num=34; %Number of tube layers per sub-bundle
heat_rod=1/2; %Number of heater rods in each tube row (1 every 2 rows)
loops=3; %Number of times tubes loop around bundle
D_curve_inner=1.324; %Average tube bundle inside diameter [m]
spacers=2; %Number of spacer gaps in each sub-bundle (allows for air mixing)
spacer_width=0.038; %Width of each spacer gap based on tie rod diameter [m]
% tube_row=3; %Number of tubes per row per manifold pipe
tube_col=entry*(tube_row*2)*loops; %Number of tube columns, inlcuding heating rods & accounting for staggered arrangement
tubes=bundles*row_num*entry*(tube_row-heat_rod); %Total number of tubes in CTGH
bank_depth=tube_col*SL*D_out+spacers*spacer_width; %Depth of tube bank based on tubes, tube spacing, and spacer gaps [m]
D_curve_outer=D_curve_inner+2*bank_depth; %Average tube bundle outside diameter [m]
D_curv_avg=(D_curve_outer+D_curve_inner)/2; %Diameter of the middle of the tube bundle [m]
H=D_out*ST*((row_num+1)/2); %Height of sub-bundle, excluding spacer disk [m]
H_bank=(H+disk_thick)*bundles; %Height of entire tube bank, including spacer disk [m]
Area_avg=pi*D_curv_avg*(H_bank-disk_thick*bundles); %Average/Middle of bundle cross sectional area
L_tube=loops*pi*D_curv_avg; %Average length of each tube in bundle
Area_surf=pi*D_out*L_tube*tubes; %Surface area of tubes, based on outside diameter
%% Analysis of Air Cross Flow
v_g_avg=Vol_flow_g/Area_avg; %Average velocity of gas based on bundle flow area [m/s]
v_g_max=v_g_avg*ST/(ST-1); %Maximum gas velocity between tubes [m/s]
Re_g=rho_g*v_g_max*D_out/mu_g; %Gas Reynolds number
SL_list=[0.90,1,1.125,1.25,1.5,2,3];%Table 7.5 values from Incropera
ST_list=[1.25,1.5,2.0,3.0];%Table 7.5 values from Incropera
C1_table=[NaN,NaN,NaN,0.5180,0.4510,0.4040,0.3100;
          NaN,0.4970,0.5010,0.5050,0.4600,0.4160,0.3560;
          0.4460,0.4602,0.4780,0.5190,0.4520,0.4820,0.4400;
          0.4010,0.4530,0.5180,0.5220,0.4880,0.4490,0.4280];%Table 7.5 values from Incropera
m_table=[NaN,NaN,NaN,0.556,0.568,0.572,0.592;
          NaN,0.558,0.556,0.554,0.562,0.568,0.580;
          0.571,0.5683,0.565,0.556,0.568,0.556,0.562;
          0.581,0.5717,0.560,0.562,0.568,0.570,0.574];%Table 7.5 values from Incropera
C1=interp2(SL_list,ST_list,C1_table,SL,ST); %Interpolated from table values
m=interp2(SL_list,ST_list,m_table,SL,ST); %Interpolated from table values
Nu_g=C1*Re_g^m; %Gas Nusselt Number calculation from Incropera
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
Q_tot=m_g*cp_g*(T_g_out-T_g_in)/10^6; %Air thermal power [MW]
% F=0.90; % Effectiveness factor (assumed)
A_ideal=Q_tot*10^6/(U*LMTD); %Ideal (F=1) surface area based on outer diameter[m^2]
F=A_ideal/Area_surf; %Reguired F factor needed to obtain heat transfer, Q_tot
L_ideal=A_ideal/(pi*D_out); %Ideal total tube length
tubes_ideal=round(L_ideal/L_tube,0); %Estimated ideal number of tubes, assuming average tube length
C_min=min(m_g*cp_g,m_l*cp_l);
Q_max=C_min*(T_l_in-T_g_in)/10^6;
e1=Q_tot/Q_max;
%% Pressure drop across CTGH
% Air Pressure Drop
f=0.4; %Air friction factor from Incropera Fig. 7.14
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
fname2=sprintf('Optimization Program/Optimization_Files/Outputs/Output%d.mat',i);
save(fname2,'tubes','D_curve_outer','H_bank','Area_surf','v_g_max','Re_g','U','A_ideal','F','deltaP_g','deltaP_l','bank_depth');
range_output=sprintf('I%d:T%d',i+1,i+1);
B=[tubes,D_curve_outer,H_bank,Area_surf,v_g_max,Re_g,U,A_ideal,F,deltaP_g,deltaP_l,bank_depth];
xlswrite('Optimization Program/Optimization_Files/Optimization_Results.xlsx',B,range_output);