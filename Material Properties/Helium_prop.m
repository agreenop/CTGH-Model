%This function will calculate and interpolate the thermophysical properties
%of the helium based on the data for helium from NIST using their RefProp
%program. For the inputs, temperature is Celsius and pressure is bar.
function [rho_g,Cp_g,mu_g,k_g,Pr_g] = Helium_prop(T_g_avg,P_g_avg) 
T_g_avg_K=T_g_avg+273.15; % [K] Converts temp. to Kelvin
atomic_mass=4.002602; %Atomic weight of Helium [AMU] or [g/mol]
rho_g=0.060965*(800/(T_g_avg_K))*(P_g_avg/1.0132); %[kg/m^3] Uses ideal gas law and density at 800K and atm press.
helium=csvread('helium_data.csv',1,0); %Air data from Incropera
Cp_g_mol=interp1(helium(:,1),helium(:,4),T_g_avg_K); %[J/mol*K]
Cp_g=Cp_g_mol/(atomic_mass/1000); %[J/kg*K]
mu_g=interp1(helium(:,1),helium(:,5),T_g_avg_K); %[Pa*s]
k_g=interp1(helium(:,1),helium(:,6),T_g_avg_K); %[W/m*K]
Pr_g=Cp_g*mu_g/k_g; %dimensionless