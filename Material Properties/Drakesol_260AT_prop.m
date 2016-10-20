%This function will calculate the thermophysical properties of Drakesol 
%260AT mineral oil at a given temperature in Celsius and given pressure in bar.  
%It is based on Scarlat and Zweibaum, "Temperature-Dependent Thermophysical Properties
%for Fluoride Salts and Simulant Fluids", CIET-DESIGN-001-03, Rev. 03, 05-2-2013
function [mu_l,Cp_l,k_l,rho_l,Pr_l] = Drakesol_260AT_prop(T_l_avg) %Temp is in Celsius
nu_l=8.403*10^-6*exp(-0.0187*T_l_avg); %[m^2/s] Kinematic Viscosity
Cp_l=2000; %J/(kg*K)
k_l=0.098; % W/(m*K)
rho_l=805-0.5*T_l_avg; % kg/(m^3)
mu_l=rho_l*nu_l; %[kg/(m*s)] Dynamic viscosity
Pr_l=Cp_l*mu_l/k_l; %dimensionless