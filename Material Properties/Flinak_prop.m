%This function will calculate the thermophysical properties of FLiNaK salt
%at a given temperature in Celsius and given pressure in bar. All
%correlations for FLiNaK (LiF - NaF - KF) obtained from 'Engineering
%Database of Liquid Salt Thermophysical and Thermochemical Properties' by
%Manohar S. Sohal, Matthias A. Ebner, Piyush Sabharwall, and Phil Sharp,
%March 2010,  INL/ETX-10-18297.
function [mu_l,Cp_l,k_l,rho_l,Pr_l] = Flinak_prop(T_l_avg) %Temp is in Celsius
T_l_avg_K=T_l_avg+273.15; %Converts temperature from Celsius to Kelvin
mu_l=2.487*10^(-5)*exp(4478.62/T_l_avg_K); %[kg/m·s]
Cp_l=(976.78+1.0634*T_l_avg_K); %[J/kg·K]
k_l=0.36+5.6*10^(-4)*T_l_avg_K; %[W/m·K]
rho_l=2729.3-0.73*T_l_avg_K; %[kg/m^3]
Pr_l=Cp_l*mu_l/k_l; %dimensionless