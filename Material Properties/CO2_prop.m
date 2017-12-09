%This function will calculate and interpolate the thermophysical properties
%of the air based on the data for air in F. Incropera and D. De Witt,
%"Fundamentals of Heat Transfer," 5th Ed., Wiley, 1990. Table A.4, pg. 852
%based on the temperature given in Celsius and pressure in bar.
function [rho_g,Cp_g,mu_g,k_g,Pr_g] = CO2_prop(T_g_avg,P_g_avg)
T_g_avg_K=T_g_avg+273.15; % [K] Converts temp. to Kelvin
CO2=csvread('CO2_data.csv',1,0); %Air data from Incropera
rho_g=interp1(CO2(:,1),CO2(:,3),800)*((800+273.15)/(T_g_avg_K))*(P_g_avg/interp1(CO2(:,1),CO2(:,2),800)); %[kg/m^3] Uses ideal gas law and density at 800C and given pressure.
Cp_g=interp1(CO2(:,1),CO2(:,4),T_g_avg_K)*10^3; %[J/kg*K]
mu_g=interp1(CO2(:,1),CO2(:,6),T_g_avg_K)*10^-6; %[N*s/m^2]
k_g=interp1(CO2(:,1),CO2(:,5),T_g_avg_K)*10^-3; %[W/m*K]
Pr_g=Cp_g*mu_g/k_g; %dimensionless