%This function will calculate and interpolate the thermophysical properties
%of the air based on the data for air in F. Incropera and D. De Witt,
%"Fundamentals of Heat Transfer," 5th Ed., Wiley, 1990. Table A.4, pg. 852
%based on the temperature given in Celsius and pressure in bar.
function [rho_g,Cp_g,mu_g,k_g,Pr_g] = Air_prop(T_g_avg,P_g_avg)
T_g_avg_K=T_g_avg+273.15; % [K] Converts temp. to Kelvin
rho_g=0.4354*(800/(T_g_avg_K))*(P_g_avg/1.013); %[kg/m^3] Uses ideal gas law and density at 800K and atm press.
air=csvread('air_data.csv',1,0); %Air data from Incropera
Cp_g=1000*interp1(air(:,1),air(:,2),T_g_avg_K); %[J/kg*K]
mu_g=interp1(air(:,1),air(:,3),T_g_avg_K); %[N*s/m^2]
k_g=interp1(air(:,1),air(:,4),T_g_avg_K); %[W/m*K]
Pr_g=Cp_g*mu_g/k_g; %dimensionless