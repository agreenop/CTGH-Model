%This function will find the thermophysical properties of the liquid, gas,
%and the tube material based on the temperature of the liquid and the
%temperature and pressure of the gas. Input temperatures are in deg C and
%input pressure is in bar.
function [Cp_l,Cp_g,mu_l,k_l,rho_l,Pr_l,rho_g,mu_g,k_g,Pr_g,k_t,rho_t,Cp_t]=Material_prop(liquid,gas,tube_material,T_l_avg,T_g_avg,P_g_avg)
%% Liquid Properties
switch liquid %Liquid properties depending on type of liquid
    case 'FLiBe Salt'
        [mu_l,Cp_l,k_l,rho_l,Pr_l] = Flibe_prop(T_l_avg);
    case 'Water'
        [mu_l,Cp_l,k_l,rho_l,Pr_l] = Water_prop(T_l_avg);
    case 'Drakesol 260AT'
        [mu_l,Cp_l,k_l,rho_l,Pr_l] = Drakesol_260AT_prop(T_l_avg);
    case 'FLiNaK Salt'
        [mu_l,Cp_l,k_l,rho_l,Pr_l] = Flinak_prop(T_l_avg);
    case 'Solar Salt (NaNO3-KNO3)'
        [mu_l,Cp_l,k_l,rho_l,Pr_l] = Solar_Nitrate_Salt_prop(T_l_avg);
    case 'Sodium'    
        [mu_l,Cp_l,k_l,rho_l,Pr_l] = Sodium_prop(T_l_avg);
end
%% Gas Properties
switch gas %Gas properties depending on type of gas
    case 'Air'
        [rho_g,Cp_g,mu_g,k_g,Pr_g] = Air_prop(T_g_avg,P_g_avg);
    case 'Helium'
        [rho_g,Cp_g,mu_g,k_g,Pr_g] = Helium_prop(T_g_avg,P_g_avg);
    case 'Supercritical CO2'
        [rho_g,Cp_g,mu_g,k_g,Pr_g] = CO2_prop(T_g_avg,P_g_avg);
end
%% Tube Material Properties
switch tube_material
    case '316 Stainless Steel'
        [k_t,rho_t,Cp_t] = SS_316_prop;
%         k_t=13.40; %316 SS thermal conductivity [W/m*K]
%         rho_t=8238; %316 SS density [kg/m^3]
%         Cp_t=468; %316 SS specific heat [J/kg*K]
end