%This function establishes the inlet conditions for the CTAH, such as
%pressure, temperature and mass flow rate.
function [T_g_in,T_l_in,P_g_in,P_l_in,m_g,m_l]=inlet_cond
T_g_in=418.6; %Gas inlet temperature in degrees Celsius
T_l_in=700; %Liquid inlet temperature in degrees Celsius
P_g_in=18.76; %Gas inlet pressure in bar
P_l_in=2; %Liquid inlet pressure in bar
m_g=418.5; %Gas mass flow rate in kg/s
m_l=480.2; %Liquid mass flow rate in kg/s
