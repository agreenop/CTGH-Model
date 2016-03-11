% This program performs a sensitivity analysis on the CTGH 0-D code.  It
% changes the geometry of the CTGH in order to minimize the necessary
% effectiveness for the CTGH along with the salt and liquid pressure drops.
%Default Values for Mk1 model:
%These are the inputs that will change for the sensitivity analysis
% D_out=0.25*0.0254;
% entry=4;
% SL=1.45;
% ST=1.256;
% t=.035*.0254;
%% Tube outside diameters and thicknesses obtained from McMaster-Carr website
clc;clear;
%These inputs will not be changed by the sensitivity analysis:
m_g=418.5;
m_l=480.2;
P_g_in=18.76;
P_l_in=3.3;
T_g_in=418.6;
T_l_in=700;
% for i=1:3000
Diameters=[1/16:1/32:5/16,3/8,7/16,1/2,5/8,11/16,3/4,7/8,1,1.25,1.375,1.5,1.75,2,2.125,2.25,2.375,2.5,3]; %Diameters in inches
D_out=Diameters(randi(numel(Diameters))); %Randomly selects an outisde diameter for the design
switch D_out %Assign thicknesses based on diameter in inches
    case 1/16
        thickness=[0.02,0.028];
    case 3/32
        thickness=0.02;
    case 1/8
        thickness=[0.02,0.028,0.035,0.049];
    case 5/32
        thickness=[0.035,0.049];
    case 3/16
        thickness=[0.02,0.028,0.035,0.049];
    case 7/32
        thickness=0.035;
    case 1/4
        thickness=[0.01,0.016,0.02,0.028,0.035,0.049,0.065,0.083];
    case 9/32
        thickness=0.028;
    case 5/16
        thickness=[0.016,0.028,0.035,0.049,0.065,0.083];
    case 3/8
        thickness=[0.01,0.02,0.028,0.035,0.049,0.065,0.083];
    case 7/16
        thickness=0.049;
    case 1/2
        thickness=[0.02,0.028,0.035,0.049,0.065,0.083,0.12];
    case 5/8
        thickness=[0.02,0.028,0.035,0.049,0.065,0.083,0.12];
    case 11/16
        thickness=0.065;
    case 3/4
        thickness=[0.035,0.049,0.065,0.083,0.12];
    case 7/8
        thickness=[0.035,0.049,0.065,0.12];
    case 1
        thickness=[0.035,0.049,0.065,0.083,0.12];
    case 1.25
        thickness=[0.035,0.049,0.065,0.083,0.12];
    case 1.375
        thickness=[0.035,0.049,0.065];
    case 1.5
        thickness=[0.049,0.065,0.095,0.12];
    case 1.75
        thickness=0.065;
    case 2
        thickness=[0.065,0.12];
    case 2.125
        thickness=0.065;
    case 2.25
        thickness=0.065;
    case 2.375
        thickness=0.065;
    case 2.5
        thickness=0.065;
    case 3
        thickness=0.065;
end
t=thickness(randi(numel(thickness))); %Randomly picks one of the thicknesses associated with the chosen diameter
D_out=D_out*0.0254; %Converts diameters to meters
t=t*0.0254; %Converts thickeness to meters
%% Randomly choose other inlet parameters
entry=randi(8);
% SL=1.45;
SL=(3-0.6)*rand+0.6; %Generates a random value for SL between 0.6 & 3.0
if SL<1
    ST=(3-1.6)*rand+1.6; %If SL is less than 1, ST must be greater 1.6
else
    ST=(3-1.25)*rand+1.25; %Otherwise ST is between 1.25 and 3.0
end
i=1;
fname=sprintf('Optimization_Files/Inputs/Input%d.mat',i);
save(fname,'m_g','m_l','P_g_in','P_l_in','T_g_in','T_l_in','D_out','t','SL','ST','entry');
CTGH_0D_calculations(i)
% end