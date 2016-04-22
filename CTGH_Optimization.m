% This program performs a sensitivity analysis on the CTGH 0-D code.  It
% changes the geometry of the CTGH in order to minimize the necessary
% effectiveness for the CTGH along with the salt and liquid pressure drops.
%Default Values for Mk1 model:
%These are the inputs that will change for the sensitivity analysis:
% D_out=0.25*0.0254;
% entry=4;
% SL=1.45;
% ST=1.256;
% t=.035*.0254;
%% Initialize the Excel Spreadsheet and assign static inputs
clc;clear;
delete('Optimization_Files/Optimization_Results.xlsx') %Delete previous results
%These inputs will not be changed by the optimization code:
m_g=418.5;
m_l=480.2;
P_g_in=18.76;
P_l_in=3.3;
T_g_in=418.6;
T_l_in=700;
T={'Run','Outside Diameter [in]','Thickness [in]','Longitudinal Pitch',...
   'Transverse Pitch','Number of Inlets','Total Tubes',...
   'Outer Diameter of Curvature [m]','Height of Tube Bank [m]',...
   'Total Surface Area [m^2]','Max Air Velocity [m/s]',...
   'Air Reynolds Number','U Coefficient [W/m^2*K]',...
   'Ideal Surface Area [m^2]','Min F Factor','Air Pressure Drop [bar]',...
   'Liquid Pressure Drop [bar]','Tube Bank Depth [m]'}; %Column labels in the Excel Spreadsheet
xlswrite('Optimization_Files/Optimization_Results.xlsx',T,'A1:R1') 
%% Tube outside diameters and thicknesses obtained from McMaster-Carr website
for i=1:3000 %Run 3000 different scenarios
Diameters=[1/16:1/32:5/16,3/8,7/16,1/2,5/8,11/16];%,3/4,7/8,1,1.25,1.375,1.5,1.75,2,2.125,2.25,2.375,2.5,3]; %Diameters in inches
D_out=Diameters(randi(numel(Diameters))); %Randomly selects a tube outisde diameter for the design
switch D_out %Assign tube thicknesses based on diameter in inches
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
%     case 3/4
%         thickness=[0.035,0.049,0.065,0.083,0.12];
%     case 7/8
%         thickness=[0.035,0.049,0.065,0.12];
%     case 1
%         thickness=[0.035,0.049,0.065,0.083,0.12];
%     case 1.25
%         thickness=[0.035,0.049,0.065,0.083,0.12];
%     case 1.375
%         thickness=[0.035,0.049,0.065];
%     case 1.5
%         thickness=[0.049,0.065,0.095,0.12];
%     case 1.75
%         thickness=0.065;
%     case 2
%         thickness=[0.065,0.12];
%     case 2.125
%         thickness=0.065;
%     case 2.25
%         thickness=0.065;
%     case 2.375
%         thickness=0.065;
%     case 2.5
%         thickness=0.065;
%     case 3
%         thickness=0.065;
end
t=thickness(randi(numel(thickness))); %Randomly picks one of the thicknesses associated with the chosen diameter
D_out=D_out*0.0254; %Converts diameters to meters
t=t*0.0254; %Converts thickeness to meters
%% Randomly choose other inlet parameters
%These other inlet parameters are restricted so that the maximum depth of
%the tube bundle is 0.75 m (without sloping tubes).  The maximum cap is
%calculated at 0.0225 assuming that there are always 2 spacers, the spacer
%widths are constant (0.038 m), each tube loops 3 times around the CTGH,
%and each tube bank from each manifold pipe has 5 tubes per row.
cap=1; %Start loop
while cap > 0.0375 
entry=randi([2,6]);
SL=(2.0-1.0)*rand+1.0; %Generates a random value for SL between 1.0 & 2.0
if SL<1.25
    ST=(2.0-1.5)*rand+1.5; %If SL is less than 1.25, ST must be between 1.5 and 2.0
else
    ST=(2.0-1.25)*rand+1.25; %Otherwise ST is between 1.25 and 2.0
end
cap=entry*SL*D_out;
end
fname=sprintf('Optimization_Files/Inputs/Input%d.mat',i);
save(fname,'m_g','m_l','P_g_in','P_l_in','T_g_in','T_l_in','D_out','t','SL','ST','entry');
CTGH_0D_calculations(i)
range_input=sprintf('A%d:F%d',i+1,i+1);
A=[i,D_out/0.0254,t/0.0254,SL,ST,entry];
xlswrite('Optimization_Files/Optimization_Results.xlsx',A,range_input);
end
