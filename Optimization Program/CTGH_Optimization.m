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
delete('0-D Model & Optimization Program/Optimization_Files/Optimization_Results.xlsx') %Delete previous results
%These inputs will not be changed by the optimization code:
m_g=418.5;
m_l=480.2;
P_g_in=18.76;
P_l_in=3.3;
T_g_in=418.6;
T_l_in=700;
T_l_out=600; %Liquid Outlet temperature [degC] 
T_g_out=670; %Gas outlet temperature [degC]
T={'Run','Outside Diameter [in]','Thickness [in]','Longitudinal Pitch',...
   'Transverse Pitch','Number of Inlets','Tube Layers per sub-bundle',...
   'Number of Tube per Layer per Manifold','Total Tubes',...
   'Outer Diameter of Curvature [m]','Height of Tube Bank [m]',...
   'Total Surface Area [m^2]','Max Air Velocity [m/s]',...
   'Air Reynolds Number','U Coefficient [W/m^2*K]',...
   'Ideal Surface Area [m^2]','Min F Factor','Air Pressure Drop [bar]',...
   'Liquid Pressure Drop [bar]','Tube Bank Depth [m]'}; %Column labels in the Excel Spreadsheet
xlswrite('Optimization_Files/Optimization_Results.xlsx',T,'A1:T1') 
%% Tube outside diameters and thicknesses obtained from McMaster-Carr website
for i=1:3000 %Run 3000 different scenarios
Diameters=[3/16:1/32:1/4];%5/16];%,3/8,7/16,1/2,5/8,11/16];%,3/4,7/8,1,1.25,1.375,1.5,1.75,2,2.125,2.25,2.375,2.5,3]; %Diameters in inches
D_out=Diameters(randi(numel(Diameters))); %Randomly selects a tube outisde diameter for the design
switch D_out %Assign tube thicknesses based on diameter in inches
%     case 1/16
%         thickness=[0.02,0.028];
%     case 3/32
%         thickness=0.02;
%     case 1/8
%         thickness=[0.02,0.028,0.035,0.049];
%     case 5/32
%         thickness=[0.035,0.049];
    case 3/16
        thickness=[0.02,0.028,0.035,0.049];
    case 7/32
        thickness=0.035;
    case 1/4
        thickness=[0.01,0.016,0.02,0.028,0.035,0.049,0.065,0.083];
%     case 9/32
%         thickness=0.028;
%     case 5/16
%         thickness=[0.016,0.028,0.035,0.049,0.065,0.083];
%     case 3/8
%         thickness=[0.01,0.02,0.028,0.035,0.049,0.065,0.083];
%     case 7/16
%         thickness=0.049;
%     case 1/2
%         thickness=[0.02,0.028,0.035,0.049,0.065,0.083,0.12];
%     case 5/8
%         thickness=[0.02,0.028,0.035,0.049,0.065,0.083,0.12];
%     case 11/16
%         thickness=0.065;
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
%% Randomly choose other inlet parameters within constraints
%Tube Bundle dimension Constraints
width=1; %Start loop
width_limit=0.09572; %Corresponds to maximum allowable bank width of 0.6502 m (bank outer diameter=2.6245 m)
height_limit=0.2192; %Corresponds to maximum allowable bank height of 8 m

%Minimum and maximum values for each parameter
SL_min=1.25; %In order for SL to be less than 1.25, ST>1.5.  In order to have the SD=ST, SL=ST*cosd(30). 
              %At SL=1.25, ST=1.4434.  So, it is impossible for SL<1.25.
SL_tot_max=2.0; %Absolute max value of SL
entry_min=2; %Minimum # of manifold inlets
entry_max=4; %Maximum # of manifold inlets
tube_row_min=3; %Minimum # of tubes per layer per manifold inlet
tube_row_tot_max=6; %Absolute maximum allowable number of tubes in each row
row_num_max=42;

while width > width_limit 
entry=randi([entry_min,entry_max]);
% entry=4;
tube_row_max=floor(width_limit/(1.25*entry*D_out)); %Maximum # of tubes with given parameters
if tube_row_max>tube_row_tot_max
    tube_row_max=tube_row_tot_max;
end
if tube_row_max>tube_row_min %If tube_row_max<tube_row_min, chooses new # of inlet manifolds
tube_row=randi([tube_row_min,tube_row_max]);
% tube_row=3;
% row_num=40;
SL_max_1=width_limit/(entry*D_out*tube_row); %Restrict SL on width of bundle
SL_max_2=2*height_limit*cosd(30)/((row_num_max+1)*D_out); %Restrict SL on height of bundle with a max # of tube layers
SL_max=min(SL_max_1,SL_max_2); %Max SL possible given parameters
if SL_max>SL_tot_max 
    SL_max=SL_tot_max; %Make sure SL stays below absolute maximum SL
end
if SL_max>SL_min %If SL_max>SL_min, chooses new parameters with same tube diameter & thickness.
    SL=(SL_max-SL_min)*rand+SL_min;
    ST=SL/cosd(30); %In order to have the SD=ST, SL=ST*cosd(30). Forms equilateral triangle.
    row_num_min=ceil(100/tube_row); %Have at least 100 tubes per manifold pipe in each sub-bundle
    row_num=randi([row_num_min,row_num_max]);
    width=entry*SL*D_out*tube_row;
end
end
end
fname=sprintf('0-D Model & Optimization Program/Optimization_Files/Inputs/Input%d.mat',i);
save(fname,'m_g','m_l','P_g_in','P_l_in','T_g_in','T_l_in','T_l_out','T_g_out','D_out','t','SL','ST','entry','row_num','tube_row');
range_input=sprintf('A%d:H%d',i+1,i+1);
A=[i,D_out/0.0254,t/0.0254,SL,ST,entry,row_num,tube_row];
xlswrite('0-D Model & Optimization Program/Optimization_Files/Optimization_Results.xlsx',A,range_input);
CTGH_0D_calculations(i)
end
