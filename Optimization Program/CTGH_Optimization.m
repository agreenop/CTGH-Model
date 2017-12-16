% This program performs a sensitivity analysis on the CTGH 0-D code.  It
% changes the geometry of the CTGH in order to minimize the necessary
% effectiveness for the CTGH along with the salt and liquid pressure drops.
%% Initialize the Excel Spreadsheet and assign static inputs
clc;clear;
results_location=sprintf('Optimization Program/Optimization_Files/Optimization_Results.xlsx');%Store file name as variable
delete(results_location); %Delete previous results
%These inputs will not be changed by the optimization code:
% load('THEEM_Input_Optimization.mat');
m_g=1360.5;
m_l=1267;
P_g_in=199.5;
P_l_in=3;
T_g_in=367;
T_l_in=528;
T_l_out=373; %Liquid Outlet temperature [degC] 
T_g_out=516.6; %Gas outlet temperature [degC]
gas="Supercritical CO2";
liquid="Sodium";
heat_rod=0.5;
tube_material="316 Stainless Steel";
tube_slope=0.0030;
THEEM_model='Optimization';
spacers=2;
spacer_width=0.0380;
disk_thick=.003;
T={'Run','Outside Diameter [in]','Thickness [in]','Longitudinal Pitch',...
   'Transverse Pitch','Number of Inlets','Tube Layers per sub-bundle',...
   'Number of Tubes per Layer per Manifold','Inner Radius of Bundle [m]',...
   'Number of Loops', 'Number of Sub-Bundles','Total Tubes',...
   'Outer Diameter of Curvature [m]','Height of Tube Bank [m]',...
   'Total Surface Area [m^2]','Max Air Velocity [m/s]',...
   'Air Reynolds Number','U Coefficient [W/m^2*K]',...
   'Ideal Surface Area [m^2]','Min F Factor','Air Pressure Drop [bar]',...
   'Liquid Pressure Drop [bar]','Tube Bank Depth [m]'}; %Column labels in the Excel Spreadsheet
xlswrite(results_location,T,'A1:W1')
%% Tube outside diameters and thicknesses obtained from McMaster-Carr website
for run=1:3000 %Run 3000 different scenarios
Diameters=[3/16:1/32:5/16,3/8,7/16,1/2]; %,5/8,11/16];%,3/4,7/8,1,1.25,1.375,1.5,1.75,2,2.125,2.25,2.375,2.5,3]; %Diameters in inches
D_out_in=Diameters(randi(numel(Diameters))); %Randomly selects a tube outisde diameter for the design
switch D_out_in %Assign tube thicknesses based on diameter in inches
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
t_in=thickness(randi(numel(thickness))); %Randomly picks one of the thicknesses associated with the chosen diameter (inches)
D_out=D_out_in*0.0254; %Converts diameters to meters
t=t_in*0.0254; %Converts thickeness to meters
%% Choose inner radius of tube bundle
R_ci_min=0.25; %Minimum inside radius of bundle is 0.25 m 
R_ci_max=1.0; %Maximum inside radius of bundle is 1.0 m
R_ci=R_ci_min+(R_ci_max-R_ci_min)*rand(1,1); %Randomly chooses an inner radius between min and max values
%% Choose number of bundles
bundles_min=20; %Minimum number of sub-bundles
bunldes_max=50; %Maximum number of sub-bundles
bundles=randi([bundles_min,bunldes_max],1,1); %Randomly chooses the number of sub-bundles
%% Randomly choose other inlet parameters within constraints
%Tube Bundle dimension Constraints
r_vessel=2; %Maximum radius of inner surface of pressure vessel (2 m)
width_limit=sqrt(r_vessel^2-R_ci^2)-spacers*spacer_width; %Maximum width of bundle, excluding tie rod spacers, based on center of bundle having same area as outer annulus between bundle and pressure vessel.   
width=2*width_limit; %Start loop with dummy variable
%bank_depth=entry*(tube_layer*2)*loops*SL*D_out+spacers*spacer_width;
height_limit_total=12; %Maximum height of tube bundle is 12 m
height_limit_sub=height_limit_total/bundles-disk_thick; %Calculates maximum height of each sub-bundle based on number of sub-bundles
% H=D_out*ST*((layer_num+1)/2); %Height of sub-bundle, excluding spacer disk [m]
%Minimum and maximum values for each parameter
SL_min=1.25; %In order for SL to be less than 1.25, ST>1.5.  In order to have the SD=ST, SL=ST*cosd(30). 
              %At SL=1.25, ST=1.4434.  So, it is impossible for SL<1.25.
SL_tot_max=2.0; %Absolute max value of SL
entry_min=2; %Minimum # of manifold inlets
entry_max=4; %Maximum # of manifold inlets
tube_layer_min=2; %Minimum # of tubes per layer per manifold inlet
tube_layer_tot_max=7; % Absolute maximum allowable number of tubes in each layer
layer_num_min=20;
layer_num_max=60; % Absolute maximum allowable number of tubes layers per sub-bundle
loops_min=2;
loops_max_tot=10;
entry=randi([entry_min,entry_max]); %Choose number of liquid manifolds
while width > width_limit 
    loops_max_1=floor(width_limit/(entry*2*tube_layer_min*SL_min*D_out));
    loops_max=min(loops_max_1,loops_max_tot);
    loops=randi([loops_min,loops_max]); %Choose number of times tubes loop around bank
    tube_layer_max_1=floor(width_limit/(SL_min*entry*D_out*loops*2)); %Maximum # of tubes with given parameters
    tube_layer_max=min(tube_layer_max_1,tube_layer_tot_max);
    if tube_layer_max>tube_layer_min %If tube_row_max<tube_row_min, chooses new # of inlet manifolds
        tube_layer=randi([tube_layer_min,tube_layer_max]);
        layer_num=randi([layer_num_min,layer_num_max]);
        SL_max_1=width_limit/(entry*D_out*tube_layer*loops*2); %Restrict SL on width of bundle
        SL_max_2=2*height_limit_sub*cosd(30)/((layer_num_max+1)*D_out); %Restrict SL on height of bundle with a max # of tube layers
        SL_max=min([SL_max_1,SL_max_2,SL_tot_max]); %Max SL possible given parameters
        if SL_max>SL_min %If SL_max>SL_min, chooses new parameters with same tube diameter & thickness.
            SL=(SL_max-SL_min)*rand+SL_min;
            ST=SL/cosd(30); %In order to have the SD=ST, SL=ST*cosd(30). Forms equilateral triangle.
            width=entry*SL*D_out*(tube_layer*2)*loops;
        else
            SL=SL_min;
            ST=SL/cosd(30); %In order to have the SD=ST, SL=ST*cosd(30). Forms equilateral triangle.
            width=entry*SL*D_out*(tube_layer*2)*loops;
            height_sub=D_out*ST*((layer_num+1)/2);
            while width>width_limit || height_sub>height_limit_sub
                if width>width_limit
                    if loops>loops_min
                        loops=loops-1;
                    elseif entry>entry_min
                        entry=entry-1;
                    elseif tube_layer>tube_layer_min
                        tube_layer=tube_layer-1;
                    elseif 0.95*R_ci>=R_ci_min
                        R_ci=0.95*R_ci;
                        width_limit=sqrt(r_vessel^2-R_ci^2)-spacers*spacer_width;
                    end
                    width=entry*SL*D_out*(tube_layer*2)*loops;
                elseif height_sub>height_limit_sub
                    if layer_num>layer_num_min
                       layer_num=layer_num-1;
                    elseif D_out_in>Diameters(1)
                        Diam_pos=find(Diameters==D_out_in);
                        D_out_in=Diameters(Diam_pos-1);
                        D_out=D_out_in/.0254;
                    elseif bundles>bundles_min
                        bundles=bundles-1;
                        height_limit_sub=height_limit_total/bundles-disk_thick;
                    end
                    height_sub=D_out*ST*((layer_num+1)/2);
                end
            end
        end
    else
        if entry>entry_min
            entry=entry-1;
        elseif 0.95*R_ci>R_ci_min
            R_ci=0.95*R_ci;
        end
    end
end
input_name=sprintf('Optimization Program/Optimization_Files/Inputs/Input%d.mat',run);
save(input_name,'m_g','m_l','P_g_in','P_l_in','T_g_in','T_l_in','T_l_out','T_g_out','D_out','t','SL','ST','entry','layer_num','tube_layer','bundles','gas','heat_rod','liquid','loops','R_ci','tube_material','tube_slope','spacer_width','spacers','results_location');
range_input=sprintf('A%d:K%d',run+1,run+1);
A=[run,D_out_in,t_in,SL,ST,entry,layer_num,tube_layer,R_ci,loops,bundles];
xlswrite(results_location,A,range_input);
CTGH_0D(THEEM_model,run)
end
