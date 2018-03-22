% This program performs a sensitivity analysis on the CTGH 0-D code.  It
% changes the geometry of the CTGH in order to minimize the necessary
% effectiveness for the CTGH along with the salt and liquid pressure drops.
%% Initialize the Excel Spreadsheet and assign static inputs
clc;clear;
results_location=sprintf('Optimization Program/Optimization_Files/Optimization_Results.xlsx');%Store file name as variable
results_location_absolute=[pwd '\' results_location];
delete(results_location); %Delete previous results
%These inputs will not be changed by the optimization code:
m_g=1360.5;
m_l=1267;
P_g_in=199.5;
P_l_in=3;
T_g_in=367;
T_l_in=528;
T_l_out=373; %Liquid Outlet temperature [degC] 
T_g_out=516.6; %Gas outlet temperature [degC]
gas='Supercritical CO2';
liquid='Sodium';
heat_rod=1/4;
tube_material='316 Stainless Steel';
tube_slope=0.0030;
spacers=2;
spacer_width=0.0380;
disk_thick=.003;
%% Physical Constraints
D_bund_out_max=5; % Max diameter of vessel [meters]
H_bund_max=10; % Max height of vessel [meters]
bund_width_min=0.10; %Minimum width of bundle as a percentage of vessel outer diameter
sub_min_tot=20; % Absolute minimum number of sub-bundles allowed.  This may increase later, but not decrease.
sub_max_tot=60; % Absolute maximum number of sub-bundles allowed.  This may decrease later, but not increase.
layer_min_tot=30;% Absolute minimum number of tube layers per sub-bundles allowed.  This may increase later, but not decrease.
layer_max_tot=100;% Absolute maximum number of tube layers per sub-bundles allowed.  This may decrease later, but not increase.
entry_min_tot=2; % Absolute minimum number of manifolds allowed.  This may increase later, but not decrease.
entry_max_tot=7; % Absolute maximum number of manifolds allowed.  This may decrease later, but not increase.
loops_min_tot=2; % Absolute minimum number of loops allowed.  This may increase later, but not decrease.
loops_max_tot=7; % Absolute maximum number of loops allowed.  This may decrease later, but not increase.
[~,~,~,S_allow] = SS_316_prop(T_g_out); % Max allowable stress for pressure vessel [bar]
t_vessel=P_g_in*(D_bund_out_max/2)/(S_allow-0.6*P_g_in); % [meters] Thickness of vessel outer wall
save('Optimization Program/Optimization_Input_Temp.mat');
%% Create New Excel File
T={'Run','Outside Diameter [in]','Thickness [in]','Longitudinal Pitch',...
   'Transverse Pitch','Number of Inlets','Tube Layers per sub-bundle',...
   'Number of Tubes per Layer per Manifold','Inner Radius of Bundle [m]',...
   'Number of Loops', 'Number of Sub-Bundles','Total Tubes',...
   'Bundle Diameter [m]','Bundle Height [m]',...
   'Total Surface Area [m^2]','Min F Factor',...
   'Air Pressure Drop [bar]','Liquid Pressure Drop [bar]',...
   'Tube Bank Depth [m]'}; %Column labels in the Excel Spreadsheet
e = actxserver('Excel.Application'); %Use Activex controls to save time. Xlswrite opens and closes in every iteration. Activex keeps it open for entire script.
eWorkbook = e.Workbooks.Add;
eSheets1 = e.ActiveWorkbook.Sheets.get('Item',1);
eSheets1.Activate;
eActivesheetRange = e.Activesheet.get('Range','A1:S1');
eActivesheetRange.Value = T;
eWorkbook.SaveAs(results_location_absolute);
%% Start Runs
run_total=10000; %Total number of runs code will perform
h=waitbar(0,sprintf('Progress = %2.2f%%',0));
clearvars -except run_total h e eWorkbook eSheets1
for run=1:run_total
    load('Optimization Program/Optimization_Input_Temp.mat');
    h=waitbar((run-1)/run_total,h,sprintf('Progress = %2.2f%%',(run-1)/run_total*100));
%Tube Selection
if P_g_in<125 %Based on ASME BPVC Division 1
    D_out_in=1/4;
    t_in=0.035;
else %Based on S-CO2 calculations
    tube_index=randi([1,5]); %Randomly select one of the 5 tubes
    [D_out,t,tube_index]=tube_selection(tube_index);
    D_out_in=D_out/0.0254;
    t_in=t/0.0254;
end
% Pitch-to-Diameter Ratios
SL=1.256; % SL & ST hav little effect on pressure drops and effectiveness, so just assigned.  Can be reduced if bundle too large
ST=SL/cosd(30); %In order to have the SD=ST, SL=ST*cosd(30). Forms equilateral triangle.
% Tubes per layer
tube_layer=randi([3,7]); %Has little effect, so just chosen randomly between 3-6. Can be reduced if bundle too large
% Choose inner radius of tube bundle
R_ci_min=0.25; %Minimum inside radius of bundle is 0.25 m 
R_ci_max=0.5*(sqrt(-(D_bund_out_max/2)^2*bund_width_min^2+2*(D_bund_out_max/2)^2-4*(D_bund_out_max/2)*t_vessel+2*t_vessel^2)-(D_bund_out_max/2)*bund_width_min);% Maximum inside radius of bundle so that bundle makes up the minimum percentage of the vessel
R_ci=R_ci_min+(R_ci_max-R_ci_min)*rand(1,1); %Randomly chooses an inner radius between min and max values
%% Choose higher priority parameters within constraints
%Tube Bundle dimension Constraints
width_limit_max=sqrt((D_bund_out_max/2-t_vessel)^2-R_ci^2)-R_ci-spacers*spacer_width; %Maximum width of bundle, excluding tie rod spacers, based on center of bundle having same area as outer annulus between bundle and pressure vessel.   
width_limit_min=bund_width_min*(D_bund_out_max/2)-spacers*spacer_width; %The bundle should take up at least a specified percentage of the vessel
% bank_depth=entry*(tube_layer*2)*loops*SL*D_out+spacers*spacer_width;
% H=D_out*ST*((layer_num+1)/2); %Height of sub-bundle, excluding spacer disk [m]
%Minimum and maximum values for each parameter
entry_max_geom_exact=sqrt(width_limit_max/(2*tube_layer*D_out*SL));
entry_max_geom_int=floor(entry_max_geom_exact);
bank_width1=(entry_max_geom_int)*(entry_max_geom_int+1)*(2*tube_layer*D_out*SL);
if bank_width1<width_limit_max
    entry_max_geom=entry_max_geom_int+1;
else
    entry_max_geom=entry_max_geom_int;
end
loops_max_geom=entry_max_geom_int;
entry_min_geom=sqrt(width_limit_min/(2*tube_layer*D_out*SL));
loops_min_geom=entry_min_geom;
entry_min=max(ceil(entry_min_geom),entry_min_tot);
entry_max=min(floor(entry_max_geom),entry_max_tot);
loops_min=max(ceil(loops_min_geom),loops_min_tot);
loops_max=min(floor(loops_max_geom),loops_max_tot);
counter_width=1;
while loops_max<loops_min || entry_max<entry_min %If the geometry is impossible, the inner radius will be shrunk until geometry is possible within width constraint
    if mod(counter_width,2)==1 && tube_index>1
        tube_index=tube_index-1; %Select the next smallest diameter tube
        [D_out,t,tube_index]=tube_selection(tube_index);
    elseif R_ci>R_ci_min
        R_ci=max(0.9*R_ci,R_ci_min);
    elseif R_ci==R_ci_min && tube_index==1
        Quit(e)
        delete(e)
        delete(h)
        error('Geometric constraints are not physically possible.  Please adjust the constraints and try again.')
    end
    width_limit_max=D_bund_out_max/2-t_vessel-2*R_ci-spacers*spacer_width; %Calculates new width limit
    entry_max_geom_exact=sqrt(width_limit_max/(2*tube_layer*D_out*SL));
    entry_max_geom_int=floor(entry_max_geom_exact);
    bank_width1=(entry_max_geom_int)*(entry_max_geom_int+1)*(2*tube_layer*D_out*SL);
    if bank_width1<width_limit_max
        entry_max_geom=entry_max_geom_int+1;
    else
        entry_max_geom=entry_max_geom_int;
    end
    loops_max_geom=entry_max_geom_int;
    entry_min_geom=sqrt(width_limit_min/(2*tube_layer*D_out*SL));
    loops_min_geom=entry_min_geom;
    entry_min=max(ceil(entry_min_geom),entry_min_tot);
    entry_max=min(floor(entry_max_geom),entry_max_tot);
    loops_min=max(ceil(loops_min_geom),loops_min_tot);
    loops_max=min(floor(loops_max_geom),loops_max_tot);
    counter_width=counter_width+1;
end %Repeats process until geometry works
entry=randi([entry_min,entry_max]);
loops=randi([loops_min,loops_max]);
bundles_max_geom=sqrt(2*H_bund_max/(D_out*ST)); %Equations obtained from optimization problem. 
layer_max_geom=bundles_max_geom-1-disk_thick/(D_out*ST); %They maximize these variables while keeping them within the height constraint.
bundles_min=sub_min_tot;
bundles_max=min(floor(bundles_max_geom),sub_max_tot);
layer_min=layer_min_tot;
layer_max=min(floor(layer_max_geom),layer_max_tot);
while bundles_max<bundles_min || layer_max<layer_min %If the geometry is impossible, the tube diameter will be shrunk until geometry is possible within width constraint
        if tube_index>1
        tube_index=tube_index-1; %Select the next smallest diameter tube
        [D_out,t,tube_index]=tube_selection(tube_index);
        else
            Quit(e) %Close excel before throwing error
            delete(e)
            delete(h) %Close waitbar before throwing error
            error('Geometric constraints are not physically possible.  Please adjust the constraints and try again.');
        end
end
bundles=randi([bundles_min,bundles_max]); %Randomly chooses the number of sub-bundles
layer_num=randi([layer_min,layer_max]); %Randomly chooses the number of tube layers per sub-bundle
input_name=sprintf('Optimization Program/Optimization_Files/Inputs/Input%d.mat',run);
THEEM_model='0D';
save(input_name,'m_g','m_l','P_g_in','P_l_in','T_g_in','T_l_in','T_l_out','T_g_out','D_out','t','SL','ST','entry','layer_num','tube_layer','bundles','gas','heat_rod','liquid','loops','R_ci','tube_material','tube_slope','spacer_width','spacers','THEEM_model');
range_input=sprintf('A%d:K%d',run+1,run+1);
A=[run,D_out/0.0254,t/0.0254,SL,ST,entry,layer_num,tube_layer,R_ci,loops,bundles];
eActivesheetRange = e.Activesheet.get('Range',range_input);
eActivesheetRange.Value = A;
THEEM_model='Optimization';
CTGH_0D(THEEM_model,run,e)
if mod(run,100)==0
    eWorkbook.Save; %Save excel file every 100 runs
end
clearvars -except run_total h e eWorkbook eSheets1
end
eWorkbook.Save;
Quit(e)
delete(e)
h=waitbar(100,h,sprintf('Progress = %2.2f%%',100));
delete(h)
delete('Optimization Program/Optimization_Input_Temp.mat')