%This program performs a parametric study analyzing the effects of changing
%the input variables for the CTGH.
function CTGH_Parametric_Study
load('THEEM_Input_Parametric.mat');
[input_variable,input_min_value,input_max_value,input_step_value]=Choose_Input_Variable; %User specifies what variable to change and over what range
file_location='Optimization Program/Parametric Study/Results.xlsx'; %Generates file and relative file path
delete(file_location) %Delete previous results
clear(input_variable,'THEEM_model') %Deletes THEEM_model variable that is not relevant to calculations
input_data=whos; %Lists all local workspace variables in 1 structure
data_names=strtrim(string(char(input_data.name))); %List all names of variables as vector of strings (Converts to characters, then to strings, then eliminates whitespace)
data_values=cell(size(data_names));
for data_index=1:length(data_names)
data_values{data_index}=feval(@()evalin('caller',data_names(data_index)));
end
Static_names_range=sprintf('A2:A%d',length(data_values));
Static_values_range=sprintf('B2:B%d',length(data_values));
xlswrite(file_location,data_names,Static_names_range)
xlswrite(file_location,data_values,Static_values_range)
Results_headers={input_variable,'Gas Pressure Drop','Liquid Pressure Drop','Effectiveness','Bundle Diameter','Outlet Gas Temperature','Outlet Liquid Temperature'};
xlswrite(file_location,Results_headers,2,'A1:G1')
%% Change Excel File Worksheet Names
full_path=which(file_location);
e = actxserver('Excel.Application'); % open Activex server
ewb = e.Workbooks.Open(full_path); % open file 
ewb.Worksheets.Item(1).Name = 'Static Variables'; % rename 1st sheet
ewb.Worksheets.Item(2).Name = 'Dynamic Variable Results'; % rename 2nd sheet
ewb.Save % save to the same file
ewb.Close(false)
e.Quit
THEEM_model='Parametric Study';
step_counter=1;
for x=input_min_value:input_step_value:input_max_value
    feval(@()assignin('caller',input_variable,x));
            save('Optimization Program/Parametric Study/THEEM_Input_temp.mat','gas','liquid','T_g_in','T_l_in',...
            'T_g_out','T_l_out','P_g_in','P_l_in','m_g','m_l','tube_material',...
            'D_out','t','SL','ST','entry','THEEM_model','loops','tube_layer',...
            'layer_num','bundles','spacers','spacer_width','R_ci','tube_slope','heat_rod');
    Dynamic_variable_range=sprintf('A%d',step_counter+1);
    xlswrite(file_location,x,2,Dynamic_variable_range)
    step_counter=step_counter+1;
end
delete('Optimization Program/Parametric Study/THEEM_Input_temp.mat')
