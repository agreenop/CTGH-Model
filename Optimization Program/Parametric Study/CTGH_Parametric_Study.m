%This program performs a parametric study analyzing the effects of changing
%the input variables for the CTGH.
function CTGH_Parametric_Study
copyfile('Optimization Program/Parametric Study/THEEM_Input_Parametric.mat','Optimization Program/Parametric Study/THEEM_Input_temp_0D.mat')
tube_holders=[];
save('THEEM_Input_temp_0D.mat','tube_holders','-append')
copyfile('Optimization Program/Parametric Study/THEEM_Input_Parametric.mat','Optimization Program/Parametric Study/THEEM_Input_temp_2D.mat')
T_l_out=[];
T_g_out=[];
save('THEEM_Input_temp_2D.mat','T_l_out','T_g_out','-append')
clear;
load('THEEM_Input_Parametric.mat');
input_data=whos; %Lists all local workspace variables in 1 structure
data_names=strtrim(string(char(input_data.name))); %List all names of variables as vector of strings (Converts to characters, then to strings, then eliminates whitespace)
[input_variable,input_min_value,input_max_value,input_step_value]=Choose_Input_Variable; %User specifies what variable to change and over what range
clear(input_variable,'THEEM_model') %Deletes THEEM_model variable that is not relevant to calculations as well as dynamic variable's original value
data_names=data_names(data_names~=input_variable&data_names~='THEEM_model'); %Delete THEEM_model and dynamic variable from list of static variables
file_location=sprintf('Optimization Program/Parametric Study/Results_%s.xlsx',input_variable); %Name of file and relative file path
delete(file_location) %Delete previous results
data_values=cell(size(data_names));
for data_index=1:length(data_names)
data_values{data_index}=feval(@()evalin('caller',data_names(data_index)));
end
Static_names_range=sprintf('A2:A%d',length(data_values));
Static_values_range=sprintf('B2:B%d',length(data_values));
Static_headers={'Variable','Value'};
xlswrite(file_location,Static_headers,'A1:B1')
xlswrite(file_location,data_names,Static_names_range)
xlswrite(file_location,data_values,Static_values_range)
Results_headers={input_variable,'Gas Pressure Drop [bar]','Liquid Pressure  [bar]','Effectiveness','Bundle Outer Diameter [m]','Bundle Height [m]','Liquid Outlet Temperature [deg C]','Gas Outlet Temperature [deg C]','Total Heat Transfer [W]'};
xlswrite(file_location,Results_headers,2,'A1:F1')
THEEM_model='Parametric Study';
step_counter=1;
for x=input_min_value:input_step_value:input_max_value
    feval(@()assignin('caller',input_variable,x));
    save('THEEM_Input_temp_0D.mat',input_variable,'-append')
    save('THEEM_Input_temp_2D.mat',input_variable,'-append')
    CTGH_0D(THEEM_model);
    CTGH_2D(THEEM_model);
    output_data_0D=load('THEEM_Output_temp_0D.mat');
    output_data_2D=load('THEEM_Output_temp_2D.mat');
    Relevant_output=[output_data_2D.deltaP_g,output_data_2D.deltaP_l,output_data_2D.epsilon,output_data_0D.D_curve_outer,output_data_0D.H_bank,output_data_2D.T_l_outlet,output_data_2D.T_g_outlet,output_data_2D.Q_tot];
%     Relevant_output=[output_data_0D.deltaP_g,output_data_0D.deltaP_l,output_data_0D.epsilon,output_data_0D.D_curve_outer,output_data_0D.H_bank,output_data_0D.T_l_out,output_data_0D.T_g_out,output_data_0D.Q_tot];
    Dynamic_variable_range=sprintf('A%d',step_counter+1);
    Output_data_range=sprintf('B%d:F%d',step_counter+1,step_counter+1);
    xlswrite(file_location,x,2,Dynamic_variable_range)
    xlswrite(file_location,Relevant_output,2,Output_data_range)
    step_counter=step_counter+1;
end
Input_range=sprintf('A2:A%d',step_counter);
Output_range=sprintf('B2:F%d',step_counter);
Input_vector=xlsread(file_location,2,Input_range);
Output_matrix=xlsread(file_location,2,Output_range);
Results_titles=regexprep(Results_headers,' [.*','');
for plot_number=1:size(Output_matrix,2)
    figure(plot_number)
    plot(Input_vector,Output_matrix(:,plot_number))
    title(sprintf('%s vs %s',Results_titles{1}, Results_titles{plot_number+1}));
    xlabel(Results_headers(1));
    ylabel(Results_headers(plot_number+1));
end
%% Change Excel File Worksheet Names & Add Figures
full_path=which(file_location); %Generates full file path for Activex application
e = actxserver('Excel.Application'); % open Activex server
ewb = e.Workbooks.Open(full_path); % open file 
ewb.Worksheets.Item(1).Name = 'Static Variables'; % rename 1st sheet
ewb.Worksheets.Item(2).Name = 'Dynamic Variable Results'; % rename 2nd sheet
ewb.Worksheets.Add([],ewb.Worksheets.Item(2),size(Output_matrix,2));
for i=1:size(Output_matrix,2)
    fig_name=sprintf('Figure %d',i);
    img=sprintf('figure%d.png',i);
    fig_num=sprintf('-f%d',i);
    ewb.Worksheets.Item(i+2).Name = fig_name; % rename figure worksheets
    print(fig_num, img, '-dpng');
    ewb.Worksheets.Item(i+2).invoke('Pictures').Insert([pwd '\' img]);
    delete([pwd '\' img])
end
ewb.Save % save to the same file
ewb.Close(false)
e.Quit
delete('Optimization Program/Parametric Study/THEEM_Input_temp_0D.mat')
delete('Optimization Program/Parametric Study/THEEM_Input_temp_2D.mat')
delete('Optimization Program/Parametric Study/THEEM_Output_temp_0D.mat')
delete('Optimization Program/Parametric Study/THEEM_Output_temp_2D.mat')
