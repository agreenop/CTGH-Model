%This program performs a parametric study analyzing the effects of changing
%the input variables for the CTGH.
function CTGH_Parametric_Study
copyfile('Optimization Program/Parametric Study/THEEM_Input_Parametric.mat','Optimization Program/Parametric Study/THEEM_Input_temp_0D.mat') %Makes parametric input file the temporary 0D mat file used by the 0D calculations
tube_holders=[]; %Makes tube_holders an empty variable.  It is not used in the 0D calculations
save('THEEM_Input_temp_0D.mat','tube_holders','-append') %Adds empty variable to 0D mat file
copyfile('Optimization Program/Parametric Study/THEEM_Input_Parametric.mat','Optimization Program/Parametric Study/THEEM_Input_temp_2D.mat') %Makes parametric input file the temporary 2D mat file used by the 2D calculations
T_l_out=[]; %Outlet temperatures become empty variables. They are not used in 2D program.
T_g_out=[];
save('THEEM_Input_temp_2D.mat','T_l_out','T_g_out','-append') %Adds empty variable to 2D mat file
clear; %Clear workspace after modifying mat files
load('THEEM_Input_Parametric.mat');
clear('THEEM_model') %Deletes THEEM_model variable that is not relevant to calculations
input_data=whos; %Lists all local workspace variables, which are variables from input file, in 1 structure
data_names=strtrim(string(char(input_data.name))); %List all names of variables as vector of strings (Converts to characters, then to strings, then eliminates whitespace)
[cancel_program,input_variable,input_min_value,input_max_value,input_step_value]=Choose_Input_Variable(data_names); %User specifies what variable to change/ make dynamic, the range of its values, and its step size
if strcmp(cancel_program,'True')==1
    return %Ends the program if user does not select variable or cicks "Cancel" Button
end
clear(input_variable) %Deletes dynamic variable's original value from workspace
data_names=data_names(data_names~=input_variable&data_names~='THEEM_model'); %Delete THEEM_model and dynamic variable from list of static variables
file_location=sprintf('Optimization Program/Parametric Study/Results_%s.xlsx',input_variable); %Creates name of excel file with relative file path where results are stored
delete(file_location) %Delete previous results file
data_values=cell(size(data_names)); %Creates cell for static variables values
for data_index=1:length(data_names) 
data_values{data_index}=feval(@()evalin('caller',data_names(data_index))); %Store static variable values in a cell
end
Static_names_range=sprintf('A2:A%d',length(data_values)); %Range of cells to hold static variable names in excel
Static_values_range=sprintf('B2:B%d',length(data_values)); %Range of cells to hold static variable values in excel
Static_headers={'Variable','Value'}; %Headers for static variables in excel
xlswrite(file_location,Static_headers,'A1:B1') %Writes headers to excel
xlswrite(file_location,data_names,Static_names_range) %Writes static variable names to excel
xlswrite(file_location,data_values,Static_values_range) %Writes static variable values to excel
Results_headers={input_variable,'Gas Pressure Drop [bar]','Liquid Pressure  [bar]','Effectiveness','Bundle Outer Diameter [m]','Bundle Height [m]','Liquid Outlet Temperature [deg C]','Gas Outlet Temperature [deg C]','Total Heat Transfer [W]'}; %Headers for dynamic variable and results worksheet
xlswrite(file_location,Results_headers,2,'A1:I1') %Writes results headers to excel
THEEM_model='Parametric Study'; %Reassigns THEEM_model variable that was previously deleted
step_counter=1; %Start counter
for x=input_min_value:input_step_value:input_max_value %Loops through the range of values given by the user
    feval(@()assignin('caller',input_variable,x)); % Assigns dynamic variable's current value 
    save('THEEM_Input_temp_0D.mat',input_variable,'-append') % Updates 0D mat file with new value
    save('THEEM_Input_temp_2D.mat',input_variable,'-append') % Updates 2D mat file with new value
    CTGH_0D(THEEM_model); %Runs 0D model, which is useful for geometric calculations, with current value
    CTGH_2D(THEEM_model); %Runs 2D model, which is useful for everything else, with current value
    output_data_0D=load('THEEM_Output_temp_0D.mat');
    output_data_2D=load('THEEM_Output_temp_2D.mat');
    Relevant_output=[output_data_2D.deltaP_g,output_data_2D.deltaP_l,output_data_2D.epsilon,output_data_0D.D_curve_outer,output_data_0D.H_bank,output_data_2D.T_l_outlet,output_data_2D.T_g_outlet,output_data_2D.Q_tot]; %Saves desired output from both 0D and 2D codes to a vector
%     Relevant_output=[output_data_0D.deltaP_g,output_data_0D.deltaP_l,output_data_0D.epsilon,output_data_0D.D_curve_outer,output_data_0D.H_bank,output_data_0D.T_l_out,output_data_0D.T_g_out,output_data_0D.Q_tot];
    Dynamic_variable_range=sprintf('A%d',step_counter+1); % Excel cell where dynamic variable variable is written (Starts with A2)
    Output_data_range=sprintf('B%d:I%d',step_counter+1,step_counter+1); % Excel cells where results are written (Starts with B2)
    xlswrite(file_location,x,2,Dynamic_variable_range) %Writes dynamic variable in excel on 2nd sheet
    xlswrite(file_location,Relevant_output,2,Output_data_range) %Writes results in excel on 2nd sheet
    step_counter=step_counter+1; % Moves counter forward
end
Input_range=sprintf('A2:A%d',step_counter); %Range of cells that excel reads for dynamic variable values
Output_range=sprintf('B2:I%d',step_counter); %Range of cells that excel reads for results values
Input_vector=xlsread(file_location,2,Input_range); % Assigns dynamic variable range of values in excel to a Matlab vector
Output_matrix=xlsread(file_location,2,Output_range); % Assigns results range of values in excel to a Matlab matrix
Results_titles=regexprep(Results_headers,' [.*',''); % Removes units in brackets from results headers. Ex: [bar]
for plot_number=1:size(Output_matrix,2) %Plots each result against the dynamic variable
    figure(plot_number) 
    plot(Input_vector,Output_matrix(:,plot_number))
    title(sprintf('%s vs %s',Results_titles{1}, Results_titles{plot_number+1})); %Plot title will not have units in brackets
    xlabel(Results_headers(1));
    ylabel(Results_headers(plot_number+1)); %Y axis label will have units in brackets
end
%% Change Excel File Worksheet Names & Add Figures
full_path=which(file_location); %Generates full file path of excel workbook for Activex application instead of relative path.
e = actxserver('Excel.Application'); % Open Activex server.  Uses VBA language for most of this application.
ewb = e.Workbooks.Open(full_path); % Open excel workbook 
ewb.Worksheets.Item(1).Name = 'Static Variables'; % Rename 1st sheet
ewb.Worksheets.Item(2).Name = 'Dynamic Variable Results'; % Rename 2nd sheet
ewb.Worksheets.Add([],ewb.Worksheets.Item(2),size(Output_matrix,2)); %Adds a worksheet for each Matlab figure
for i=1:size(Output_matrix,2)
    fig_name=sprintf('Figure %d',i); %Name of figure for worksheet names
    fig_num=sprintf('-f%d',i); %Generates figure number of desired figured to be used
    ewb.Worksheets.Item(i+2).Activate %Changes active worksheet
    ewb.Worksheets.Item(i+2).Name = fig_name; % Rename figure worksheets
    print(fig_num,'-clipboard', '-dmeta') %Saves desired Matlab plot to system clipboard
%     ewb.Worksheets.Item(i+2).PasteSpecial;
    ewb.Activesheet.PasteSpecial; %Pastes desired figure into active worksheet     
    ewb.Save % Save excel workbook
end
ewb.Close(false) %Closes excel workbook
e.Quit %Exits out of Activex server
delete('Optimization Program/Parametric Study/THEEM_Input_temp_0D.mat') %Delete temporary input and output files generated by program
delete('Optimization Program/Parametric Study/THEEM_Input_temp_2D.mat')
delete('Optimization Program/Parametric Study/THEEM_Output_temp_0D.mat')
delete('Optimization Program/Parametric Study/THEEM_Output_temp_2D.mat')
