%This program performs a parametric study analyzing the effects of changing
%the input variables for the CTGH.
function CTGH_Parametric_Study
[input_variable,input_min_value,input_max_value,input_step_value]=Choose_Input_Variable; %User specifies what variable to change and over what range
load('THEEM_Input_Parametric.mat');
for x=input_min_value:input_step_value:input_max_value
    feval(@()assignin('caller',input_variable,x));
            save('Optimization Program/Parametric Study/THEEM_Input_temp.mat','gas','liquid','T_g_in','T_l_in',...
            'T_g_out','T_l_out','P_g_in','P_l_in','m_g','m_l','tube_material',...
            'D_out','t','SL','ST','entry','THEEM_model','loops','tube_layer',...
            'layer_num','bundles','spacers','spacer_width','R_ci','tube_slope','heat_rod');
end
delete('Optimization Program/Parametric Study/THEEM_Input_temp.mat')
