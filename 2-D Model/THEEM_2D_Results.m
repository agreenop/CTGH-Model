%This script will load an existing 2-D output file and give relevant outlet
%conditions and plot relevant figures.  This script essentially presents
%and plots the same information as the CTGH_2D function.  This just skips
%all of the calculations and uses an existing output file.
close all;clear;clc;
[FileName,PathName] = uigetfile('*.mat','Select 2D THEEM Output File'); %User selects existing output file
if FileName~=0 %A .mat file needs to be selected
    load([PathName,FileName]);
    if exist('T_l','var')==1 && exist('P_g','var')==1 && exist('h_l_matrix','var')==1 && exist('T_s_out_matrix','var')==1 %Checks to see if certain matrices are present to check if 2D output file
        Q_total=section*bundles*Q_actual; %Assuming all 2-D cross sections are equal in the bundle, this calculates the overall heat transfer in the bundle.
        fprintf('The effectiveness of this heat exchanger Is %4.4f.\n',epsilon)
        fprintf('The %s outlet temperature is %1.1f%cC.\n',liquid,T_l_outlet,char(176))
        fprintf('The %s outlet temperature is %1.1f%cC.\n',gas,T_g_outlet,char(176))
        fprintf('The %s pressure drop is %1.2f bar.\n',liquid,P_l_in-P_l_outlet)
        fprintf('The %s pressure drop is %1.4f bar.\n',gas,P_g_in-P_g_outlet)
        fprintf('The estimated overall heat transfer is %1.3e W.\n',Q_total)
        CTGH_plot(T_l,T_g,Q,P_l,P_g,UA_matrix,Re_g_matrix,h_g_matrix,h_l_matrix,Re_l_matrix,U_matrix,T_s_in_matrix,T_s_out_matrix,gas,liquid,R_ci,vol_cells_gap,vol_wid,spacer_width) %Plots the values
    else
       box=errordlg({'This is not a 2D output file.','Please select another file.'},'Output File Error');
       uiwait(box);
       THEEM_2D_Results
    end
else
    clear; %Ends program and clears workspace if no file is selected
end