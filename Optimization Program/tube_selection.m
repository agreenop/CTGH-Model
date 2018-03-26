function [D_out,t,tube_index]=tube_selection(tube_index)
    switch tube_index
        case 1
            D_out_in=3/16;
            t_in=0.049;            
        case 2
            D_out_in=1/4;
            t_in=0.065;
        case 3
            D_out_in=5/16;
            t_in=0.065;
        case 4
            D_out_in=3/8;
            t_in=0.083;
        case 5
            D_out_in=1/2;
            t_in=0.12;
        otherwise %This makes sure that a number higher than expected was selected
            tube_index=6;
            D_out_in=5/8;
            t_in=0.12;
    end
D_out=D_out_in*0.0254; %Converts diameters to meters
t=t_in*0.0254; %Converts thickeness to meters
% Tube outside diameters and thicknesses obtained from McMaster-Carr website
% Diameters=[3/16:1/32:5/16,3/8,7/16,1/2]; %,5/8,11/16];%,3/4,7/8,1,1.25,1.375,1.5,1.75,2,2.125,2.25,2.375,2.5,3]; %Diameters in inches
% D_out_in=Diameters(randi(numel(Diameters))); %Randomly selects a tube outisde diameter for the design
% switch D_out_in %Assign tube thicknesses based on diameter in inches
%     case 1/16
%         thickness=[0.02,0.028];
%     case 3/32
%         thickness=0.02;
%     case 1/8
%         thickness=[0.02,0.028,0.035,0.049];
%     case 5/32
%         thickness=[0.035,0.049];
%     case 3/16
%         thickness=[0.02,0.028,0.035,0.049];
%     case 7/32
%         thickness=0.035;
%     case 1/4
%         thickness=[0.01,0.016,0.02,0.028,0.035,0.049,0.065,0.083];
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
% end