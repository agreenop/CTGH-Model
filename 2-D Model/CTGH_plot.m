%This function will map the temperatures, pressures, and heat transfer
%outputs for the liquid coolant and gas medium for the CTGH.
function CTGH_plot(T_l,T_g,Q,P_l,P_g,UA_matrix,Re_g_matrix,h_g_matrix,h_l_matrix,Re_l_matrix,U_matrix,gas,liquid,R_ci,vol_cells_gap,vol_wid,spacer_width)
%Plot liquid coolant temperatures
rho=zeros(1,size(T_l,1));
rho(1)=R_ci+vol_wid/2;
space_count=1;
for x=2:size(T_l,1)
    space_count=space_count+1;
    if space_count==vol_cells_gap+1
        rho(x)=rho(x-1)+(spacer_width+vol_wid)/2;
    elseif space_count==vol_cells_gap+2
        rho(x)=rho(x-1)+(spacer_width+vol_wid)/2;
        space_count=1;
    else
        rho(x)=rho(x-1)+vol_wid;
    end
end
theta = linspace(90,-270,size(T_l,2))*pi/180;
[th, r] = meshgrid(theta, rho);
figure(1)
surf(r.*cos(th),r.*sin(th),T_l,'linestyle','none'); 
view(2);
axis equal tight;
title([liquid ' Temperature Distribution'])
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
title(colorbar,'Degrees Celsius');
%Plot gas temperatures
rho_g=zeros(1,size(T_g,1));
rho_g(1)=R_ci;
space_count_g=1;
for y=2:size(T_g,1)
    space_count_g=space_count_g+1;
    if space_count_g==vol_cells_gap+2
        rho_g(y)=rho_g(y-1)+spacer_width;
        space_count_g=1;
    else
        rho_g(y)=rho_g(y-1)+vol_wid;
    end
end
theta_g = linspace(90,-270,size(T_g,2))*pi/180;
[th_g, r_g] = meshgrid(theta_g, rho_g);
figure(2)
surf(r_g.*cos(th_g),r_g.*sin(th_g),T_g,'linestyle', 'none'); 
view(2);
axis equal tight;
title([gas ' Temperature Distribution'])
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
title(colorbar,'Degrees Celsius');
%Plot heat transfer
figure(3)
surf(r.*cos(th),r.*sin(th),Q,'linestyle', 'none');
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
view(2);
axis equal tight;
title('Heat Transfer Distribution')
colorbar;
title(colorbar,'Watts');
%Plot coolant pressure drop
figure(4)
surf(r.*cos(th),r.*sin(th),P_l,'linestyle', 'none'); 
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
view(2);
axis equal tight;
title([liquid ' Pressure Distribution'])
colorbar;
title(colorbar,'Units of Bar');
%Plot gas pressure drop
figure(5)
surf(r_g.*cos(th_g),r_g.*sin(th_g),P_g,'linestyle', 'none'); 
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
view(2);
axis equal tight;
title([gas ' Pressure Distribution'])
colorbar;
title(colorbar,'Units of Bar');
%Plot UA of each volume cell
figure(6)
surf(r.*cos(th),r.*sin(th),UA_matrix,'linestyle', 'none'); 
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
view(2);
axis equal tight;
title('UA Distribution')
colorbar;
title(colorbar,'Watts/Kelvin');
%Plot Gas Reynolds Number
figure (7)
surf(r.*cos(th),r.*sin(th),Re_g_matrix,'linestyle', 'none'); 
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
view(2);
axis equal tight;
title([gas ' Reynolds Number Distribution'])
colorbar;
%Plot Liquid Reynolds Number
figure (8)
surf(r.*cos(th),r.*sin(th),Re_l_matrix,'linestyle', 'none'); 
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
view(2);
axis equal tight;
title([liquid ' Reynolds Number Distribution'])
colorbar;
%Plot Gas Convective Heat Transfer Coefficient
figure (9)
surf(r.*cos(th),r.*sin(th),h_g_matrix,'linestyle', 'none'); 
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
view(2);
axis equal tight;
title([gas ' Heat Transfer Coefficient (h) Distribution'])
colorbar;
title(colorbar,'W/(m^2*K)');
%Plot Liquid Convective Heat Transfer Coefficient
figure (10)
surf(r.*cos(th),r.*sin(th),h_l_matrix,'linestyle', 'none'); 
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
view(2);
axis equal tight;
title([liquid ' Heat Transfer Coefficient (h) Distribution'])
colorbar;
title(colorbar,'W/(m^2*K)');
%Plot Overall Heat Transfer Coefficient, U
figure (11)
surf(r.*cos(th),r.*sin(th),U_matrix,'linestyle', 'none'); 
view(2);
axis equal tight;
title(' U distribution')
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
title(colorbar,'W/(m^2*K)');