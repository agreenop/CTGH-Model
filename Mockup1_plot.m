%This function will map the temperatures, pressures, and heat transfer
%outputs for the liquid coolant and gas medium for the CTGH.
function Mockup1_plot(T_l,T_g,Q,P_l,P_g,UA_matrix,Re_g_matrix,h_g_matrix,Re_l_matrix,gas,liquid)
%Plot liquid coolant temperatures
rho=10:size(T_l,1)+9;
theta = linspace(90,-270,size(T_l,2))*pi/180;
[th, r] = meshgrid(theta, rho);
figure(1)
surf(r.*cos(th)/10,r.*sin(th)/10,T_l,'linestyle','none'); 
view(2);
axis equal tight;
title([liquid ' Temperature Distribution'])
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
title(colorbar,'Degrees Celsius');
%Plot gas temperatures
rho_g=10:size(T_g,1)+9;
theta_g = linspace(90,-270,size(T_g,2))*pi/180;
[th_g, r_g] = meshgrid(theta_g, rho_g);
figure(2)
surf(r_g.*cos(th_g)/10,r_g.*sin(th_g)/10,T_g,'linestyle', 'none'); 
view(2);
axis equal tight;
title([gas ' Temperature Distribution'])
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
title(colorbar,'Degrees Celsius');
%Plot heat transfer
figure(3)
surf(r.*cos(th)/10,r.*sin(th)/10,Q,'linestyle', 'none'); 
view(2);
axis equal tight;
title('Heat Transfer Distribution')
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
title(colorbar,'Watts');
%Plot coolant pressure drop
figure(4)
surf(r.*cos(th)/10,r.*sin(th)/10,P_l,'linestyle', 'none'); 
view(2);
axis equal tight;
title([liquid ' Pressure Distribution'])
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
title(colorbar,'Units of Bar');
%Plot gas pressure drop
figure(5)
surf(r_g.*cos(th_g)/10,r_g.*sin(th_g)/10,P_g,'linestyle', 'none'); 
view(2);
axis equal tight;
title([gas ' Pressure Distribution'])
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
title(colorbar,'Units of Bar');
%Plot UA of each cell
figure(6)
surf(r.*cos(th)/10,r.*sin(th)/10,UA_matrix,'linestyle', 'none'); 
view(2);
axis equal tight;
title('UA Distribution')
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
title(colorbar,'Watts/Kelvin');
figure (7)
surf(r.*cos(th)/10,r.*sin(th)/10,Re_g_matrix,'linestyle', 'none'); 
view(2);
axis equal tight;
title([gas ' Reynolds Number Distribution'])
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
figure (8)
surf(r.*cos(th)/10,r.*sin(th)/10,h_g_matrix,'linestyle', 'none'); 
view(2);
axis equal tight;
title([gas ' Heat Transfer Coefficient (h) Distribution'])
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;
figure (9)
surf(r.*cos(th)/10,r.*sin(th)/10,Re_l_matrix,'linestyle', 'none'); 
view(2);
axis equal tight;
title([liquid ' Reynolds Number Distribution'])
xlabel('Position on X-axis, m')
ylabel('Position on Y-axis, m')
colorbar;