clc; clear all; close all; format long; format compact;
% this code uses the Imperial units (in,lbf,in-lbf)
% To use this code, please only Run Section
% Run this first section for all data calculations after inputing correct
% variables
% Then run either sections 2-4 to plot the graphs for the specific shaft
% plots

% inputs
sut = 68*10^3; % [psi] Ultimate Tensile of 1020 CD Steel
sy = 57*10^3; % [psi] Yielding of 1020 CD Steel
E = 29*10^6; % [psi] Young's Modulus of 1020 CD Steel
Dc = 1; Df = 1; Dg = 1; Dj = 1; % Diameters of gears
GR_CF = 16/73; % gear ratio = (Number of teeth of Gear 2)/(Number of teeth of Gear 3)
GR_JG = 16/73; % gear ratio = (Number of teeth of Gear 4)/(Number of teeth of Gear 5)

% Assumptions
% Ft=Ct; Jt = Gt;
theta = pi/9; % Pressure Angle in radians estimated to be 20
theta2 = pi/9; % pressure angle of the 2nd pair
P1 = 50.25*3; % lbs
V1 = 260*3; % lbs
T1 = 38.6*3; % in-lbs

%% Shear, Moment and Torque Calculations for Input Shaft
radius_C = Dc/2; % radius of Gear C
Tc = T1;
Cz = Tc/radius_C;
Cy = Cz*tan(theta);
Ay = (8.625*V1 + 1.875*Cy)/4.125; % length
By = V1 + Cy - Ay;
Bx = P1;
Az = -1.875/4.125*Cz; % length
Bz = -2.25/4.125*Cz; % length
% calculating
x1 = [0,4, 4.5+2.25, 4.5+2.25+1.875]; % length
Force1_x = [P1,P1,P1,Bx];
Force1_y = [-V1,-V1+Ay,-V1+Ay-Cy,-V1+Ay-Cy+By];
Force1_z = [0,Az,Az+Cz,Az+Cz+Bz];
Torque1 = [T1,T1,T1-Tc,T1-Tc];
Shaft1 = linspace(0,8.625,863); % initializing variables for plotting length
Force_Y_Shaft = zeros(size(Shaft1)); Force_Z_Shaft1 = zeros(size(Shaft1)); 
Moment_Y_Shaft1 = zeros(size(Shaft1)); Moment_Z_Shaft1 = zeros(size(Shaft1));
Torque_Shaft1 = zeros(size(Shaft1)); Force_X_Shaft1 = zeros(size(Shaft1));
% inputing the values of Forces as a function of x
Force_X_Shaft1(1:(x1(4)*100)) = Force1_x(1);
Torque_Shaft1(1:(x1(3)*100)) = Torque1(1);
Force_Y_Shaft1(1:(x1(2)*100)) = Force1_y(1);
Force_Y_Shaft1((x1(2)*100+1):round(x1(3)*100)) = Force1_y(2); 
Force_Y_Shaft1((round(x1(3)*100)+1):(x1(4)*100)) = Force1_y(3);
Force_Y_Shaft1(round(x1(4)*100)) = Force1_y(4);
Force_Z_Shaft1(1:(x1(2)*100)) = Force1_z(1);
Force_Z_Shaft1((x1(2)*100+1):(x1(3)*100)) = Force1_z(2); 
Force_Z_Shaft1((x1(3)*100+1):((x1(4)*100))) = Force1_z(3);
Force_Z_Shaft1(round(x1(4)*100)) = Force1_z(4);
% integrating Shear forces for Moments
for n = 2:813
    Moment_Y_Shaft1(n) = Moment_Y_Shaft1(n-1) + (Force_Z_Shaft1(n-1)+Force_Z_Shaft1(n))/2 * 0.01;
    Moment_Z_Shaft1(n) = Moment_Z_Shaft1(n-1) + (Force_Y_Shaft1(n-1)+Force_Y_Shaft1(n))/2 * 0.01;
end

%% Shear, Moment and Torque Calculations for Intermediate Shaft
% Starting Analysis for Shaft 2 (Intermediate Shaft)
r_g = Dg/2; % radius of gear G
Fy = Cy;
Fz = Cz;
Tf = Tc*GR_CF;
Tg = Tf;
Gz = Tg/r_g;
Gy = Gz*tan(theta2);
Dy = (-2*Gy-7.75*Fy)/10; % Length
Ey = -Gy-Fy-Dy;
Dz = (-2*Gz-7.75*Fz)/10; % Length
Ez = -Gz-Fz-Dz; % testing
x2 = [0,2.25,8,10]; % Length
Force2_y = [Dy,Dy+Fy,Dy+Fy+Gy,Dy+Fy+Gy+Ey];
Force2_z = [Dz,Dz+Fz,Dz+Fz+Gz,Dz+Fz+Gz+Ez];
Torque2 = [0,Tf,Tf-Tg,Tf-Tg];
Shaft2 = linspace(0,10,1001); % initializing variables for plotting
Force_Y_Shaft2 = zeros(size(Shaft2)); Force_Z_Shaft2 = zeros(size(Shaft2)); 
Moment_Y_Shaft2 = zeros(size(Shaft2)); Moment_Z_Shaft2 = zeros(size(Shaft2));
Torque_Shaft2 = zeros(size(Shaft2));

% inputing the values of Forces as a function of x
Torque_Shaft2((x2(2)*100+1):(x2(4)*100)) = Torque2(2);
Force_Y_Shaft2(1:(x2(2)*100)) = Force2_y(1);
Force_Y_Shaft2((x2(2)*100+1):(x2(3)*100)) = Force2_y(2); 
Force_Y_Shaft2((x2(3)*100+1):((x2(4)*100))) = Force2_y(3);
Force_Y_Shaft2(x2(4)*100+1) = Force2_y(4);
Force_Z_Shaft2(1:(x2(2)*100)) = Force2_z(1);
Force_Z_Shaft2((x2(2)*100+1):(x2(3)*100)) = Force2_z(2); 
Force_Z_Shaft2((x2(3)*100+1):((x2(4)*100))) = Force2_z(3);
Force_Z_Shaft2(x2(4)*100+1) = Force2_z(4);
% integrating Shear forces for Moments
for n = 2:1001
    Moment_Y_Shaft2(n) = Moment_Y_Shaft2(n-1) + (Force_Z_Shaft2(n-1)+Force_Z_Shaft2(n))/2 * 0.01;
    Moment_Z_Shaft2(n) = Moment_Z_Shaft2(n-1) + (Force_Y_Shaft2(n-1)+Force_Y_Shaft2(n))/2 * 0.01;
end

%% Shear, Moment and Torque Calculations for Output Shaft
% Starting Analysis for Shaft 3 (Output Shaft)
Tj = GR_JG*Tg;
Jy = Gy;
Jz = Gz;
Tk = Tj;
Iy = 1.625/3.625*Jy;
Hy = Jy - Iy;
Iz = -1.625/3.625*Jz; % Length
Hz = -Iz - Jz;
x3 = [0,1.625,3.625,5.625]; % Length
Force3_y = [Hy,Hy-Jy,Hy-Jy+Iy,Hy-Jy+Iy];
Force3_z = [Hz,Hz+Jz,Hz+Jz+Iz,Hz+Jz+Iz];
Torque3 = [0,Tj,Tj,Tj-Tk];
Shaft3 = linspace(0,5.625,563); % initializing variables for plotting length
Force_Y_Shaft3 = zeros(size(Shaft3)); Force_Z_Shaft3 = zeros(size(Shaft3)); 
Moment_Y_Shaft3 = zeros(size(Shaft3)); Moment_Z_Shaft3 = zeros(size(Shaft3));
Torque_Shaft3 = zeros(size(Shaft3));
% inputing the values of Forces as a function of x
Torque_Shaft3((x3(2)*100+1):(x3(4)*100)) = Torque3(2);
Torque_Shaft3(round(x3(4)*100)) = Torque3(4);
Force_Y_Shaft3(1:round(x3(2)*100)) = Force3_y(1);
Force_Y_Shaft3((round(x3(2)*100)+1):round(x3(3)*100)) = Force3_y(2); 
Force_Y_Shaft3((round(x3(3)*100)+1):(round(x3(4)*100))) = Force3_y(3);       
Force_Y_Shaft3(round(x3(4)*100)) = Force3_y(4);
Force_Z_Shaft3(1:round(x3(2)*100)) = Force3_z(1);
Force_Z_Shaft3((round(x3(2)*100)+1):round(x3(3)*100)) = Force3_z(2); 
Force_Z_Shaft3((round(x3(3)*100)+1):(round(x3(4)*100))) = Force3_z(3);
Force_Z_Shaft3(round(x3(4)*100)) = Force3_z(4);
% integrating Shear forces for Moments
for n = 2:563
    Moment_Y_Shaft3(n) = Moment_Y_Shaft3(n-1) + (Force_Z_Shaft3(n-1)+Force_Z_Shaft3(n))/2 * 0.01;
    Moment_Z_Shaft3(n) = Moment_Z_Shaft3(n-1) + (Force_Y_Shaft3(n-1)+Force_Y_Shaft3(n))/2 * 0.01;
end
%% Fatigue Analysis
% Fatigue Calculations for Input Shaft
% Can vary the lengths here for moment calculations
M1 = sqrt((Bz*1.875)^2+(By*1.875)^2); % Square sum of Moment for shaft 1
T1 = Tc; % Torque at Gear 
[dd1, DD1, r1, ~, ~, ~, ~, ~, ny1] = shaft(M1, T1, sut, sy);
fprintf('Shaft 1 Yielding FOS = %f\n',ny1)
fprintf('Shaft 1 Fatigue FOS = %f\n', 2)
fprintf('Shaft 1 filet radius = %f [in]\n', r1)

% Fatigue Calculations for Intermediate Shaft
% Can vary the lengths here for moment calculations
M2 = sqrt((Fz*2.25)^2+(Fy*2.25)^2); % Square sum of Moment for shaft 2
T2 = Tf; % Torque at Gear
[dd2, DD2, r2, ~, ~, ~, ~, ~, ny2] = shaft(M2, T2, sut, sy);
fprintf('Shaft 2 Yielding FOS = %f\n',ny2)
fprintf('Shaft 2 Fatigue FOS = %f\n', 2)
fprintf('Shaft 2 filet radius = %f [in]\n', r2)

% Fatigue Calculations for Output Shaft
% Can vary the lengths here for moment calculations
M3 = sqrt((Jz*1.625)^2+(Jy*1.625)^2); % Square sum of Moment for shaft 3
T3 = Tj; % Torque at Gear
[dd3, DD3, r3, ~, ~, ~, ~, ~, ny3] = shaft(M3, T3, sut, sy);
fprintf('Shaft 3 Yielding FOS = %f\n',ny3)
fprintf('Shaft 3 Fatigue FOS = %f\n', 2)
fprintf('Shaft 3 filet radius = %f [in]\n', r3)

%% Deflection calculations
Bearing_Slope_max = 0.0005; % radians
SpurGear_Deflect_max = 0.003; % inches

Moment_avg_Shaft1 = sqrt(Moment_Y_Shaft1.^2 + Moment_Z_Shaft1.^2);
Moment_avg_Shaft2 = sqrt(Moment_Y_Shaft2.^2 + Moment_Z_Shaft2.^2);
Moment_avg_Shaft3 = sqrt(Moment_Y_Shaft3.^2 + Moment_Z_Shaft3.^2);

% Can change the distances of the shaft here, denoted as x
d1 = [dd1, DD1, dd1, dd1/1.1]; x_shaft1 = [0, 2, 5, 6, 8.625];
d2 = [dd2/1.1, dd2, DD2, dd2, dd2/1.1]; x_shaft2 = [0,2,5,6,8.625];
d3 = [dd3/1.1, dd3, DD3, dd3]; x_shaft3 = [0,2,5,6,8.625];

Inertia1 = pi/16 .* d1; Inertia2 = pi/16 .* d2; Inertia3 = pi/16 .* d3;

Deflection1 = zeros(size(Shaft1)); Theta1 = zeros(size(Shaft1)); Inertia_Shaft1 = ones(size(Shaft1)) * Inertia1(1);
Deflection2 = zeros(size(Shaft2)); Theta2 = zeros(size(Shaft2)); Inertia_Shaft2 = ones(size(Shaft1)) * Inertia2(1);
Deflection3 = zeros(size(Shaft3)); Theta3 = zeros(size(Shaft3)); Inertia_Shaft3 = ones(size(Shaft1)) * Inertia3(1);

for i = 1:size(d1)
    Inertia_Shaft1(round((x_shaft1(i)*100)+1):(x_shaft1(i+1)*100)) = Inertia1(i);
end
for i = 1:size(d2)
    Inertia_Shaft2(round((x_shaft2(i)*100)+1):(x_shaft2(i+1)*100)) = Inertia2(i);
end
for i = 1:size(d3)
    Inertia_Shaft3(round((x_shaft3(i)*100)+1):(x_shaft3(i+1)*100)) = Inertia3(i);
end

for n = 2:size(Shaft1) 
    Theta1(n) = Theta1(n-1) + 1./(E*Inertia1(n))*(Moment_avg_Shaft1(n-1) + Moment_avg_Shaft1(n))/2 * 0.01;
end
for n = 2:size(Shaft2) 
    Theta2(n) = Theta2(n-1) + 1./(E*Inertia2(n))*(Moment_avg_Shaft2(n-1) + Moment_avg_Shaft2(n))/2 * 0.01;
end
for n = 2:size(Shaft3) 
    Theta3(n) = Theta3(n-1) + 1./(E*Inertia3(n))*(Moment_avg_Shaft3(n-1) + Moment_avg_Shaft3(n))/2 * 0.01;
end

for n = 2:size(Shaft1) 
    Deflection1(n) = Deflection1(n-1) + (Theta1(n-1) + Theta1(n))/2 * 0.01;
end
for n = 2:size(Shaft2) 
    Deflection2(n) = Deflection2(n-1) + (Theta2(n-1) + Theta2(n))/2 * 0.01;
end
for n = 2:size(Shaft3) 
    Deflection3(n) = Deflection3(n-1) + (Theta3(n-1) + Theta3(n))/2 * 0.01;
end

Critical_Gear_deflection = 0.003; % inches
Critical_Bearing_Slope = 0.0005; % radians

% Bearing A: right - 4 in, left - 5 in 
% Bearing B: right - 8.25 in
% Bearing D: left - 0.5 in
% Bearing E: right - 9.5 in
% Bearing H: left - 0.375
% Bearing I: right - 3.125. left - 4.125
Bear_left_A = false; Bear_right_A = false; Bear_right_B = false;
Bear_left_D = false; Bear_right_E = false; Bear_left_H = false;
Bear_right_I = false; Bear_left_I = false;

if (Theta1(4*100+1) < Critical_Bearing_Slope)
    Bear_right_A = true; end
if (Theta1(5*100+1) < Critical_Bearing_Slope)
    Bear_left_A = true; end
if (Theta1(8.25*100+1) < Critical_Bearing_Slope)
    Bear_right_B = true; end
if (Theta2(0.5*100+1) < Critical_Bearing_Slope)
    Bear_left_D = true; end
if (Theta2(9.5*100+1) < Critical_Bearing_Slope)
    Bear_right_E = true; end
if (Theta3(round(0.375*100)) < Critical_Bearing_Slope)
    Bear_left_H = true; end
if (Theta3(round(3.125*100)) < Critical_Bearing_Slope)
    Bear_right_I = true; end
if (Theta3(round(4.125*100)) < Critical_Bearing_Slope)
    Bear_left_I = true; end

% Gear C: right - 1.25, left - 3.25
% Gear F: right - 1.25, left - 3.25
% Gear G: right - 7.25, left - 8.75
% Gear J: right - 0.875, left - 2.375
Gear_right_C = false; Gear_left_C = false;
Gear_right_F = false; Gear_left_F = false;
Gear_right_G = false; Gear_left_G = false;
Gear_right_J = false; Gear_left_J = false;

if (Deflection1(1.25*100+1) < Critical_Gear_deflection)
    Gear_right_C = true; end
if (Deflection1(3.25*100+1) < Critical_Gear_deflection)
    Gear_left_C = true; end
if (Deflection2(1.25*100+1) < Critical_Gear_deflection)
    Gear_right_C = true; end
if (Deflection2(3.25*100+1) < Critical_Gear_deflection)
    Gear_left_C = true; end
if (Deflection2(7.25*100+1) < Critical_Gear_deflection)
    Gear_right_C = true; end
if (Deflection2(8.75*100+1) < Critical_Gear_deflection)
    Gear_left_C = true; end
if (Deflection3(round(0.875*100)) < Critical_Gear_deflection)
    Gear_right_C = true; end
if (Deflection3(round(2.375)) < Critical_Gear_deflection)
    Gear_left_C = true; end

disp("Shaft 1:")
fprintf("Bearing A: left: " + Bear_left_A + " right: " + Bear_right_A + "\n")
fprintf("Bearing B: right: " + Bear_right_B + "\n")
fprintf("Gear C:    left: " + Gear_left_C + " right: " + Gear_right_C + "\n")
fprintf("Rotational speed of Shaft 1 = %f [rpm]\n", 200);

disp("Shaft 2:")
fprintf("Bearing D: left: " + Bear_left_D + "\n")
fprintf("Bearing E: right: " + Bear_right_E + "\n")
fprintf("Gear F: left: " + Gear_left_F + " right: " + Gear_right_F + "\n")
fprintf("Gear G: left: " + Gear_left_G + " right: " + Gear_right_G + "\n")
fprintf("Rotational speed of Shaft 2 = %f [rpm]\n", 200/GR_CF);

disp("Shaft 3:")
fprintf("Bearing H: left: " + Bear_left_H + "\n")
fprintf("Bearing I: left: " + Bear_left_I + " right: " + Bear_right_I + "\n")
fprintf("Gear J:    left: " + Gear_left_J + " right: " + Gear_right_J + "\n")
fprintf("Rotational speed of Shaft 3 = %f [rpm]\n", 200/GR_JG/GR_JG);
%% Plots for Shaft 1
close all;
figure(1)
plot(Shaft1,Torque_Shaft1,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Torque (in-lbs)")
title("Torque from right to left")

figure(2)
plot(Shaft1,Force_X_Shaft1,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Axial Force (lbs)")
title("Axial Force from right to left")

figure(3)
plot(Shaft1,Force_Y_Shaft1,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Shear Force Y(lbs)")
title("Shear Force Y from right to left")

figure(4)
plot(Shaft1,Force_Z_Shaft1,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Shear Force Z(lbs)")
title("Shear Force Z from right to left")

figure(5)
plot(Shaft1,Moment_Y_Shaft1,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Moment Y(in-lbs)")
title("Moment Y from right to left")

figure(6)
plot(Shaft1,Moment_Z_Shaft1,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Moment Z(in-lbs)")
title("Moment Z from right to left")

%% Plots for Shaft 2
close all;
figure(1)
plot(Shaft2,Torque_Shaft2,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Torque (in-lbs)")
title("Torque from right to left")

figure(2)
plot(Shaft2,Force_Y_Shaft2,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Shear Force Y(lbs)")
title("Shear Force Y from right to left")

figure(3)
plot(Shaft2,Force_Z_Shaft2,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Shear Force Z(lbs)")
title("Shear Force Z from right to left")

figure(4)
plot(Shaft2,Moment_Y_Shaft2,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Moment Y(in-lbs)")
title("Moment Y from right to left")

figure(5)
plot(Shaft2,Moment_Z_Shaft2,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Moment Z(in-lbs)")
title("Moment Z from right to left")

%% Plots for Shaft 3
close all;
figure(1)
plot(Shaft3,Torque_Shaft3,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Torque (in-lbs)")
title("Torque from right to left")

figure(2)
plot(Shaft3,Force_Y_Shaft3,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Shear Force Y(lbs)")
title("Shear Force Y from right to left")

figure(3)
plot(Shaft3,Force_Z_Shaft3,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Shear Force Z(lbs)")
title("Shear Force Z from right to left")

figure(4)
plot(Shaft3,Moment_Y_Shaft3,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Moment Y(in-lbs)")
title("Moment Y from right to left")

figure(5)
plot(Shaft3,Moment_Z_Shaft3,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Moment Z(in-lbs)")
title("Moment Z from right to left")
