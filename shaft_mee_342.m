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
Dc = 6.08; Df = 1.33; Dg = 6.08; Dj = 1.33; % Diameters of gears
GR_CF = 16/73; % gear ratio = (Number of teeth of Gear 2)/(Number of teeth of Gear 3)
GR_JG = 16/73; % gear ratio = (Number of teeth of Gear 4)/(Number of teeth of Gear 5)

% Assumptions
% Ft=Ct; Jt = Gt;
theta = pi/9; % Pressure Angle in radians estimated to be 20
theta2 = pi/9; % pressure angle of the 2nd pair
P1 = 50.25*3; % lbs
V1 = 260*3; % lbs
T1 = 38.6*3; % in-lbs

% Shear, Moment and Torque Calculations for Input Shaft
radius_C = Dc/2; % radius of Gear C
Tc = T1;
Cz = 456.8646;
Cy = Cz * tan(theta);
Ay = (7.875*V1 + 1.5*Cy)/3.0625; % Momennt at B
By = Cy - Ay + V1; % Sum of Shear Y
Bx = P1;
Az = -1.5/3.0625*Cz; % moment at B
Bz = -1.5625/3.0625*Cz; % moment at A
% calculating
x1 = [0,4.8125, 4.8125+1.5625, 4.8125+1.5625+1.5]; % length
Force1_x = [P1,P1,P1,Bx];
Force1_y = [-V1,-V1+Ay,-V1+Ay-Cy,-V1+Ay-Cy+By];
Force1_z = [0,Az,Az+Cz,Az+Cz+Bz];
Torque1 = [T1,T1,T1-Tc,T1-Tc];
Shaft1 = linspace(0,7.875,788); % initializing variables for plotting length
Force_Y_Shaft = zeros(size(Shaft1)); Force_Z_Shaft1 = zeros(size(Shaft1)); 
Moment_Y_Shaft1 = zeros(size(Shaft1)); Moment_Z_Shaft1 = zeros(size(Shaft1));
Torque_Shaft1 = zeros(size(Shaft1)); Force_X_Shaft1 = zeros(size(Shaft1));
% inputing the values of Forces as a function of x
Force_X_Shaft1(1:(x1(4)*100)) = Force1_x(1);
Torque_Shaft1(1:(x1(3)*100)) = Torque1(1);
Force_Y_Shaft1(1:(x1(2)*100)) = Force1_y(1);
Force_Y_Shaft1((x1(2)*100+1):round(x1(3)*100+1)) = Force1_y(2); round(x1(3)*100) % ueahgewadfawkudhwakudhwakudh
Force_Y_Shaft1((round(x1(3)*100)+1):(x1(4)*100)) = Force1_y(3);
Force_Y_Shaft1(round(x1(4)*100)) = Force1_y(4);
Force_Z_Shaft1(1:(x1(2)*100)) = Force1_z(1);
Force_Z_Shaft1((x1(2)*100+1):(x1(3)*100)) = Force1_z(2); 
Force_Z_Shaft1((x1(3)*100+1):((x1(4)*100))) = Force1_z(3);
Force_Z_Shaft1(round(x1(4)*100)) = Force1_z(4);
% integrating Shear forces for Moments
for n = 2:788
    Moment_Y_Shaft1(n) = Moment_Y_Shaft1(n-1) + (Force_Z_Shaft1(n-1)+Force_Z_Shaft1(n))/2 * 0.01;
    Moment_Z_Shaft1(n) = Moment_Z_Shaft1(n-1) + (Force_Y_Shaft1(n-1)+Force_Y_Shaft1(n))/2 * 0.01;
end

% Shear, Moment and Torque Calculations for Intermediate Shaft
% Starting Analysis for Shaft 2 (Intermediate Shaft)
r_g = Dg/2; % radius of gear G
Fy = Cy;
Fz = Cz;
Tf = Tc*GR_CF;
Tg = Tf;
Gz = 100.19;
Gy = Gz * tan(theta2);
Dy = (-1.5625*Gy-6.0625*Fy)/7.3125; % Moment at E
Ey = -Gy-Fy-Dy;
Dz = (-1.5625*Gz-6.0625*Fz)/7.3125; % Length
Ez = -Gz-Fz-Dz; % testing
x2 = [0,1.25,1.25+4.5,1.25+4.5+1.5625]; % Length
Force2_y = [Dy,Dy+Fy,Dy+Fy+Gy,Dy+Fy+Gy+Ey];
Force2_z = [Dz,Dz+Fz,Dz+Fz+Gz,Dz+Fz+Gz+Ez];
Torque2 = [0,Tf,Tf-Tg,Tf-Tg];
Shaft2 = linspace(0,7.3125,732); % initializing variables for plotting
Force_Y_Shaft2 = zeros(size(Shaft2)); Force_Z_Shaft2 = zeros(size(Shaft2)); 
Moment_Y_Shaft2 = zeros(size(Shaft2)); Moment_Z_Shaft2 = zeros(size(Shaft2));
Torque_Shaft2 = zeros(size(Shaft2));

% inputing the values of Forces as a function of x
Torque_Shaft2((x2(2)*100+1):(x2(4)*100)) = Torque2(2);
Force_Y_Shaft2(1:(x2(2)*100)) = Force2_y(1);
Force_Y_Shaft2((x2(2)*100+1):(x2(3)*100)) = Force2_y(2); 
Force_Y_Shaft2((x2(3)*100+1):((x2(4)*100))) = Force2_y(3);
Force_Y_Shaft2(round(x2(4)*100)) = Force2_y(4);
Force_Z_Shaft2(1:(x2(2)*100)) = Force2_z(1);
Force_Z_Shaft2((x2(2)*100+1):(x2(3)*100)) = Force2_z(2); 
Force_Z_Shaft2((x2(3)*100+1):((x2(4)*100))) = Force2_z(3);
Force_Z_Shaft2(round(x2(4)*100)) = Force2_z(4);
% integrating Shear forces for Moments
for n = 2:732
    Moment_Y_Shaft2(n) = Moment_Y_Shaft2(n-1) + (Force_Z_Shaft2(n-1)+Force_Z_Shaft2(n))/2 * 0.01;
    Moment_Z_Shaft2(n) = Moment_Z_Shaft2(n-1) + (Force_Y_Shaft2(n-1)+Force_Y_Shaft2(n))/2 * 0.01;
end

% Shear, Moment and Torque Calculations for Output Shaft
% Starting Analysis for Shaft 3 (Output Shaft)
Tj = GR_JG*Tg;
Jy = Gy;
Jz = Gz;
Tk = Tj;
Iy = 1.5625/3.0625*Jy;
Hy = Jy - Iy;
Iz = -1.5625/3.0625*Jz; % Length
Hz = -Iz - Jz;
x3 = [0,1.5625,1.5625+1.5, 1.5625+1.5+2.75]; % Length
Force3_y = [Hy,Hy-Jy,Hy-Jy+Iy,Hy-Jy+Iy];
Force3_z = [Hz,Hz+Jz,Hz+Jz+Iz,Hz+Jz+Iz];
Torque3 = [0,Tj,Tj,Tj-Tk];
Shaft3 = linspace(0,58125,582); % initializing variables for plotting length
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
for n = 2:582
    Moment_Y_Shaft3(n) = Moment_Y_Shaft3(n-1) + (Force_Z_Shaft3(n-1)+Force_Z_Shaft3(n))/2 * 0.01;
    Moment_Z_Shaft3(n) = Moment_Z_Shaft3(n-1) + (Force_Y_Shaft3(n-1)+Force_Y_Shaft3(n))/2 * 0.01;
end
% Fatigue Analysis
% Fatigue Calculations for Input Shaft
% Can vary the lengths here for moment calculations
M1 = sqrt((Bz*1.5)^2+(By*1.5)^2); % Square sum of Moment for shaft 1
T1 = Tc; % Torque at Gear 
[dd1, DD1, r1, ~, ~, ~, ~, ~, ny1] = shaft(M1, T1, sut, sy);
fprintf('Shaft 1 Yielding FOS = %f\n',ny1)
fprintf('Shaft 1 Fatigue FOS = %f\n', 2)
fprintf('Shaft 1 filet radius = %f [in]\n', r1)

% Fatigue Calculations for Intermediate Shaft
% Can vary the lengths here for moment calculations
M2 = sqrt((Fz*1.25)^2+(Fy*1.25)^2); % Square sum of Moment for shaft 2
T2 = Tf; % Torque at Gear
[dd2, DD2, r2, ~, ~, ~, ~, ~, ny2] = shaft(M2, T2, sut, sy);
fprintf('Shaft 2 Yielding FOS = %f\n',ny2)
fprintf('Shaft 2 Fatigue FOS = %f\n', 2)
fprintf('Shaft 2 filet radius = %f [in]\n', r2)

% Fatigue Calculations for Output Shaft
% Can vary the lengths here for moment calculations
M3 = sqrt((Jz*1.5625)^2+(Jy*1.5625)^2); % Square sum of Moment for shaft 3
T3 = Tj; % Torque at Gear
[dd3, DD3, r3, ~, ~, ~, ~, ~, ny3] = shaft(M3, T3, sut, sy);
fprintf('Shaft 3 Yielding FOS = %f\n',ny3)
fprintf('Shaft 3 Fatigue FOS = %f\n', 2)
fprintf('Shaft 3 filet radius = %f [in]\n', r3)

% Deflection calculations
Bearing_Slope_max = 0.0005; % radians
SpurGear_Deflect_max = 0.003; % inches

%Moment_avg_Shaft1 = sqrt(Moment_Y_Shaft1.^2 + Moment_Z_Shaft1.^2);
%Moment_avg_Shaft2 = sqrt(Moment_Y_Shaft2.^2 + Moment_Z_Shaft2.^2);
%Moment_avg_Shaft3 = sqrt(Moment_Y_Shaft3.^2 + Moment_Z_Shaft3.^2);

% Can change the distances of the shaft here, denoted as x
d1 = [dd1/1.1, dd1, DD1, dd1]; x_shaft1 = [0, 1.5, 3, 4.625, 8.625];
d2 = [dd2/1.1, dd2, DD2, dd2, dd2/1.1]; x_shaft2 = [0,1,2.5,5.5,7,8.625];
d3 = [dd3, DD3, dd3, dd3/1.1]; x_shaft3 = [0,1.625,3.125,4.625,6.625];

Inertia1 = pi/16 .* d1.^4; Inertia2 = pi/16 .* d2.^4; Inertia3 = pi/16 .* d3.^4;

Deflection1_Y = zeros(size(Shaft1)); Deflection1_Z = zeros(size(Shaft1)); Inertia_Shaft1 = ones(size(Shaft1)) * Inertia1(1);
Deflection2_Y = zeros(size(Shaft2)); Deflection2_Z = zeros(size(Shaft2)); Inertia_Shaft2 = ones(size(Shaft2)) * Inertia2(1);
Deflection3_Y = zeros(size(Shaft3)); Deflection3_Z = zeros(size(Shaft3)); Inertia_Shaft3 = ones(size(Shaft3)) * Inertia3(1);

Theta1_Y = zeros(size(Shaft1)); Theta1_Z = zeros(size(Shaft1));
Theta2_Y = zeros(size(Shaft1)); Theta2_Z = zeros(size(Shaft1));
Theta3_Y = zeros(size(Shaft1)); Theta3_Z = zeros(size(Shaft1));

for i = 1:size(d1,2)
    Inertia_Shaft1(round((x_shaft1(i)*100)+1):(x_shaft1(i+1)*100)) = Inertia1(i);
end
for i = 1:size(d2,2)
    Inertia_Shaft2(round((x_shaft2(i)*100)+1):(x_shaft2(i+1)*100)) = Inertia2(i);
end
for i = 1:size(d3,2)
    Inertia_Shaft3(round((x_shaft3(i)*100)+1):(x_shaft3(i+1)*100)) = Inertia3(i);
end

for n = 2:size(Shaft1,2) 
    Theta1_Y(n) = Theta1_Y(n-1) + 1./(E*Inertia_Shaft1(n))*(Moment_Z_Shaft1(n-1) + Moment_Z_Shaft1(n))/2 * 0.01;
    Theta1_Z(n) = Theta1_Z(n-1) + 1./(E*Inertia_Shaft1(n))*(Moment_Y_Shaft1(n-1) + Moment_Y_Shaft1(n))/2 * 0.01;
end
for n = 2:size(Shaft2,2) 
    Theta2_Y(n) = Theta2_Y(n-1) + 1./(E*Inertia_Shaft2(n))*(Moment_Z_Shaft2(n-1) + Moment_Z_Shaft2(n))/2 * 0.01;
    Theta2_Z(n) = Theta2_Z(n-1) + 1./(E*Inertia_Shaft2(n))*(Moment_Y_Shaft2(n-1) + Moment_Y_Shaft2(n))/2 * 0.01;
end
for n = 2:size(Shaft3,2) 
    Theta3_Y(n) = Theta3_Y(n-1) + 1./(E*Inertia_Shaft3(n))*(Moment_Z_Shaft3(n-1) + Moment_Z_Shaft3(n))/2 * 0.01;
    Theta3_Z(n) = Theta3_Z(n-1) + 1./(E*Inertia_Shaft3(n))*(Moment_Y_Shaft3(n-1) + Moment_Y_Shaft3(n))/2 * 0.01;
end

for n = 2:size(Shaft1,2) 
    Deflection1_Y(n) = Deflection1_Y(n-1) + (Theta1_Y(n-1) + Theta1_Y(n))/2 * 0.01;
    Deflection1_Z(n) = Deflection1_Z(n-1) + (Theta1_Z(n-1) + Theta1_Z(n))/2 * 0.01;
end
for n = 2:size(Shaft2,2) 
    Deflection2_Y(n) = Deflection2_Y(n-1) + (Theta2_Y(n-1) + Theta2_Y(n))/2 * 0.01;
    Deflection2_Z(n) = Deflection2_Z(n-1) + (Theta2_Z(n-1) + Theta2_Z(n))/2 * 0.01;
end
for n = 2:size(Shaft3,2) 
    Deflection3_Y(n) = Deflection3_Y(n-1) + (Theta3_Y(n-1) + Theta3_Y(n))/2 * 0.01;
    Deflection3_Z(n) = Deflection3_Z(n-1) + (Theta3_Z(n-1) + Theta3_Z(n))/2 * 0.01;
end

Theta1 = sqrt(Theta1_Y.^2 + Theta2_Z.^2);
Theta2 = sqrt(Theta2_Y.^2 + Theta2_Z.^2);
Theta3 = sqrt(Theta3_Y.^2 + Theta3_Z.^2);

Deflection1 = sqrt(Deflection1_Y.^2 + Deflection1_Z.^2);
Deflection2 = sqrt(Deflection2_Y.^2 + Deflection2_Z.^2);
Deflection3 = sqrt(Deflection3_Y.^2 + Deflection3_Z.^2);

Critical_Gear_deflection = 0.003; % inches
Critical_Bearing_Slope = 0.0005; % radians

if Deflection1(size(Deflection1,2)) < Critical_Gear_deflection
    fprintf("Shaft 1 has a Maximum Deflection less than " + Critical_Gear_deflection + " in at all points\n")
end
if Deflection2(size(Deflection2,2)) < Critical_Gear_deflection
    fprintf("Shaft 2 has a Maximum Deflection less than " + Critical_Gear_deflection + " in at all points\n")
end
if Deflection3(size(Deflection3,2)) < Critical_Gear_deflection
    fprintf("Shaft 3 has a Maximum Deflection less than " + Critical_Gear_deflection + " in at all points\n")
end

%if Theta1(size(Theta1,2)) < Critical_Bearing_Slope
%    fprintf("Shaft 1 has a Maximum Slope less than " + Critical_Bearing_Slope + " radians at all points\n")
%end
%if Theta1(size(Theta2,2)) < Critical_Bearing_Slope
%    fprintf("Shaft 2 has a Maximum Slope less than " + Critical_Bearing_Slope + " radians at all points\n")
%end
%if Theta1(size(Theta3,2)) < Critical_Bearing_Slope
%    fprintf("Shaft 3 has a Maximum Slope less than " + Critical_Bearing_Slope + " radians at all points\n")
%end

%% Plots for Shaft 1
close all;
figure(1)
plot(Shaft1,Torque_Shaft1,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Torque (in-lbs)")
title("Torque")

figure(2)
plot(Shaft1,Force_X_Shaft1,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Axial Force (lbs)")
title("Axial Force")

figure(3)
plot(Shaft1,Force_Y_Shaft1,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Shear Force Y(lbs)")
title("Shear Force Y")

figure(4)
plot(Shaft1,Force_Z_Shaft1,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Shear Force Z(lbs)")
title("Shear Force Z")

figure(5)
plot(Shaft1,Moment_Y_Shaft1,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Moment Y(in-lbs)")
title("Moment Y")

figure(6)
plot(Shaft1,Moment_Z_Shaft1,"k-","LineWidth",5)
xlabel("Inches from Input to Output")
ylabel("Moment Z(in-lbs)")
title("Moment Z")

%% Plots for Shaft 2
close all;
figure(1)
plot(Shaft2,Torque_Shaft2,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Torque (in-lbs)")
title("Torque")

figure(2)
plot(Shaft2,Force_Y_Shaft2,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Shear Force Y(lbs)")
title("Shear Force Y")

figure(3)
plot(Shaft2,Force_Z_Shaft2,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Shear Force Z(lbs)")
title("Shear Force Z")

figure(4)
plot(Shaft2,Moment_Y_Shaft2,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Moment Y(in-lbs)")
title("Moment Y")

figure(5)
plot(Shaft2,Moment_Z_Shaft2,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Moment Z(in-lbs)")
title("Moment Z")

%% Plots for Shaft 3
close all;
figure(1)
plot(Shaft3,Torque_Shaft3,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Torque (in-lbs)")
title("Torque")

figure(2)
plot(Shaft3,Force_Y_Shaft3,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Shear Force Y(lbs)")
title("Shear Force Y")

figure(3)
plot(Shaft3,Force_Z_Shaft3,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Shear Force Z(lbs)")
title("Shear Force Z")

figure(4)
plot(Shaft3,Moment_Y_Shaft3,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Moment Y(in-lbs)")
title("Moment Y")

figure(5)
plot(Shaft3,Moment_Z_Shaft3,"k-","LineWidth",5)
xlabel("Inches")
ylabel("Moment Z(in-lbs)")
title("Moment Z")
