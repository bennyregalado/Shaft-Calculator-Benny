function [dn, D, rn, kt, kts, Kf, Kfs, se,ny] = shaft(M, T, sut, sy)
% All parameters, inputs and outputs are in English Units
% Input:
% M: Maximum Moment on the shaft 
% T: Maximum Torque on the shaft 
% sut: Ultimate Tensile Strength of the material chosen 
% Output:
% D: Bigger diameter 
% d: smaller diameter
% r: fillet radius 
% Assumptions:
% a: factor for calculating ka, for machined/cold-drawn material 
% b: factor for calculating ka, for machined/cold-drawn material 
% n: factor of safety is 2
% D/d = 1.1
% r/d = .05
% Note: In the absence of information regarding 𝑘𝑑, 𝑘𝑒 and 𝑘𝑓, we approximate 
% them to be equal to 1.
% ka: surface condition modification factor 
% kb: size modification factor, assuming the diameter will be 0.11 <= d <= 2
% kc: load modification factor, assuming bending
% kd: temperature modification factor, assuming that the maximum temperature 
% is experienced for shaft (~482 F or 250 C)
% ke: reliability factor, assuming 99% reliability 
% kf: miscellaneous-effects modification factor, assumption of 1, is only added 
% in special cases. It comes from experience of design engineer.
% se_p: rotary-beam test specimen endurance limit, assuming sut <= 200 kpsi
% kt: stress concentration factor, using Table 7-1 an assumption is made
% kts: stress concentration factor for torsion, using Table 7-1 an assumption is made
% Kf: stress concentration due to fatigue from Moment
% Kfs: stress concentration due to fatigue from Torsion
% Fatigue stress concentration calculation factors
% sqrt_a
% sqrt_as
% q
% qs
n = 2;
a = 2.7;
b = -0.265;
ka = a*(sut*10^-3)^b;
do = .5; % This will be the initial diameter that will be used as a guess
kb = 0.897*(do)^(-0.107);
kc = 1;
TF = 482;
kd = 0.975+0.432*(10^-3)*TF-0.115*(10^-5)*TF^2+0.104*(10^-8)*TF^3-0.595*(10^-12)*TF^4;
ke = 0.814;
kf = 1;
se_p = 0.5*sut*10^-3;
% Calculation for the Endurance Limit
se = ka * kb * kc * kd * ke * kf * se_p;
kt = 1.7; % Used as an assumption
kts = 1.5; % Used as an assumption
Kf = kt;
Kfs = kts;
% Below is from the Distortion Energy Goodman Equation
sigma_ap = Kf*(32*M*10^(-3)/pi);
sigma_pp = sqrt(3)*Kfs*(32*T*10^(-3)/pi);
do = ((sigma_ap/(se)+sigma_pp/(sut))*n)^(1/3);
ro = .05*do;
D = 1.1*do;
% While loop for final calculations
z = 1 - (1 / 1.1);
sqrt_a = 0.246 - (3.08 * 10^-3) * (sut * 10^-3) + (1.51 * 10^-5) * (sut * 10^-3)^2 - (2.67 * 10^-8) * (sut * 10^-3)^3;
sqrt_as = 0.19 - (2.51 * 10^-3) * (sut * 10^-3) + (1.35 * 10^-5) * (sut * 10^-3)^2 - (2.67 * 10^-8) * (sut * 10^-3)^3;
error = 1;
error_1 = 0.001;
% While loop to calculate the final calculations
while error > error_1
    y = ro/do;
    x = 0.05/y;
    rn = y*do;
    % Bending Stress Concentration Factor
    c1 = 0.947+1.206*sqrt(x)-0.131*x;
    c2 = 0.022-3.405*sqrt(x)+0.915*x;
    c3 = 0.869-1.777*sqrt(x)-0.555*x;
    c4 = 0.810-0.422*sqrt(x)-0.260*x;
    kt = c1+c2*(z)+c3*(z)^2+c4*(z)^3;
    % Torsion Stress Concentration Factor
    c11 = 0.905+0.783*sqrt(x)-0.075*x;
    c22 = -0.437-1.969*sqrt(x)+0.553*x;
    c33 = 1.557+1.073*sqrt(x)-0.578*x;
    c44 = -1.061+0.171*sqrt(x)+0.086*x;
    kts = c11+c22*(z)+c33*(z)^2+c44*(z)^3;
    % Calculate fatigue stress concentration factors
    q = 1/(1+(sqrt_a/sqrt(rn)));
    qs = 1/(1+(sqrt_as/sqrt(rn)));
    % Calculate fatigue stress concentration
    Kf = 1+q*(kt-1);
    Kfs = 1+qs*(kts-1);
    kb = 0.897*(do)^(-0.107);
    % Calculation for the Endurance Limit
    se = ka*kb*kc*kd*ke*kf*se_p;
    % Calculate new diameter using the Distortion Energy Goodman Theory
    sigma_ap = Kf*(32*M*10^(-3)/pi);
    sigma_pp = sqrt(3)*Kfs*(32*T*10^(-3)/pi);
    dn = ((sigma_ap/(se)+sigma_pp/(sut))*n)^(1/3);
    D = 1.1*dn;
    rn = .05*dn;
    error = abs((do-dn)/dn);
    do = dn;
    ro=rn;
    ny = ((sy*10^-3)/(sigma_ap+sigma_pp));
end
end
