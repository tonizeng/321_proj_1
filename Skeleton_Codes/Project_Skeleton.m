%% Project Skeleton Code

clear; clc; close all;

%% Part 2 - Force and Moment Calculation

%%initial parameter: unit: m, degree, rad/sec
r1 = 0.084; 
r2 = 0.36;
r3 = 1.2;
r6 = 0.6;
m2 = 0.07634;
m6 = 0.12723;
m3 = 0.25447;
m4 = 5;
m5 = 5;
b2 = r2/2;
b6 = r6/2;
b3 = r3/2;

M12_list = [];
theta2_list = [];
Fs_list = [];  % shaking force
alpha_s_list = []; % direction of a shaking force
Ms_list =[]; % Shaking moment
%Forces list
F23_list = [];
F12_list = [];
F34_list = [];
F16_list = [];
F56_list = [];
F14_list = [];
N35_list = [];
%angles list
F12_alpha = [];
F23_alpha = [];
F34_alpha = [];
F16_alpha = [];
F56_alpha = [];

for theta2 = 0.01:0.1:4*pi

dtheta2 = 2;
ddtheta2 = 0; 

%% Part 1- Calculations for kinematic variables, caculated based on loop closure eqn

%equations for kinematic values
theta3 = pi - asin((r2.*sin(theta2)-r1)./r3);
r4 = r2.*cos(theta2) - r3.*cos(theta3);
theta5 = theta3; % eqn 15
theta6 = pi - asin((r1.*cos(theta5)-r4.*sin(theta5))./r6) + theta5; %eqn 18
r5 = (r6.*cos(theta6)-r4)./cos(theta5);

%% Take time derivative of loop eqn (d/dt) 
% and solve them for dtheta3, dtheta5 & dr6
% and the same for the second derivatives. 
% add r5dot, theta6dot, r5 double dot, theta6 double dot

% First derivates of r4 and theta 3
dtheta3 = (r2.*dtheta2.*cos(theta2))./(r3.*cos(theta3));
dr4 = (r2.*dtheta2.*sin(theta3 - theta2))./cos(theta3);

% Second derivates of r4 and theta 3
ddtheta3 = ((r3.*(dtheta3.^2).*sin(theta3))-(r2.*(dtheta2.^2).*sin(theta2)))./(r3.*cos(theta3));
ddr4 = (r3.*ddtheta3.*sin(theta3)) + (r3.*(dtheta3.^2).*cos(theta3)) - (r2.*(dtheta2.^2).*cos(theta2));

% First derivatives of r5 and theta 6
dtheta5 = dtheta3;
dtheta6 = ((r5.*dtheta5)-dr4.*sin(theta5))./(r6.*cos(theta5-theta6));
dr5 = ((r6.*dtheta6.*cos(theta6))-(r5.*dtheta5.*cos(theta5)))./(sin(theta5));

% Second derivatives of r5 and theta 6
ddtheta5 = ddtheta3;
ddtheta6 = ((-ddr4.*sin(theta5))+(2.*dr5.*dtheta5)+(r5.*ddtheta5)-(r6.*(dtheta6.^2).*sin(theta5-theta6)))./(r6.*cos(theta5-theta6));
ddr5 = (((-2.*dr5.*dtheta5-r5.*ddtheta5).*cos(theta5))+(r5.*(dtheta5.^2).*sin(theta5))+(r6.*ddtheta6.*cos(theta6))-(r6.*(dtheta6.^2).*sin(theta6)))./(sin(theta5));

a_coriolis = abs((2.*(r6.*dtheta6.*cos(theta6)-r5.*(dtheta5).*cos(theta5)).*(r2.*dtheta2.*cos(theta2)))./(r3.*cos(theta3).*sin(theta5)));

%given parameters

beta3 = 2*pi-theta3;
dbeta3 = -1*theta3;
ddbeta3 = -1*ddtheta3;

%inertias
inertiaG3 = (1/12)*m3*r3.^2;
inertiaA2 = (1/3)*m2*r2.^2;
inertiaA6 = (1/3)*m6*r6.^2;
%shortcuts
rGC = r5 - b3;
rGB = r3./2;
rDG = r3./2;
%trig substitutions

A33 = -r2.*sin(theta2);
A43 = r2.*cos(theta2);
A134 = sin(theta6);
A135 = cos(theta6);
A36 = rGB.*sin(beta3);
A46 = rGB.*cos(beta3);
A1311 = -sin(theta6);
A1312 = -cos(theta6);
A136 = -rGC.*cos(theta6 + beta3);
A66 = -rDG.*cos(beta3);
A56 = -rDG.*sin(beta3);
A99 = -r6.*sin(theta6);
A109 = r6.*(cos(theta6));

%accelerations for forces
ag2x = b2.*(-dtheta2.^2).*cos(theta2);
ag2y = b2.*(-dtheta2.^2)*sin(theta2);

ag3x = (b3.*(-dbeta3.^2).*cos(2*pi - beta3)) + ag2x; % huh
ag3y = (b3.*(ddbeta3).*sin(2*pi - beta3)) + ag2y; % huh

ag6x = b6.*(-ddtheta6.*sin(theta6)-(dtheta6.^2).*cos(theta6));
ag6y = b6.*(ddtheta6.*cos(theta6)-(dtheta6.^2).*sin(theta6));

ag4x = ddr4;
ag4y = 0;

ag5x = ag6x.*2;
ag5y = ag6y.*2;

% and so on    

    B = [m2.*ag2x;
        m2.*ag2y;
        0;
        m3.*ag3x;
        m3.*ag3y;
        inertiaG3.*ddtheta3; % wtf
        m6.*ag6x;
        m6.*ag6y;
        inertiaA6.*ddtheta6;
        -m4.*ag4x;
        -m5.*ag5x;
        -m5.*ag5y;
        0;
    ];
    
    A = [
    -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, A33, A43, 0, 0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, A134;
    0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, A135;
    0, 0, A36, A46, A56, A66, 0, 0, 0, 0, 0, 0, A136;
    0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, A99, A109, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, A1311;
    0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, A1312;
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0;
    ];

    x = A\B % Ax = B, solution for x; note that in MATLAB: A\B = B/A
    
    
    % M12:
    M12 = x(12);
    M12_list = [M12_list; M12];
    
    F12x = x(1);
    F12y = x(2);
    F23x = x(3);
    F23y = x(4);
    F34x = x(5);
    F34y = x(6);
    F16x = x(7);
    F16y = x(8);
    F56x = x(9);
    F56y = x(10);
    F14y = x(11);
    N35 = x(13);
    
    % Magnitudes of all forces: 
    % Atan is defined on [-pi/2; pi/2]. 
    % This if clause will help to adjust the value of the angle 
    % to its true value:	
    F23_list = [F23_list; sqrt(F23x.^2+F23y.^2)];
    F12_list = [F12_list; sqrt(F12x.^2+F12y.^2)];
    F34_list = [F34_list; sqrt(F34x.^2+F34y.^2)];
    F16_list = [F16_list; sqrt(F16x.^2+F16y.^2)];
    F56_list = [F56_list; sqrt(F56x.^2+F56y.^2)];
    F14_list = [F14_list; F14y];
    N35_list = [N35_list; N35];

    % Directions of all forces:   
    %F23
    fx = F23x;
    fy = F23y;
    alpha_23 = atan(fy/fx);
    if fx < 0
        alpha_23 = alpha_23 + pi;
    end 
    F23_alpha = [F23_alpha; alpha_23];

    %F12
    fx = F12x;
    fy = F12y;
    alpha_12 = atan2(fy, fx);
    if fx < 0
        alpha_12 = alpha_12 + pi;
    end 
    F12_alpha = [F12_alpha; alpha_12];

    %F34
    fx = F34x;
    fy = F34y;
    alpha_34 = atan(fx\fy);
    if fx < 0
        alpha_34 = alpha_34 + pi;
    end 
    F34_alpha = [F34_alpha; alpha_34];

    %F16
    fx = F16x;
    fy = F16y;
    alpha_16 = atan(fx\fy);
    if fx < 0
        alpha_16 = alpha_16 + pi;
    end 
    F16_alpha = [F16_alpha; alpha_16];

    %F56
    fx = F56x;
    fy = F56y;
    alpha_56 = atan(fx\fy);
    if fx < 0
        alpha_56 = alpha_56 + pi;
    end 
    F56_alpha = [F56_alpha; alpha_56];

  
    % Collecting the values of theta2:
    theta2_list = [theta2_list; theta2];
     
   %}
    
end


% Regular and Polar plots:
% Might have to transpose the Force vectors for polar plot. Do so if needed
% Polar plot only works with radians so will have to do it accordingly
M12_list = M12_list(:);

figure (1)
plot(F23_list)
grid on;
title('M_{12} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('M12   unit: N-m')
hold on;

figure (2)
plot(theta2_list,F23_list)
grid on;
title('F_1_2 vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('M12   unit: N-m')
hold on;

% Convert degrees to the radians
theta2_rad = deg2rad(theta2_list);

%figure (4)
%polarplot(Fij_alpha,Fij_list)
%grid on;
%title('F_{ij} polar plot')

% and so on ...