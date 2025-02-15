%% The following code plots all kinetic equations
clear; clc; close all;

%% Part 2 - Force and Moment Calculation
%% Lists to store calculated values for different parameters
% Moment & Angle lists
M12_list = []; % List to store moment M12 values
theta2_list = []; % List to store angular positions of link 2 (theta2)

% Shaking Force & Moment
Fs_list = [];  % List to store magnitudes of shaking forces
Fs_alpha = []; % List to store angles (directions) of shaking forces
Ms_list =[]; % List to store shaking moments

% Initialization orce components acting at different joints
F23_list = []; % List to store force between links 2 and 3
F12_list = []; % List to store force between base 1 and link 2
F34_list = []; % List to store force between links 3 and slider 4
F16_list = []; % List to store force between base 1 and link 6
F56_list = []; % List to store force between links 5 and 6
F14_list = []; % List to store force between base 1 and link 4
N35_list = []; % List to store normal force at slider 5

% Angles corresponding to the forces
F12_alpha = []; % Angle of force F12
F23_alpha = []; % Angle of force F23
F34_alpha = []; % Angle of force F34
F16_alpha = []; % Angle of force F16
F56_alpha = []; % Angle of force F56

%% Givens / Initial Parameters 
% Link Lengths (m)
r1 = 0.084; 
r2 = 0.36;
r3 = 1.2;
r6 = 0.6;

% Link Masses (kg)
m2 = 0.07634;
m6 = 0.12723;
m3 = 0.25447;

% Slider Masses (kg)
m4 = 5;
m5 = 5;

% Center of Gravity Distances (m)
b2 = r2/2;
b6 = r6/2;
b3 = r3/2;

%Inertias
inertiaG3 = (1/12)*m3.*r3.^2;
inertiaA2 = (1/3)*m2.*r2.^2;
inertiaA6 = (1/3)*m6.*r6.^2;

% Theta 2
dtheta2 = 2;
ddtheta2 = 0; 

for theta2 = 0:0.01:2*pi

% Calculations of kinematic variables from Part 1
theta3 = pi - asin((r2.*sin(theta2)-r1)./r3);
r4 = r2.*cos(theta2) - r3.*cos(theta3);
theta5 = theta3;% eqn 15
theta6 = pi - asin((r1.*cos(theta5)-r4.*sin(theta5))./r6) + theta5; %eqn 18
r5 = (r6.*cos(theta6)-r4)./cos(theta5);

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

% Beta 3
beta3 = 2*pi-theta3;
dbeta3 = -1*dtheta3;
ddbeta3 = -1*ddtheta3;

% Shortcuts for Link Components 
rGC = r5 - b3;
rGB = r3./2;
rDG = r3./2;

% Coefficients for the A matrix 
% Naming convention of the coefficients follows the x,y coordinates of the
% location in the A matrix
A33 = -r2.*sin(theta2);
A43 = r2.*cos(theta2);
A134 = -sin(beta3);
A135 = -cos(beta3);
A36 = rGB.*sin(beta3);
A46 = rGB.*cos(beta3);
A1311 = sin(beta3);
A1312 = cos(beta3);
A136 = rGC;
A66 = -rDG.*cos(beta3);
A56 = -rDG.*sin(beta3);
A99 = -r6.*sin(theta6);
A109 = r6.*(cos(theta6));

% Acceleration Equations for Forces
ag2x = b2.*-(dtheta2.^2).*cos(theta2);
ag2y = b2.*-(dtheta2.^2).*sin(theta2);

ag3x = (b3.*-(dbeta3.^2).*cos(2*pi - beta3)) + ag2x;
ag3y = (b3.*(ddbeta3).*sin(2*pi - beta3)) + ag2y; 

ag6x = b6.*(-ddtheta6.*sin(theta6)-(dtheta6.^2).*cos(theta6));
ag6y = b6.*(ddtheta6.*cos(theta6)-(dtheta6.^2).*sin(theta6));

ag4x = ddr4;
ag4y = 0;

ag5x = ag6x.*2;
ag5y = ag6y.*2;


% B Vector
    B = [m2.*ag2x;
        m2.*ag2y;
        0;
        m3.*ag3x;
        m3.*ag3y;
        inertiaG3.*ddtheta3; 
        m6.*ag6x;
        m6.*ag6y;
        inertiaA6.*ddtheta6;
        m4.*ag4x;
        m5.*ag5x;
        m5.*ag5y;
        0;
    ];
    
% A Matrix
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

    x = A\B; % Ax = B, solution for x;
    
    
    % Force & Moment Results from Solving the Matrix
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
    M12 = x(12);
    N35 = x(13);
    Fsx = F12x + F16x;
    Fsy = F14y + F12y + F16y;
    Ms = -M12 + F14y.*r4;

    % Store Force & Moment Magnitudes into List
    F23_list = [F23_list; sqrt(F23x.^2+F23y.^2)];
    F12_list = [F12_list; sqrt(F12x.^2+F12y.^2)];
    F34_list = [F34_list; sqrt(F34x.^2+F34y.^2)];
    F16_list = [F16_list; sqrt(F16x.^2+F16y.^2)];
    F56_list = [F56_list; sqrt(F56x.^2+F56y.^2)];
    F14_list = [F14_list; abs(F14y)];
    M12_list = [M12_list; abs(M12)];
    N35_list = [N35_list; abs(N35)];
    Fs_list = [Fs_list; sqrt(Fsx.^2+Fsy.^2)];
    Ms_list = [Ms_list; abs(Ms)];

    % Directions of Forces & Moment:   
    % F23
    fx = F23x;
    fy = F23y;
    alpha_23 = atan(fy/fx);
    if fx < 0
        alpha_23 = alpha_23 + pi;
    end 
    F23_alpha = [F23_alpha; alpha_23];

    % F12
    fx = F12x;
    fy = F12y;
    alpha_12 = atan2(fy, fx);
    if fx < 0
        alpha_12 = alpha_12 + pi;
    end 
    F12_alpha = [F12_alpha; alpha_12];

    % F34
    fx = F34x;
    fy = F34y;
    alpha_34 = atan(fx\fy);
    if fx < 0
        alpha_34 = alpha_34 + pi;
    end 
    F34_alpha = [F34_alpha; alpha_34];

    % F16
    fx = F16x;
    fy = F16y;
    alpha_16 = atan(fx\fy);
    if fx < 0
        alpha_16 = alpha_16 + pi;
    end 
    F16_alpha = [F16_alpha; alpha_16];

    % F56
    fx = F56x;
    fy = F56y;
    alpha_56 = atan(fx\fy);
    if fx < 0
        alpha_56 = alpha_56 + pi;
    end 
    F56_alpha = [F56_alpha; alpha_56];

    % Fs
    fx = Fsx;
    fy = Fsy;
    alpha_s = atan(fx\fy);
    if fx < 0
        alpha_s = alpha_s + pi;
    end 
    Fs_alpha = [Fs_alpha; alpha_s];

  
    % Collecting the values of theta2:
    theta2_list = [theta2_list; theta2];
end

% Plots
% M12 
subplot(3,4,1);
plot(theta2_list.*(180/pi), M12_list)
grid on;
title('M_{12} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('M12   unit: N-m')
xlim([0, 360])   
xticks(0:90:360) 
hold on;

% F23 
subplot(3,4,2);
plot(theta2_list.*(180/pi),F23_list)
grid on;
title('F_{23} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('F_{23}   unit: N')
xlim([0, 360])   
xticks(0:90:360) 
hold on;

% F12 
subplot(3,4,3);
plot(theta2_list.*(180/pi), F12_list)
grid on;
title('F_{12} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('F_{12}   unit: N')
xlim([0, 360])   
xticks(0:90:360) 
hold on;

% F34 
subplot(3,4,4);
plot(theta2_list.*(180/pi), F34_list)
grid on;
title('F_{34} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('F_{34}   unit: N')
xlim([0, 360])   
xticks(0:90:360) 
hold on;

% F16 
subplot(3,4,5);
plot(theta2_list.*(180/pi), F16_list)
grid on;
title('F_{16} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('F_{16}   unit: N')
xlim([0, 360])   
xticks(0:90:360) 
hold on;

% F56 
subplot(3,4,6);
plot(theta2_list.*(180/pi), F56_list)
grid on;
title('F_{56} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('F_{56}   unit: N')
xlim([0, 360])   
xticks(0:90:360) 
hold on;

% F14 
subplot(3,4,7);
plot(theta2_list.*(180/pi), F14_list)
grid on;
title('F_{14} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('F_{14}   unit: N')
xlim([0, 360])   
xticks(0:90:360) 
hold on;

% F35 
subplot(3,4,8);
plot(theta2_list.*(180/pi), N35_list)
grid on;
title('F_{35} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('F_{35}   unit: N')
xlim([0, 360])   
xticks(0:90:360) 
hold on

% Shaking Force Fs
subplot(3,4,9);
plot(theta2_list.*(180/pi), Fs_list)
grid on
title('F_{s} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('F_{s}   unit: N')
xlim([0, 360])   
xticks(0:90:360) 
hold on

% Shaking Moment Ms
subplot(3,4,10);
plot(theta2_list.*(180/pi), Ms_list)
grid on
title('M_{s} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('M_{s}   unit: N')
xlim([0, 360])   
xticks(0:90:360) 
hold on