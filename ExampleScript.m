%% Example - compute aerodynamic force and moment for a Cessna 152, using QuadAir
clear all
%% 1. Define Aicraft geometry and discretization: run Cessna152Example.m
Cessna152Example

%% 2. Run QuadAir to pre-compute the aerodynamic matrices
[Aircraft, AeroMatrices] = QuadAir1_3(Aircraft, geo_disc);

%% 3. Define flight condition
V = 20;
alfa = 5*pi/180;
beta = 0*pi/180;

rho  = 1.225; % (kg/m^3) air density
u_CG = V*cos(alfa)*cos(beta);    % (m/s) Body X velocity of CG, in Standard Frame (X: front, Y: right side, Z: down)
v_CG = V*sin(beta);     % (m/s) Body X velocity of CG, in Standard Frame (X: front, Y: right side, Z: down)
w_CG = V*sin(alfa)*cos(beta);     % (m/s) Body Z velocity of CG, in Standard Frame (X: front, Y: right side, Z: down)
p    = 0.0;   % (rad/s) roll rate
q    = 0.0;   % (rad/s) pitch rate
r    = 0.0;   % (rad/s) yaw rate
da   =  0*pi/180; % (rad) aileron command deflection
de   = -2*pi/180; % (rad) elevator command deflection
dr   =  0*pi/180; % (rad) rudder command deflection
STATE   = [u_CG v_CG w_CG p q r];  % vehicle state vector
CONTROL = [da de dr];              % vehicle control surface vector,

%% 4. Calculate force and moment at flight condition
[Force, Moment] = Force_Moment(rho, Aircraft, AeroMatrices, STATE, CONTROL);

fprintf('\n**********************')
fprintf('\nAerodynamic Force and Moment in Standard NED frame\n')

fprintf('\n***********************')
fprintf('\nNet aerodynamic force: (N)')
fprintf('\n***********************\n')
Fx = Force(1)
Fy = Force(2)
Fz = Force(3)


fprintf('\n***********************')
fprintf('\nNet aerodynamic moment about CG: (N.m)')
fprintf('\n***********************\n')
Mx_CG = Moment(1)
My_CG = Moment(2)
Mz_CG = Moment(3)



%% 5. Post-process force and moment to calculate derived data 
% Examples: non-dimensional coefficients, dertivatives (by doing multiple evaluations at perturbed conidtions), etc.
S = 14.9;
b = 10.2;
c = 1.5;

V = norm([u_CG, v_CG, w_CG]);
qS = 0.5*rho*V^2*S;

fprintf('\n***********************')
fprintf('\nExamples of non-dimensional coefficients')
fprintf('\n***********************\n')
CL = (-Fz*cos(alfa) + Fx*sin(alfa))/(qS)
CD = (-Fz*sin(alfa) - Fx*cos(alfa))/(qS)
Cl = Mx_CG/(qS*b)
Cm = My_CG/(qS*c)
Cn = Mz_CG/(qS*b)
CL_CD = CL/CD

% Showing aircraft geometry
close all
theta_z = -30; theta_x = 20; 
Aircraft_plotterV(Aircraft, theta_z, theta_x);

theta_z = 0; theta_x = 0; 
Aircraft_plotterV(Aircraft, theta_z, theta_x);

theta_z = 90; theta_x = 0; 
Aircraft_plotterV(Aircraft, theta_z, theta_x);

theta_z = 0; theta_x = 90; 
Aircraft_plotterV(Aircraft, theta_z, theta_x);
commandwindow