%% Example - compute aerodynamic force and moment for a Cessna 152, using QuadAir

%% 1.Define aircraft geometry: run Cessna152.m
Cessna152Example

%% 2. Run QuadAir to pre-compute the P & Q aerodynamic matrices
[Aircraft, AeroMatrices] = QuadAir1_3(Aircraft, geo_disc);

fprintf('\nThese are the aerodynamic matrices, only need to be computed once!\n')
AeroMatrices.Px;
AeroMatrices.Py;
AeroMatrices.Pz;
AeroMatrices.Qx;
AeroMatrices.Qy;
AeroMatrices.Qz;

%% 3. Define flight condition
rho  = 1.225; % (kg/m^3) air density
u_CG = 20;    % (m/s) Body X velocity of CG, in Standard Frame (X: front, Y: right side, Z: down)
v_CG = 0;     % (m/s) Body X velocity of CG, in Standard Frame (X: front, Y: right side, Z: down)
w_CG = 2;     % (m/s) Body Z velocity of CG, in Standard Frame (X: front, Y: right side, Z: down)
p    = 0.0;   % (rad/s) roll rate
q    = 0.0;   % (rad/s) pitch rate
r    = 0.0;   % (rad/s) yaw rate
da   =  0*pi/180; % (rad) aileron command deflection
de   = 0*pi/180; % (rad) elevator command deflection
dr   =  10*pi/180; % (rad) rudder command deflection
STATE   = [u_CG v_CG w_CG p q r];  % vehicle state vector
CONTROL = [da de dr];              % vehicle control surface vector,

%% 4. Make vector X_tilde
X_tilde = Make_X_tilde_QuadAir1_3(STATE, CONTROL, 2, 1, Aircraft);

%% 5.a) Compute net aerodynamic force
fprintf('\n***********************')
fprintf('\nNet aerodynamic force: (N)')
fprintf('\n***********************\n')
Fx = rho*X_tilde'*AeroMatrices.Px*X_tilde    % (N)
Fy = rho*X_tilde'*AeroMatrices.Py*X_tilde    % (N)
Fz = rho*X_tilde'*AeroMatrices.Pz*X_tilde    % (N)

%% 5.b) Compute net aerodynamic moment about CG
fprintf('\n***********************')
fprintf('\nNet aerodynamic moment about CG: (N.m)')
fprintf('\n***********************\n')
Mx_CG = rho*X_tilde'*AeroMatrices.Qx*X_tilde % (N.m)
My_CG = rho*X_tilde'*AeroMatrices.Qy*X_tilde % (N.m)
Mz_CG = rho*X_tilde'*AeroMatrices.Qz*X_tilde % (N.m)



% % Showing aircraft geometry
% close all
% theta_z = -30; theta_x = 20; 
% Aircraft_plotterV(Aircraft, theta_z, theta_x);
% 
% theta_z = 0; theta_x = 0; 
% Aircraft_plotterV(Aircraft, theta_z, theta_x);
% 
% theta_z = 90; theta_x = 0; 
% Aircraft_plotterV(Aircraft, theta_z, theta_x);
% 
% theta_z = 0; theta_x = 90; 
% Aircraft_plotterV(Aircraft, theta_z, theta_x);
% commandwindow




