
V = 20;
alfa = 0*pi/180;
beta = 0*pi/180;
%% 1. Define flight condition
rho  = 1.225; % (kg/m^3) air density
u_CG = V*cos(alfa)*cos(beta);    % (m/s) Body X velocity of CG, in Standard Frame (X: front, Y: right side, Z: down)
v_CG = V*sin(beta);     % (m/s) Body X velocity of CG, in Standard Frame (X: front, Y: right side, Z: down)
w_CG = V*sin(alfa)*cos(beta);     % (m/s) Body Z velocity of CG, in Standard Frame (X: front, Y: right side, Z: down)
p    = 0.0;   % (rad/s) roll rate
q    = 0.0;   % (rad/s) pitch rate
r    = 0.0;   % (rad/s) yaw rate
da   =  0*pi/180; % (rad) aileron command deflection
de   = 0*pi/180; % (rad) elevator command deflection
dr   =  0*pi/180; % (rad) rudder command deflection
STATE   = [u_CG v_CG w_CG p q r];  % vehicle state vector
CONTROL = [da de dr];              % vehicle control surface vector,

[Force, Moment] = Force_Moment(rho, Aircraft, AeroMatrices, STATE, CONTROL)


S = 14.9;
V = norm([u_CG, v_CG, w_CG]);
qS = 0.5*rho*V^2*S;

CL = -Force(3)/(qS)
CD = (-Force(3)*sin(alfa)-Force(1)*cos(alfa))/(qS)

CL_CD = CL/CD