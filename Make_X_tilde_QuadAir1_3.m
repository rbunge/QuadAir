function X_tilde = Make_X_tilde_QuadAir1_3(STATE, CONTROL, StateOption, ControlOption, acraft_dsgn)

%% Constructs the vector X_tilde given the current state and control deflections, as
%% required by the P and Q matrices outputted by QuadAir1.3

% Allows for 2 different ways of inputting states and control deflections
% 
% StateOption: 1- Alfa & Beta wind angles 
%              2- Velocity components of CG
% 
% ControlOption: 1- CONTROL indicates deflection of groups of surfaces, according to 
% acraft_dsgn.con_surf_group 
%                2- CONTROL defines deflection of individial surfaces
% 
% u_CG = x-component wind relative velocity of CG expressed Body Fixed Standard Axis
% v_CG = y-component wind relative velocity of CG expressed Body Fixed Standard Axis
% w_CG = z-component wind relative velocity of CG expressed Body Fixed Standard Axis
% p = x-component of angular velocity expressed in Body Fixed Standard Axis
% q = y-component of angular velocity expressed in Body Fixed Standard Axis
% r = z-component of angular velocity expressed in Body Fixed Standard Axis

if StateOption == 1
    Vel_inf = STATE(1,1);
    alfa    = STATE(1,2);
    beta    = STATE(1,3);
    u_CG    = Vel_inf*cos(alfa)*cos(beta);    
    v_CG    = Vel_inf*sin(beta);             
    w_CG    = Vel_inf*sin(alfa)*cos(beta);   
else if StateOption == 2
        u_CG    = STATE(1,1);    
        v_CG    = STATE(1,2);   
        w_CG    = STATE(1,3);    
    end
end

p      = STATE(1,4);                    
q      = STATE(1,5);                     
r      = STATE(1,6);                         

X = [u_CG; v_CG; w_CG; p; q; r]; 
X_tilde(1:6,1) = X;

if ControlOption == 1
    for i = 1:length(acraft_dsgn.con_surf_group)
        group_nbr = abs(acraft_dsgn.con_surf_group(i));
        Sign = sign(acraft_dsgn.con_surf_group(i));
        con_surf_dflct(i) = CONTROL(1,group_nbr)*Sign;
    end
else if ControlOption == 2
        con_surf_dflct = CONTROL;
    end
end

% Going over each of the control surfaces
for j = 1:length(acraft_dsgn.con_surf_group)
    c_th = cos(con_surf_dflct(j));
    s_th = sin(con_surf_dflct(j));
    
    X_tilde(6+12*(j-1)+1:6+12*(j-1)+12,1) = [X*c_th
                                             X*s_th];
end
