function [Force, Moment] = Force_Moment(rho, Aircraft, AeroMatrices, STATE, CONTROL)

%% Make vector X_tilde for inviscid part
X_tilde = Make_X_tilde_QuadAir1_3(STATE, CONTROL, 2, 1, Aircraft);

%% Compute inviscid aerodynamic force and momemt about CG
Fx_inv = rho*X_tilde'*AeroMatrices.Px*X_tilde;    % (N)
Fy_inv = rho*X_tilde'*AeroMatrices.Py*X_tilde;    % (N)
Fz_inv = rho*X_tilde'*AeroMatrices.Pz*X_tilde;    % (N)

Mx_inv = rho*X_tilde'*AeroMatrices.Qx*X_tilde; % (N.m)
My_inv = rho*X_tilde'*AeroMatrices.Qy*X_tilde; % (N.m)
Mz_inv = rho*X_tilde'*AeroMatrices.Qz*X_tilde; % (N.m)

%% Calculate viscous force and moment
F_visc = zeros(1,3);
M_visc = zeros(1,3);
strips = size(AeroMatrices.Strips,1);
r_Bo_CG = Aircraft.r_bo_CG;
for s = 1:strips
    % point of application of viscous force
    r_P_Bo_Std = AeroMatrices.Strip.quarterChordPoint(s,:);
    r_P_CG = r_Bo_CG + r_P_Bo_Std;
    
    Ss = AeroMatrices.Strip.chord(s)*AeroMatrices.Strip.width(s);
    
    % calculate local velocity vector
    Vel_inv = (AeroMatrices.VelMat_Strips(:,:,s)*X_tilde)';
    
    % compute total circulation for strip
    Gamma_Ti = AeroMatrices.Strip.width(s)*(AeroMatrices.Gamma_total_strip(s,:)*X_tilde);
    
    % vector projections for local velocity and alfa
    x_vec = [1 0 0];
    b_vec = AeroMatrices.Strip.bound(s,:);
    n_vec = cross(x_vec,b_vec);
    n_vec = n_vec/norm(n_vec);
    b_vec = cross(n_vec, x_vec);
    
    Vs = Vel_inv - dot(Vel_inv, b_vec)*b_vec;
    Cl_s = 2*Gamma_Ti/(Ss*norm(Vs));
    
    wing = AeroMatrices.Strip.wing(s);
    part = AeroMatrices.Strip.part(s);
    
    cd_0 = Aircraft.cd_0(wing,part);
    cd_1 = Aircraft.cd_1(wing,part);
    cd_2 = Aircraft.cd_2(wing,part);
    Cd_s = cd_0 + cd_1*Cl_s + cd_2*Cl_s^2;
    
    % Calculate viscous drag force & moment
    Ds = 0.5*rho*Ss*Cd_s*norm(Vs)*Vs;
    F_visc = F_visc + Ds;
    M_visc = M_visc + cross(r_P_CG,Ds);
    
end

%% Add inviscid and viscous contriubutions to get net force and moment about CG
Force  = [Fx_inv, Fy_inv, Fz_inv] + F_visc;
Moment = [Mx_inv, My_inv, Mz_inv] + M_visc;