


Ds_total = zeros(1,3);
strips = size(AeroMatrices.Strips,1);
r_Bo_CG = Aircraft.r_bo_CG;
for s = 1:strips
    % point of application of viscous force
    r_P_Bo_Std = AeroMatrices.Strip.quarterChordPoint(s,:);
    r_P_CG = r_Bo_CG + r_P_Bo_Std;
    
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
    
    Vs = Vel_inv - dot(Vel_inv, b_vec)*Vel_inv;
    Cl_s = 2*Gamma_Ti/(AeroMatrices.Strip.chord(s)*AeroMatrices.Strip.width(s)*norm(Vs))
    Cd_s = 0.01 + 0.05*Cl_s^2;
    % Calculate viscous drag force
    Ds = 0.5*rho*AeroMatrices.Strip.chord(s)*AeroMatrices.Strip.width(s)*Cd_s*norm(Vs)*Vs;
    Ds_total = Ds_total + Ds;
end

Ds_total