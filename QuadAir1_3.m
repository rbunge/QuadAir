%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% © Copyright Roberto A. Bunge, Dept. of Aero/Astro, Stanford University, 2010-2012.  
% Praise be to God!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes for the user:
% QuadAir is matlab implementation of the Compact VLM.  It allows for rapid recalculation of steady low speed 
% inviscid aerodynamic forces and moements produced by the lifting surfaces of fixed wing aircraft.  The user should:
%     1. Define the aircraft geometry and CG location in the format required by QuadAir. See commented Example below.
%     2. Run QuadAir and store the structure AeroMatrices, which contains the matrices required to rapidly 
%        re-calculate aerodynamic forces and moments.  
%     3. During the flght simualtion, define the flight condition (air density, velocity, angular velocity & control surface deflections) 
%     4. Calculate X_tilde, using "Make_X_tilde.m", 
%     5. Calculate forces and moments by doing the following calculation:
%         Fx = rho*X_tilde'*AeroMatrices.Px*X_tilde
%         Fy = rho*X_tilde'*AeroMatrices.Py*X_tilde
%         Fz = rho*X_tilde'*AeroMatrices.Pz*X_tilde
%         Mx = rho*X_tilde'*AeroMatrices.Qx*X_tilde
%         My = rho*X_tilde'*AeroMatrices.Qy*X_tilde
%         Mz = rho*X_tilde'*AeroMatrices.Qz*X_tilde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Develeopment notes: 
% OUTPUTS:
% P & Q matrices are inviscid forces and moment matrices, where these end up being expressed in Standard Frame and moments 
% taken with respect to CG.  
% We also output P_red & Q_red which leverage on grouped surfaces, in order to make computation time faster.
% We also output A, Cj, Sj matrices which can be used to compute Xtilde
% directly.  This works for the non-reduced P & Q matrices.
%
% UPDATES:
% - Dec 23 2011: Added matrices L and D so that X can be constructed directly from u_wind_rel_CG, v_wind_rel_CG, w_wind_rel_CG and also modified Q matrices 
% so that moments are automa% tically taken with respect to the CG.  For this the user has to specify acraft_dsgn.r_bo_CG which is the position of Bo 
% relative to CG expressed in Body Fixed Standard Axis.
% - Jan 19 2012: Added strip force capabilities, by computing the relation between panels and strips, and for each strip computing its relevant data, like 
% chord length and quarter chord position.  This allows for direct computation of the inviscid force acting on a given strip.  If the chordwise discretization
% is not the same for flapped and unflapped parts, the strip capabilities won't work properly.
% - Jan 23 2012: Added capabilities to produce reduced P and Q matrices in order to reduce the number of computations both in Force Moment evaluation and 
% Gradient evaluation.  This is done based on the fact that control surfaces are somtimes grouped such that thetai = thetaj or thetai = -thetaj, and also 
% on the fact that it is possible to approximate cos(theta) = 1, when theta
% is small.  The latter won't be implemented just yet.
% - Jan 27: added more strip information, like twsit angle.
% - Feb 1: added different transformation matrices for X, to convert from flight simulation X to internal 
% X required by QuadAir.  Also the A, C, S
% matrices for contstruction of Xtilde are ok, as per checking with
% Trsistan and checking the derivations on paper.  Added normal vectors as
% outputs too.  Added velocity vector matrices at all bound vor% tices, and at closest
% bound vortex to the quarter chord of every strip.  All of the these
% velocity vectors are expressed in Standard frame, and have to be
% calculated by premultiyplying the corresponding matrix with the Xtilde
% that QuadAir uses internally, i.e. that calculated with XGeoQinfBo. Alos  
% composed of Qinf
% - Apr 24: added correspondance between strips and control surfaces
% - Sept 10: added control surface torque and force
% - Sept 21: added matrix DD to allow X to be defined directly in CG_V_Std variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INPUT FILE EXAMPLE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% Aircraft geometry parameters
% % Aircraft.r_CG_bo = [-0.375, 0, 0]; % position of CG with respect to apex of first wing, in Standard Frame (X: front, Y: right side, Z: down)
% % Aircraft.symmetry  = [1
% %                       1
% %                       0];        % Indicates which wings have symetry about the X-Z plane
% % Aircraft.n_prt_wng = [2
% %                       1
% %                       1];        % Number of partitions per wing 
% % Aircraft.flapped   = [0 1
% %                       1 0
% %                       1 0];      % Indicates which partitions are flapped
% % Aircraft.flp_frac  = [0    0.2     
% %                       0.2  0
% %                       0.2  0];    % Indicates the chord fraction occupied by the flap at eah partition.  If partition is UNFLAPPED then set to 0.
% % Aircraft.spn       = [4.2  5.7
% %                       3    0
% %                       1.75 0];      % Span (including the symmetric part if exists) of each partition of each wing
% % Aircraft.root_chrd = [1.5
% %                       1
% %                       1];      % Root chord of each wing                                 
% % Aircraft.tpr_rto   = [1      0.75
% %                       0.7    0
% %                       0.6    0];         % Taper ratio of each partition of each wing
% % Aircraft.swp_angle = [0      0
% %                       1      0
% %                       30     0]*pi/180;  % Sweep angle of each partition of each wing 
% % Aircraft.dih_angle = [0      0
% %                       0      0
% %                       90     0]*pi/180;  % Dihedral angle of each partition of each wing 
% % Aircraft.xyz_000   = [0      0       0
% %                       5      0       -0.6
% %                       5      0       -0.6];     % XYZ Position of the apex of each the wing, in Geometric Frame (X: back, Y: right side, Z: up). First one should be (0, 0, 0), to make it the reference point
% % Aircraft.twst_ang  = [+2    +2    +2
% %                       -2    -2     0
% %                        0     0     0]*pi/180; % Angle of incidence of each station of each wing (NOTE: each wing has #partitions+1 stations)
% % Aircraft.airfoil   = {0 0;         % Sets the airfoil shape for each partition of each wing.  0 sets a flat plate.  User defined airfoils can be used
% %                       0 0;         % User defined airfoils can be used, if they are in the same folder as the path.  Example 'Ag37.dat'.  
% %                       0 0};      % Airfoils are defined as an X-Y column matrix going continuously from TE to LE and all the way back to TE. 
% % 
% % Aircraft.wng_con_surf            = [0];       % Indicating which wings are full control surfaces.  The whole wing is rotated.
% % Aircraft.wng_con_surf_axis_rot   = [0 1 0];   % Specifying the axis of rotation of each full wing control surface.
% % Aircraft.con_surf_group   = [1                % Indicates grouping of control surfaces and symmetric/anti-symmetric relation
% %                              1
% %                              2
% %                             -2
% %                              3];
% %                             
% % %% Geometric Discretization Parameters
% % % The structure "geo_disc" holds the relevant geometric disretization parameters.
% % % UNFLAPPED part
% % geo_disc.spn_div(:,:,1)   = [5 5 
% %                              4 0
% %                              4 0];      % Number of span-wise divisions for the UN-FLAPPED part of each partition of each wing
% % geo_disc.chrd_div(:,:,1)  = [4 3
% %                              4 0
% %                              4 0];      % Number of chordwise-wise divisions for the UN-FLAPPED part of each partition of each wing
% % % FLAPPED part
% % geo_disc.spn_div(:,:,2)   = geo_disc.spn_div(:,:,1);      % Set to be equal in FLAPPED and UN-FLAPPED parts.
% % geo_disc.chrd_div(:,:,2)  = [0 3
% %                              3 0
% %                              3 0];      % Number of chordwise-wise divisions for the FLAPPED part of each partition of each wing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acraft_dsgn, AeroMatrices, xyzG] = QuadAir1_3(acraft_dsgn, geo_disc)
fprintf('\n*******************************************************************\n')
fprintf('\nRUNNING QUADAIR .....\n')
% % tic
% Rotation matrix that transforms vectors expressed in Geometric Axis to
% Standard Axis and vice versa
RotGeoStd = [-1 0 0; 0 1 0; 0 0 -1];
%% Passing in variables from structured fields to local names, to make it less cumbersome
% NOTE: if a whole wing is declared as a control surface, then it must have a
% a single partition and no flaps
symmetry    = acraft_dsgn.symmetry;  
n_prt_wng   = acraft_dsgn.n_prt_wng; 
flapped     = acraft_dsgn.flapped;
flp_frac    = acraft_dsgn.flp_frac;
spn         = acraft_dsgn.spn;       
root_chrd   = acraft_dsgn.root_chrd;                              
tpr_rto     = acraft_dsgn.tpr_rto;   
swp_angle   = acraft_dsgn.swp_angle; 
dih_angle   = acraft_dsgn.dih_angle; 
xyz_000     = acraft_dsgn.xyz_000;  
twst_ang    = acraft_dsgn.twst_ang;  
wng_con_surf = [acraft_dsgn.wng_con_surf 
                0];  
acraft_dsgn.r_bo_CG = -acraft_dsgn.r_CG_bo;
% wng_con_surf = [0]; 
wng_con_surf_axis_rot    = acraft_dsgn.wng_con_surf_axis_rot;  %
% wng_con_surf_axis_rot    = [0 1 0]; 
airfoil     = acraft_dsgn.airfoil;
spn_div     = geo_disc.spn_div;  
chrd_div    = geo_disc.chrd_div; 
% % toc
%% GEOMETRIC DISCRETIZATION
% % tic
N_wng = size(n_prt_wng,1);

DOF_i = 0;
DOF_f = 0;
c = 0;    % Counter of wing control surfaces
count = 0; % Counter for AVL xyz section leading edge coordinates
flag = 0; % Flag indicating that a wing is a control surface
flap = 0; % Flap counter
Wng_con_surf_DOF_pointer = [];
wng_con_surf_axis_rot_aux = [];
wng_con_surf_hinge1 = [];
wng_con_surf_hinge2 = [];
Flap_DOF_pointer   = [];
flap_axis_rot      = [];
DOFstrips_i = 1;
for wng = 1:N_wng
    N_prt        = n_prt_wng(wng);  % # of pariitions in the wing = wng
    XYZL_tip_trpz = zeros(1,3);               % Resetting this value for each wing 
    
    Root_chrd =root_chrd(wng,1);
    
    if wng == wng_con_surf(c+1,1);
        c = c + 1;  % Counter of wing control surfaces
        flag = 1;   % Flag indicating a wing control surface
    end
    count = count + 1;
    AVL.xyz(count,:) = xyz_000(wng,:);
    AVL.chord(count) = Root_chrd;
    for prt = 1:N_prt
        % Figuring out how many trapezoids the partition is made of (i.e.
        % flapped of un-flapped)
        if flapped(wng,prt) == 1
            N_trpz = 2;
        else
            N_trpz = 1;
        end
        
        if symmetry(wng) == 1
            B         = 0.5*spn(wng,prt);
        end
        if symmetry(wng) == 0
            B         = spn(wng,prt);
        end
        
        if prt > 1
            Root_chrd     = Root_chrd*tpr_rto(wng,prt-1);
        end
        % Geometric paramteres of the partition
        TR            = tpr_rto(wng,prt);
        Swp_ang       = swp_angle(wng,prt);
        Dih_ang       = dih_angle(wng,prt);
        XYZ_000       = xyz_000(wng,:);
        Root_twst_ang = twst_ang(wng,prt);     % Angle in radians
        Tip_twst_ang  = twst_ang(wng,prt+1);   % Angle in radians
        Flap_frac     = flp_frac(wng,prt);
        Airfoil       = cell2mat(airfoil(wng,prt));
        if Airfoil == 0
            X = 0;
            Y = 0;
        else
            [X, Y] = arfoil_coord_parser(Airfoil);
        end        
        
        for trpz = 1:N_trpz  % Loops over the main (fixed) trapezoid and the flapped trapezoid if it exists 
            N     = spn_div(wng,prt,trpz);   % # of span-wise divisions in partition = prt of wing = wng
            M     = chrd_div(wng,prt,trpz);  % # of chord-wise divisions in partition = prt of wing = wng
            DOF_i = DOF_f + 1;     % Indicating the start and end of the DOFs corresponding to the trapezoid
            DOF_f = DOF_i + N*M - 1;
            DOFs_trpz  = DOF_i:DOF_f;
            
            if trpz == 1
                % Generate coordinates of the four points that define the
                % trapezoid
                XYZL_root_MCL  = 0.25*Root_chrd*[1 0 0];
                aux            = B/sqrt(1 + 1/(tan(pi*0.5-Swp_ang))^2);
                XYZL_tip_MCL   = XYZL_root_MCL + aux*[1/(tan(pi*0.5-Swp_ang)) cos(Dih_ang) sin(Dih_ang)];
                xyzL_pt1       = [0 0 0];   
                xyzL_pt2       = XYZL_tip_MCL  - 0.25*Root_chrd*TR*[1 0 0];
                xyzL_pt3       = xyzL_pt1      + Root_chrd*(1 - Flap_frac)*[1 0 0];
                xyzL_pt4       = xyzL_pt2      + Root_chrd*(1 - Flap_frac)*TR*[1 0 0];
                % Calculate XYZ for AVL
                count = count + 1;
                AVL.xyz(count,:) = xyzL_pt2 + XYZL_tip_trpz + xyz_000(wng,:);
                AVL.chord(count) = Root_chrd*TR;
                XYZL_tip_next_trpz = xyzL_pt2 + XYZL_tip_trpz; 
                if flag == 1
                    Wng_con_surf_DOF_pointer = [DOF_i DOF_f];  % Storing the DOF's corresponding to each Wing control surface
                    wng_con_surf_axis_rot_aux(c,:) = wng_con_surf_axis_rot(c,:);  % Auxiliary matrix to carry the axis of rotation of wing control surfaces (to avoid error when no wing control surfaces are declared) 
                    wng_con_surf_hinge1(c,:) = xyzL_pt1 + XYZL_tip_trpz + xyz_000(wng,:);
                    wng_con_surf_hinge2(c,:) = xyzL_pt2 + XYZL_tip_trpz + xyz_000(wng,:);
                end
            else
                % Generate coordinates of the four points that define the
                % trapezoid
                xyzL_pt1       = xyzL_pt3;
                xyzL_pt2       = xyzL_pt4;
                xyzL_pt3       = xyzL_pt1      + Root_chrd*Flap_frac*[1 0 0];
                xyzL_pt4       = xyzL_pt2      + Root_chrd*Flap_frac*TR*[1 0 0];
                
                flap = flap + 1;
                Flap_DOF_pointer(flap,:) = [DOF_i DOF_f];  % Storing the DOF's corresponding to each flap
                flap_axis_rot(flap,:)    = (xyzL_pt2 - xyzL_pt1)/norm(xyzL_pt2 - xyzL_pt1);  % Establishes the axis of rotation of the flap as that on the line from pt1 ro pt2
                flap_hinge1(flap,:) = xyzL_pt1 + XYZL_tip_trpz + xyz_000(wng,:);
                flap_hinge2(flap,:) = xyzL_pt2 + XYZL_tip_trpz + xyz_000(wng,:);
            end
            % PANELIZATION
            % Generating XYZ local to the trapezoid (using the apex of the partition as origin)
            xyzL = trapezoid_discretization(xyzL_pt1, xyzL_pt2, xyzL_pt3, xyzL_pt4, N, M, trpz, Flap_frac);
            
            % Translate the local XYZ to the global coordinate system. This includes fitting one partition on the tip of the previous one
            i = 0;
            for p = DOF_i:DOF_f
                i = i + 1;
                xyzG.crnr1(p,:) = xyzL.crnr1(i,:) + XYZL_tip_trpz + xyz_000(wng,:);
                xyzG.crnr2(p,:) = xyzL.crnr2(i,:) + XYZL_tip_trpz + xyz_000(wng,:);
                xyzG.crnr3(p,:) = xyzL.crnr3(i,:) + XYZL_tip_trpz + xyz_000(wng,:);
                xyzG.crnr4(p,:) = xyzL.crnr4(i,:) + XYZL_tip_trpz + xyz_000(wng,:);
                xyzG.HS1(p,:)   = xyzL.HS1(i,:)   + XYZL_tip_trpz + xyz_000(wng,:);
                xyzG.HS2(p,:)   = xyzL.HS2(i,:)   + XYZL_tip_trpz + xyz_000(wng,:);
                xyzG.col(p,:)   = xyzL.col(i,:)   + XYZL_tip_trpz + xyz_000(wng,:);
                xyzG.midpt(p,:) = xyzL.midpt(i,:) + XYZL_tip_trpz + xyz_000(wng,:);
                xyzG.Index(p) = wng;
                xyzG.chord(p) = xyzL.chord(i);
%                 
%                 % Read y,z coord
%                 x = xyzG.col(i,1);
%                 y = xyzG.col(i,2);
%                 z = xyzG.col(i,3);
%                 
%                 if (isempty(find(Strips(:,2) == y & Strips(:,3) == z & x <= Strips(:,4) & x >= Strips(:,1))))
%                     % create new strip with y,z
%                     % Strips = xyz of first col point, x of 1st col point + chord, xyz
%                     % of quarter chord point, chord length
%                     s = s+1;
%                     Strips(s,1:3) = [x y z];
%                     Strips(s,4) = x + xyzG.chord(i);
%                     Strips(s,5:7) = (xyzG.crnr1(i) + xyzG.crnr2(i))*0.5 + [1,0,0]*xyzG.chord(i)*0.25;
%                     Strips(s,8)  = xyzG.chord(i);
%                     Strip.chord = xyzG.chord(i);
%                     Strip.chord_dir =
%                     xyzG.Panel_strip(i) = s;
%                 else
%                     I = find(Strips(:,2) == y & Strips(:,3) == z & x <= Strips(:,4) & x >= Strips(:,1));% & x > Strips(:,1) & x < Strips(:,4));
%                     xyzG.Panel_strip(i) = I;
%                 end
            end
            
            % Calculating the normal vector at the collocation points for the trapezoid
            Normal_col(DOFs_trpz,:) = normal_col_discretization(Airfoil, X, Y, N, M, Flap_frac, trpz, Root_twst_ang, Tip_twst_ang, Dih_ang);
            
            % In case of Symmetry
            if symmetry(wng) == 1
                DOF_i = DOF_f + 1;     % Indicating the start and end of the DOFs corresponding to the trapezoid
                DOF_f = DOF_i + N*M - 1;
                DOFs_trpz  = DOF_i:DOF_f;
                
                Reflect_x  = [1  0 0
                              0 -1 0
                              0  0 1];

                xyzG.crnr1(DOFs_trpz,:) = xyzG.crnr2(DOF_i-N*M:DOF_f-N*M,:)*Reflect_x;
                xyzG.crnr2(DOFs_trpz,:) = xyzG.crnr1(DOF_i-N*M:DOF_f-N*M,:)*Reflect_x;
                xyzG.crnr3(DOFs_trpz,:) = xyzG.crnr4(DOF_i-N*M:DOF_f-N*M,:)*Reflect_x;
                xyzG.crnr4(DOFs_trpz,:) = xyzG.crnr3(DOF_i-N*M:DOF_f-N*M,:)*Reflect_x;
                xyzG.HS1(DOFs_trpz,:)   = xyzG.HS2(DOF_i-N*M:DOF_f-N*M,:)*Reflect_x;
                xyzG.HS2(DOFs_trpz,:)   = xyzG.HS1(DOF_i-N*M:DOF_f-N*M,:)*Reflect_x;
                xyzG.col(DOFs_trpz,:)   = xyzG.col(DOF_i-N*M:DOF_f-N*M,:)*Reflect_x;
                xyzG.midpt(DOFs_trpz,:) = xyzG.midpt(DOF_i-N*M:DOF_f-N*M,:)*Reflect_x;
                Normal_col(DOFs_trpz,:) = Normal_col(DOF_i-N*M:DOF_f-N*M,:)*Reflect_x;
                xyzG.Index(DOFs_trpz) = wng;
                xyzG.chord(DOFs_trpz) = xyzG.chord(DOF_i-N*M:DOF_f-N*M);
                if trpz == 1
                    if flag == 1
                    Wng_con_surf_DOF_pointer(c,2) = DOF_f;  % Extending the last DOF for this wing control surface to the last DOF in the symmetric wing 
                    end
                else
                    flap = flap + 1;
                    Flap_DOF_pointer(flap,:) = [DOF_i DOF_f];  % Storing the DOF's corresponding to each flap
                    flap_axis_rot(flap,:)    = flap_axis_rot(flap-1,:)*Reflect_x;  % Establishes the axis of rotation for the symmetric flap, by reflecting the axis of the previous one
                    flap_hinge1(flap,:) = flap_hinge1(flap-1,:)*Reflect_x;
                    flap_hinge2(flap,:) = flap_hinge2(flap-1,:)*Reflect_x;
                end
            end
            
        end
            XYZL_tip_trpz = XYZL_tip_next_trpz;  % Updating the XYZ of the tip of the previos partition, so the next one can be fitted at this point
    end
    % Go over all the DOFs belonging to a given wing, and assign each one
    % of them to a corresonping strip
    DOFstrips_f = size(xyzG.col,1);
    DOFs = DOFstrips_i:DOFstrips_f;
%    xyzG.Strip_index(DOFs) = assign_strip(xyzG.col,DOFs);
    xyzG.wng(DOFs) = wng;
    DOFstrips_i = DOFstrips_f+1;
    flag = 0;  % Reset flag indicating a wing control surface
end
n_panels = length(xyzG.crnr1);   % Calculating the number of panels in the problem

% % toc
%% Assemblying the DOF pointers for the control surfaces (Wings + flaps) and axis of rotation
con_surf_DOF_pointer = [Wng_con_surf_DOF_pointer;
                        Flap_DOF_pointer];     
con_surf_axis_rot    = [wng_con_surf_axis_rot_aux;
                        flap_axis_rot];              
con_surf_hinge1 = [wng_con_surf_hinge1
                    flap_hinge1]; 
con_surf_hinge2 = [wng_con_surf_hinge2
                    flap_hinge2];
                    
N_con_surf           = size(con_surf_DOF_pointer,1);
% Passing these new variables into acrft.dsgn structure
acraft_dsgn.N_con_surf           = N_con_surf;
acraft_dsgn.con_surf_DOF_pointer = con_surf_DOF_pointer;
acraft_dsgn.con_surf_axis_rot    = con_surf_axis_rot;
xyzG.Normal_col = Normal_col;
xyzG.con_surf_hinge1 = con_surf_hinge1;
xyzG.con_surf_hinge2 = con_surf_hinge2;
acraft_dsgn.xyzG = xyzG;

% Rotate xyzG coordinates, normal collocation vectors and hinge vectors
xyzG_Std.crnr1      = xyzG.crnr1*RotGeoStd;
xyzG_Std.crnr2      = xyzG.crnr2*RotGeoStd;
xyzG_Std.crnr3      = xyzG.crnr3*RotGeoStd;
xyzG_Std.crnr4      = xyzG.crnr4*RotGeoStd;
xyzG_Std.HS1        = xyzG.HS1*RotGeoStd;
xyzG_Std.HS2        = xyzG.HS2*RotGeoStd;
xyzG_Std.col        = xyzG.col*RotGeoStd;
xyzG_Std.midpt      = xyzG.midpt*RotGeoStd;
xyzG_Std.Normal_col = xyzG.Normal_col*RotGeoStd;
xyzG_Std.con_surf_hinge1 = con_surf_hinge1*RotGeoStd;
xyzG_Std.con_surf_hinge2 = con_surf_hinge2*RotGeoStd;

if isempty(acraft_dsgn.con_surf_axis_rot)
else
acraft_dsgn.con_surf_axis_rot_Std   = acraft_dsgn.con_surf_axis_rot*RotGeoStd;
end

% %% Constructing the AIC matrix and the Force-Moment matrices
% %% QAUDAIR 1.4 version (all coordinates and vectors in Standar frame
% x_infinity = -40*abs(min(xyzG.crnr1(:,1))-max(xyzG.crnr3(:,1)));  % Finding the maximum characteris% tic distance in the x-direction, to get what we will consider as "x_infinty" or a point really far downstrea,
% 

% % tic
% size_X = 6 + 9;
% %% Construct L_tilde matrix for B.C. equation
% % Construct L_i_o, L_i_c and L_i_s
% L_o = zeros(n_panels, size_X);
% L_tilde1 = [];
% DOFs_no_con_surf = 1:n_panels;
% for i = 1:N_con_surf
%     DOFs = con_surf_DOF_pointer(i,1):con_surf_DOF_pointer(i,2);  % DOFs corresponding to i-th control surface
%     DOFs_no_con_surf = setdiff(DOFs_no_con_surf,DOFs);
%     L_i_o = zeros(n_panels, size_X); L_i_c = zeros(n_panels, size_X); L_i_s = zeros(n_panels, size_X); 
%     for k = DOFs
%         L1   = acraft_dsgn.con_surf_axis_rot_Std(i,1); L2   = acraft_dsgn.con_surf_axis_rot_Std(i,2); L3   = acraft_dsgn.con_surf_axis_rot_Std(i,3);
%         R_o = [L1^2    L1*L2 L3*L1
%             L1*L2   L2^2  L2*L3
%             L3*L1   L2*L3 L3^2];
%         
%         R_c = [1-L1^2 -L1*L2 -L3*L1
%               -L1*L2   1-L2^2 -L2*L3
%               -L3*L1  -L2*L3 1-L3^2];
%         
%         R_s = [0    L3 -L2
%               -L3   0   L1
%                L2  -L1  0];
%         
%         r_col = xyzG_Std.col(k,:);
%         
%         r_cross = [0         -r_col(3)   r_col(2)
%                    r_col(3)   0         -r_col(1)
%                   -r_col(2)   r_col(1)   0       ];
%         
%         r_tilde = [r_col        zeros(1,3)  zeros(1,3)
%                    zeros(1,3)   r_col       zeros(1,3)
%                    zeros(1,3)   zeros(1,3)  r_col   ];
%         
%         d_k = [eye(3) -r_cross -r_tilde];
%         S_n_k = xyzG_Std.Normal_col(k,:);
%         
%         L_i_o(k,:) = S_n_k*R_o*d_k;
%         L_i_c(k,:) = S_n_k*R_c*d_k;
%         L_i_s(k,:) = S_n_k*R_s*d_k;
%     end
%     L_tilde1 = [L_tilde1 L_i_c L_i_s];
%   
%     L_o = L_o + L_i_o;
% %    R11 = c_th + L1^2*(1 - c_th); R22 = c_th + L2^2*(1 - c_th); R33 = c_th + L3^2*(1 - c_th); R12 =  L3*s_th + L1*L2*(1 - c_th); R21 = -L3*s_th + L1*L2*(1 - c_th); R13 = -L2*s_th + L3*L1*(1 - c_th); R31 =  L2*s_th + L3*L1*(1 - c_th); R23 =  L1*s_th + L2*L3*(1 - c_th);  R32 = -L1*s_th + L2*L3*(1 - c_th);
% end
% 
% % Construct L^o
% for k = DOFs_no_con_surf
%     r_col = xyzG_Std.col(k,:);
%     
%     r_cross = [0         -r_col(3)   r_col(2)
%                r_col(3)   0         -r_col(1)
%               -r_col(2)   r_col(1)   0       ];
%     
%     r_tilde = [r_col        zeros(1,3)  zeros(1,3)
%                zeros(1,3)   r_col       zeros(1,3)
%                zeros(1,3)   zeros(1,3)  r_col   ];
%     d_k = [eye(3) -r_cross -r_tilde];
%     S_n_k = xyzG_Std.Normal_col(k,:);
%     L_o(k,:) = S_n_k*d_k;
% end
% L_tilde = [L_o L_tilde1];
% 
% %% Construct AIC for B.C. equation, G and H matrices for Force and Moment equations
% AIC_new = zeros(n_panels,n_panels);
% for k = 1:n_panels       % Sweeping over the collocation points
%     % For congruence with paper derivations, we add index j = k
%     j = k;
%     % Finding the midpoint on the quarter chord of the panel (goesright through the bounded vortex filament)
%     xpbound = (xyzG_Std.HS1(k,1) + xyzG_Std.HS2(k,1))*0.5;
%     ypbound = (xyzG_Std.HS1(k,2) + xyzG_Std.HS2(k,2))*0.5;
%     zpbound = (xyzG_Std.HS1(k,3) + xyzG_Std.HS2(k,3))*0.5;
%     
%     %xyzG_Std.bound(k,:) = [xpbound ypbound zpbound];
%     r_bound = [xpbound ypbound zpbound]; 
% 
%     % Finding the unit vector of the bounded vortex filament
%     U = xyzG_Std.HS2(k,:) - xyzG_Std.HS1(k,:);
%     % Matrices for cross-products
%     U_cross = -1*[0    -U(3)  U(2)
%                   U(3)  0    -U(1)
%                  -U(2)  U(1)  0];
%     r_cross = [0        -zpbound   ypbound
%                zpbound   0        -xpbound
%               -ypbound   xpbound   0];
%           
%     r_tilde = [r_bound      zeros(1,3)  zeros(1,3)
%                zeros(1,3)   r_bound     zeros(1,3)
%                zeros(1,3)   zeros(1,3)  r_bound   ];
%           
%     kappa(:,:,j) = [-eye(3) r_cross r_tilde zeros(3,2*size_X*N_con_surf)];
%     
%     Ucross_kappa = U_cross*kappa(:,:,j);
%     rcross_Ucross_kappa = r_cross*U_cross*kappa(:,:,j);
%     
%     Gx_Xtilde(k,:) = Ucross_kappa(1,:);
%     Gy_Xtilde(k,:) = Ucross_kappa(2,:);
%     Gz_Xtilde(k,:) = Ucross_kappa(3,:);
%     
%     rcross_Ucross_kappa = r_cross*U_cross*kappa(:,:,j);
%     Hx_Xtilde(k,:) = rcross_Ucross_kappa(1,:);
%     Hy_Xtilde(k,:) = rcross_Ucross_kappa(2,:);
%     Hz_Xtilde(k,:) = rcross_Ucross_kappa(3,:);
%     
%     xpcol= xyzG_Std.col(k,1);
%     ypcol= xyzG_Std.col(k,2);
%     zpcol= xyzG_Std.col(k,3);
%     for i = 1:n_panels   % Sweeping over each of the HS vor% tices
%         %% AIC matrix
%         x1 = xyzG_Std.HS1(i,1);
%         y1 = xyzG_Std.HS1(i,2);
%         z1 = xyzG_Std.HS1(i,3);
%         
%         x2 = xyzG_Std.HS2(i,1);
%         y2 = xyzG_Std.HS2(i,2);
%         z2 = xyzG_Std.HS2(i,3);
%         
%         [u, v, w] = HS_vortex(xpcol, ypcol, zpcol, x1, y1, z1, x2, y2, z2, 1, x_infinity,xyzG.Index(k),xyzG.Index(i), xyzG.chord(i));
%         
%         AIC_new(k,i) = u*xyzG_Std.Normal_col(k,1) + v*xyzG_Std.Normal_col(k,2)+ w*xyzG_Std.Normal_col(k,3); 
%         
%         %% Force-Moment matrices
%         [u, v, w] = HS_vortex(xpbound, ypbound, zpbound, x1, y1, z1, x2, y2, z2, 1, x_infinity,xyzG.Index(k),xyzG.Index(i), xyzG.chord(i));
%         
%         %% Matrices for velocity computation
%         % Induced velocity part
%         h_ji = [u v w]';
%         h_j(:,i,j) = h_ji;
%     end
%     Ucross_h_j = U_cross*h_j(:,:,j);  
%     Gx_gam_new(k,:) = Ucross_h_j(1,:);
%     Gy_gam_new(k,:) = Ucross_h_j(2,:);
%     Gz_gam_new(k,:) = Ucross_h_j(3,:);
%     
%     rcross_Ucross_h_j = r_cross*U_cross*h_j(:,:,j);
%     Hx_gam_new(k,:) = rcross_Ucross_h_j(1,:);
%     Hy_gam_new(k,:) = rcross_Ucross_h_j(2,:);
%     Hz_gam_new(k,:) = rcross_Ucross_h_j(3,:);    
% end
% 
% J_new = AIC_new\L_tilde;
% 
% AeroMatrices.Px_new = J_new'*(Gx_Xtilde + Gx_gam_new*J_new);
% AeroMatrices.Py_new = J_new'*(Gy_Xtilde + Gy_gam_new*J_new);
% AeroMatrices.Pz_new = J_new'*(Gz_Xtilde + Gz_gam_new*J_new);
% 
% rx = acraft_dsgn.r_bo_CG(1);
% ry = acraft_dsgn.r_bo_CG(2);
% rz = acraft_dsgn.r_bo_CG(3);
% 
% AeroMatrices.Qx_new = J_new'*(Hx_Xtilde + Hx_gam_new*J_new);
% AeroMatrices.Qy_new = J_new'*(Hy_Xtilde + Hy_gam_new*J_new);
% AeroMatrices.Qz_new = J_new'*(Hz_Xtilde + Hz_gam_new*J_new);
% 
% AeroMatrices.Qx_new = AeroMatrices.Qx_new + ry*AeroMatrices.Pz_new - rz*AeroMatrices.Py_new; 
% AeroMatrices.Qy_new = AeroMatrices.Qy_new + rz*AeroMatrices.Px_new - rx*AeroMatrices.Pz_new; 
% AeroMatrices.Qz_new = AeroMatrices.Qz_new + rx*AeroMatrices.Py_new - ry*AeroMatrices.Px_new;  
% 
% AeroMatrices.Px_new = (AeroMatrices.Px_new' + AeroMatrices.Px_new)/2;
% AeroMatrices.Py_new = (AeroMatrices.Py_new' + AeroMatrices.Py_new)/2;
% AeroMatrices.Pz_new = (AeroMatrices.Pz_new' + AeroMatrices.Pz_new)/2;
% AeroMatrices.Qx_new = (AeroMatrices.Qx_new' + AeroMatrices.Qx_new)/2;
% AeroMatrices.Qy_new = (AeroMatrices.Qy_new' + AeroMatrices.Qy_new)/2;
% AeroMatrices.Qz_new = (AeroMatrices.Qz_new' + AeroMatrices.Qz_new)/2;
% 
% % toc







%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% QUADAIR 1.2 version
 % % tic
 x_infinity = 40*abs(min(xyzG.crnr1(:,1))-max(xyzG.crnr3(:,1)));  % Finding the maximum characteris% tic distance in the x-direction, 
                                                              % to get what we will consider as "x_infinty" or a point really far downstrea,
 % tic
AIC = zeros(n_panels,n_panels);
for k = 1:n_panels       % Sweeping over the collocation points
    % Finding the midpoint on the quarter chord of the panel (goes
    % right through the bounded vortex filament)
    xpbound = (xyzG.HS1(k,1) + xyzG.HS2(k,1))*0.5;
    ypbound = (xyzG.HS1(k,2) + xyzG.HS2(k,2))*0.5;
    zpbound = (xyzG.HS1(k,3) + xyzG.HS2(k,3))*0.5;
    
    xyzG.bound(k,:) = [xpbound ypbound zpbound];

    % Finding the unit vector of the bounded vortex filament
    U = xyzG.HS2(k,:) - xyzG.HS1(k,:);
    % Matrices for cross-products
    U_cross = [0    -U(3)  U(2)
               U(3)  0    -U(1)
              -U(2)  U(1)  0];
    r_cross = [0        -zpbound   ypbound
               zpbound   0        -xpbound
              -ypbound   xpbound   0];
    % Matrices for force and moment calculation
    GxV(k,:) = U_cross(1:3,1)';
    GyV(k,:) = U_cross(1:3,2)';
    GzV(k,:) = U_cross(1:3,3)';
    
    %GxW = r cross -(U cross x_hat)
    GxW(k,:) = (r_cross*U_cross*[1 0 0]')';
    GyW(k,:) = (r_cross*U_cross*[0 1 0]')';
    GzW(k,:) = (r_cross*U_cross*[0 0 1]')';
    
    % HxV = u cross -(r cross x_hat)
    HxV(k,:) = U_cross*-r_cross*[1 0 0]';
    HyV(k,:) = U_cross*-r_cross*[0 1 0]';
    HzV(k,:) = U_cross*-r_cross*[0 0 1]';
    
    % HxW = r cross ( U cross (-r cross x_hat)))
    HxW(k,:) = r_cross*U_cross*-r_cross*[1 0 0]';
    HyW(k,:) = r_cross*U_cross*-r_cross*[0 1 0]';
    HzW(k,:) = r_cross*U_cross*-r_cross*[0 0 1]';
    
    xpcol= xyzG.col(k,1);
    ypcol= xyzG.col(k,2);
    zpcol= xyzG.col(k,3);
    for i = 1:n_panels   % Sweeping over each of the HS vor% tices
        %% AIC matrix
        x1 = xyzG.HS1(i,1);
        y1 = xyzG.HS1(i,2);
        z1 = xyzG.HS1(i,3);
        
        x2 = xyzG.HS2(i,1);
        y2 = xyzG.HS2(i,2);
        z2 = xyzG.HS2(i,3);
        
        [u, v, w] = HS_vortex(xpcol, ypcol, zpcol, x1, y1, z1, x2, y2, z2, 1, x_infinity,xyzG.Index(k),xyzG.Index(i), xyzG.chord(i));
        
        AIC(k,i) = u*Normal_col(k,1) + v*Normal_col(k,2)+ w*Normal_col(k,3); 
        
        %% Force-Moment matrices
        [u, v, w] = HS_vortex(xpbound, ypbound, zpbound, x1, y1, z1, x2, y2, z2, 1, x_infinity,xyzG.Index(k),xyzG.Index(i), xyzG.chord(i));
        
        h_cross = [0  -w   v
                   w   0  -u
                  -v   u   0];
        % Gx_gam = (hij cross U) dot x_hat
        Gx_gam(k,i) = [1 0 0]*h_cross*U';
        Gy_gam(k,i) = [0 1 0]*h_cross*U';
        Gz_gam(k,i) = [0 0 1]*h_cross*U';
        
        % Hx_gam = (r cross (h cross U)) dot x_hat
        Hx_gam(k,i) = [1 0 0]*r_cross*h_cross*U';
        Hy_gam(k,i) = [0 1 0]*r_cross*h_cross*U';
        Hz_gam(k,i) = [0 0 1]*r_cross*h_cross*U'; 
        
        %% Matrices for velocity computation
        % For congruence with paper derivations, we add index j = k
        j = k;
%         % Rigid body translation and rotation part
        AA(:,:,j) = [eye(3) r_cross zeros(3, 12*N_con_surf)];
%         % Induced velocity part
        h_ji = [u v w]';
        HH(:,i,j) = h_ji;
        
        
        
    
    end
end
% toc


%% Constructing the RHS's and solving them
% % tic
Normal_col_con_surf = zeros(n_panels,3,N_con_surf);
RxNx = zeros(n_panels,N_con_surf);
RxNy = zeros(n_panels,N_con_surf);
RxNz = zeros(n_panels,N_con_surf);

RyNx = zeros(n_panels,N_con_surf);
RyNy = zeros(n_panels,N_con_surf);
RyNz = zeros(n_panels,N_con_surf);

RzNx = zeros(n_panels,N_con_surf);
RzNy = zeros(n_panels,N_con_surf);
RzNz = zeros(n_panels,N_con_surf);

Normal_col_Fixed = Normal_col;
RCNx_Fixed = xyzG.col(:,2).*Normal_col(:,3) - xyzG.col(:,3).*Normal_col(:,2);
RCNy_Fixed = xyzG.col(:,3).*Normal_col(:,1) - xyzG.col(:,1).*Normal_col(:,3);
RCNz_Fixed = xyzG.col(:,1).*Normal_col(:,2) - xyzG.col(:,2).*Normal_col(:,1);
J = zeros(n_panels,6);
for i = 1:N_con_surf
    %% Setting up (34) regarding rotation matrix decomposition
    L1   = acraft_dsgn.con_surf_axis_rot(i,1); L2   = acraft_dsgn.con_surf_axis_rot(i,2); L3   = acraft_dsgn.con_surf_axis_rot(i,3);
    % Independent matrix
    A = [L1^2  L1*L2 L3*L1
                L1*L2 L2^2  L2*L3
                L3*L1 L2*L3 L3^2];
    B = [1-L1^2 -L1*L2 -L3*L1
                -L1*L2 1-L2^2 -L2*L3
                -L3*L1 -L2*L3 1-L3^2];
    C = [0   L3 -L2
                -L3 0   L1
                L2  -L1 0];    

    DOFs = con_surf_DOF_pointer(i,1):con_surf_DOF_pointer(i,2);  % DOFs corresponding to i-th control surface
    
    % Qinf terms
    Normal_col_con_surf(DOFs,:,i) = Normal_col(DOFs,:);
  
    % Angular velocity terms
    RxNx(DOFs,i) = xyzG.col(DOFs,1).*Normal_col(DOFs,1);
    RxNy(DOFs,i) = xyzG.col(DOFs,1).*Normal_col(DOFs,2);
    RxNz(DOFs,i) = xyzG.col(DOFs,1).*Normal_col(DOFs,3);
    
    RyNx(DOFs,i) = xyzG.col(DOFs,2).*Normal_col(DOFs,1);
    RyNy(DOFs,i) = xyzG.col(DOFs,2).*Normal_col(DOFs,2);
    RyNz(DOFs,i) = xyzG.col(DOFs,2).*Normal_col(DOFs,3);
    
    RzNx(DOFs,i) = xyzG.col(DOFs,3).*Normal_col(DOFs,1);
    RzNy(DOFs,i) = xyzG.col(DOFs,3).*Normal_col(DOFs,2);
    RzNz(DOFs,i) = xyzG.col(DOFs,3).*Normal_col(DOFs,3);
    
    % Putting zeros in the rows that correspond to control surfaces DOFs
    Normal_col_Fixed(DOFs,:) = 0;
    RCNx_Fixed(DOFs) = 0;
    RCNy_Fixed(DOFs) = 0;
    RCNz_Fixed(DOFs) = 0;
    
    % Making A_emf, B_emf, C_emf Eq. (35) & (36), without includeing n_fix
    % in A_emf

    A_emf = -AIC\(Normal_col_con_surf(:,:,i)*A);
    B_emf = -AIC\(Normal_col_con_surf(:,:,i)*B);
    C_emf = -AIC\(Normal_col_con_surf(:,:,i)*C);

    % Decomposing RCN for moving surface j into independent, cosine and
    % sine parts
    j = i;
    c_th = 0;
    s_th = 0;
    R11 = c_th + L1^2*(1 - c_th);      R22 = c_th + L2^2*(1 - c_th);      R33 = c_th + L3^2*(1 - c_th); R12 =  L3*s_th + L1*L2*(1 - c_th); R21 = -L3*s_th + L1*L2*(1 - c_th); R13 = -L2*s_th + L3*L1*(1 - c_th); R31 =  L2*s_th + L3*L1*(1 - c_th); R23 =  L1*s_th + L2*L3*(1 - c_th); R32 = -L1*s_th + L2*L3*(1 - c_th);
    RCN_indep(:,1,j) = R13*RyNx(:,j) + R23*RyNy(:,j) + R33*RyNz(:,j) - R12*RzNx(:,j) - R22*RzNy(:,j) - R32*RzNz(:,j);
    RCN_indep(:,2,j) = R11*RzNx(:,j) + R21*RzNy(:,j) + R31*RzNz(:,j) - R13*RxNx(:,j) - R23*RxNy(:,j) - R33*RxNz(:,j);
    RCN_indep(:,3,j) = R12*RxNx(:,j) + R22*RxNy(:,j) + R32*RxNz(:,j) - R11*RyNx(:,j) - R21*RyNy(:,j) - R31*RyNz(:,j);

    c_th = 1;
    s_th = 0;
    R11 = c_th + L1^2*(1 - c_th);      R22 = c_th + L2^2*(1 - c_th);      R33 = c_th + L3^2*(1 - c_th); R12 =  L3*s_th + L1*L2*(1 - c_th); R21 = -L3*s_th + L1*L2*(1 - c_th); R13 = -L2*s_th + L3*L1*(1 - c_th); R31 =  L2*s_th + L3*L1*(1 - c_th); R23 =  L1*s_th + L2*L3*(1 - c_th); R32 = -L1*s_th + L2*L3*(1 - c_th);
    RCN_c_th(:,1,j) = R13*RyNx(:,j) + R23*RyNy(:,j) + R33*RyNz(:,j) - R12*RzNx(:,j) - R22*RzNy(:,j) - R32*RzNz(:,j);
    RCN_c_th(:,2,j) = R11*RzNx(:,j) + R21*RzNy(:,j) + R31*RzNz(:,j) - R13*RxNx(:,j) - R23*RxNy(:,j) - R33*RxNz(:,j);
    RCN_c_th(:,3,j) = R12*RxNx(:,j) + R22*RxNy(:,j) + R32*RxNz(:,j) - R11*RyNx(:,j) - R21*RyNy(:,j) - R31*RyNz(:,j);
    RCN_c_th(:,:,j) = RCN_c_th(:,:,j) - RCN_indep(:,:,j);
    
    c_th = 0;
    s_th = 1;
    R11 = c_th + L1^2*(1 - c_th);      R22 = c_th + L2^2*(1 - c_th);      R33 = c_th + L3^2*(1 - c_th); R12 =  L3*s_th + L1*L2*(1 - c_th); R21 = -L3*s_th + L1*L2*(1 - c_th); R13 = -L2*s_th + L3*L1*(1 - c_th); R31 =  L2*s_th + L3*L1*(1 - c_th); R23 =  L1*s_th + L2*L3*(1 - c_th); R32 = -L1*s_th + L2*L3*(1 - c_th);
    RCN_s_th(:,1,j) = R13*RyNx(:,j) + R23*RyNy(:,j) + R33*RyNz(:,j) - R12*RzNx(:,j) - R22*RzNy(:,j) - R32*RzNz(:,j);
    RCN_s_th(:,2,j) = R11*RzNx(:,j) + R21*RzNy(:,j) + R31*RzNz(:,j) - R13*RxNx(:,j) - R23*RxNy(:,j) - R33*RxNz(:,j);
    RCN_s_th(:,3,j) = R12*RxNx(:,j) + R22*RxNy(:,j) + R32*RxNz(:,j) - R11*RyNx(:,j) - R21*RyNy(:,j) - R31*RyNz(:,j);
    RCN_s_th(:,:,j) = RCN_s_th(:,:,j) - RCN_indep(:,:,j);
    
    D_emf = AIC\RCN_indep(:,:,j);
    E_emf = AIC\RCN_c_th(:,:,j);
    F_emf = AIC\RCN_s_th(:,:,j);
    
    
%     6+12*(j-1)+1:6+12*(j-1)+12
    ii = 6+12*(j-1)+1;            % Sets the pointer for J
    J(:,1:3) = J(:,1:3) + A_emf;
    J(:,4:6) = J(:,4:6) + D_emf;
    J(:,ii:ii+11) = [B_emf E_emf C_emf F_emf];
    
end
sol_N_Fixed    = AIC\Normal_col_Fixed;
sol_RCNx_Fixed = AIC\RCNx_Fixed;
sol_RCNy_Fixed = AIC\RCNy_Fixed;
sol_RCNz_Fixed = AIC\RCNz_Fixed;

% Adding components from fixed normals to J
J(:,1:3) = J(:,1:3) - sol_N_Fixed;
J(:,4:6) = J(:,4:6) + [sol_RCNx_Fixed sol_RCNy_Fixed sol_RCNz_Fixed];

% Making P and Q matrcices
% The minus signs added to Px, Pz, Qx and Qz are to switch the computed
% Forces and Moments so that they are w.r.t. the Standard Frame, instead of the Geometric
% Frame.  This change was addded when QuadAir 1.0 was created from AircraftPre_computer5,
% to facilitate flight simulation code.  There is also an adjustemnt of the
% moments so that they are with respect to CG instead of Bo.

Px = -J'*( [GxV -GxW zeros(n_panels,12*N_con_surf)] + Gx_gam*J ); 
Py = J'*( [GyV -GyW zeros(n_panels,12*N_con_surf)] + Gy_gam*J ); 
Pz = -J'*( [GzV -GzW zeros(n_panels,12*N_con_surf)] + Gz_gam*J ); 

Qx = -J'*( [HxV -HxW zeros(n_panels,12*N_con_surf)] + Hx_gam*J ); 
Qy = J'*( [HyV -HyW zeros(n_panels,12*N_con_surf)] + Hy_gam*J ); 
Qz = -J'*( [HzV -HzW zeros(n_panels,12*N_con_surf)] + Hz_gam*J ); 

AeroMatrices.J = J;



%% Matrices for analy% tic gradient
n = 6 + 12*acraft_dsgn.N_con_surf; 
N = acraft_dsgn.N_con_surf;
M = max(abs(acraft_dsgn.con_surf_group));

A = zeros(n,6); 
A(1:6,:) = eye(6);

C = zeros(n,6,N); 
S = zeros(n,6,N);

K = zeros(N,M);

% We mutiply by D matrix to adjust for the fact that internally we use Qinf and write its
% components in the Geometric Frame. D = Rotation from Geo to Std, with a
% negative sign added to the velocity part.
D = [1  0  0   0  0  0;
     0 -1  0   0  0  0 
     0  0  1   0  0  0 
     0  0  0  -1  0  0 
     0  0  0   0  1  0
     0  0  0   0  0 -1];
% We multiply by L to adjust for the fact that X is internally taken with respect to Bo, and the user
% will want it with respect to CG.
rx = acraft_dsgn.r_bo_CG(1);
ry = acraft_dsgn.r_bo_CG(2);
rz = acraft_dsgn.r_bo_CG(3);
r_cross = [0   -rz   ry
           rz   0   -rx
          -ry   rx   0];

L = [eye(3) -r_cross
    zeros(3,3) eye(3)];
for j = 1:acraft_dsgn.N_con_surf
    % Matrices to compute Xtilde
    ii = 6+12*(j-1)+1;
    C(ii:ii+5,:,j) = eye(6);
    
    S(ii+6:ii+11,:,j) = eye(6);
    
    C(:,:,j) = C(:,:,j)*D*L;
    S(:,:,j) = S(:,:,j)*D*L;
    
    % Matrix to compute Theta_c
    group_nbr = abs(acraft_dsgn.con_surf_group(j));
    
    K(j,group_nbr) = sign(acraft_dsgn.con_surf_group(j));
    
end

A = A*D*L;

AeroMatrices.A = A;
AeroMatrices.C = C;
AeroMatrices.S = S;
AeroMatrices.K = K;
AeroMatrices.XFS2XQuadAir = D*L;
AeroMatrices.XFS2XStdVBo = L;
AeroMatrices.XStdVBo2XGeoQinfBo = D;
AeroMatrices.D = D;
AeroMatrices.L = L;

DD = [];
for i = 1:(size(Px,2)/6)
    DD = blkdiag(DD,AeroMatrices.XFS2XQuadAir);
end



% Adjusting Qx, Qy and Qz for the fact that the user is interested in
% moments about CG, instead of Bo

Qx = Qx + ry*Pz - rz*Py; 
Qy = Qy + rz*Px - rx*Pz; 
Qz = Qz + rx*Py - ry*Px;  

% Adjusting P and Q matrices to allow X to be defined directly in CG_V_Std variables
Px = DD'*Px*DD;
Py = DD'*Py*DD;
Pz = DD'*Pz*DD;
Qx = DD'*Qx*DD;
Qy = DD'*Qy*DD;
Qz = DD'*Qz*DD;

AeroMatrices.Px = (Px' + Px)/2;
AeroMatrices.Py = (Py' + Py)/2;
AeroMatrices.Pz = (Pz' + Pz)/2;
AeroMatrices.Qx = (Qx' + Qx)/2;
AeroMatrices.Qy = (Qy' + Qy)/2;
AeroMatrices.Qz = (Qz' + Qz)/2;






% Reduce P & Q matrices due to couples control surfaces
N_control_groups = max(acraft_dsgn.con_surf_group);
if ~isempty(acraft_dsgn.con_surf_group)
    R = zeros(6+acraft_dsgn.N_con_surf*6*2, 6 + N_control_groups*6*2);
    R(1:6,1:6) = eye(6);
    
    for j = 1:acraft_dsgn.N_con_surf
        jj = 6 + 12*(j-1) + 1;
        i = abs(acraft_dsgn.con_surf_group(j));
        ii = 6 + 12*(i-1) + 1;
        if (sign(acraft_dsgn.con_surf_group(j)) == 1)
            
            R(jj:jj+11, ii:ii+11) = [eye(6)       zeros(6,6)
                                    zeros(6,6)    eye(6)];
        else
            R(jj:jj+11, ii:ii+11) = [eye(6)       zeros(6,6)
                                    zeros(6,6)    -eye(6)];
        end
    end
    
    AeroMatrices.Px_red = R'*Px*R;
    AeroMatrices.Py_red = R'*Py*R;
    AeroMatrices.Pz_red = R'*Pz*R;
    AeroMatrices.Qx_red = R'*Qx*R;
    AeroMatrices.Qy_red = R'*Qy*R;
    AeroMatrices.Qz_red = R'*Qz*R;
    
    AeroMatrices.J_red  = J*R;
    AeroMatrices.ReduceXtilde = R;
else
    if isempty(N_control_groups)
        R = eye(6);
    end
end

% %%%%%%%%%
% % Matrix reduction for QAIR1.4
% 
% N_control_groups = max(acraft_dsgn.con_surf_group);
% if ~isempty(acraft_dsgn.con_surf_group)
%     R_new = zeros(6+acraft_dsgn.N_con_surf*size_X*2, size_X + N_control_groups*size_X*2);
%     R_new(1:size_X,1:size_X) = eye(size_X);
%     
%     for j = 1:acraft_dsgn.N_con_surf
%         jj = size_X + 2*size_X*(j-1) + 1;
%         i = abs(acraft_dsgn.con_surf_group(j));
%         ii = size_X + 2*size_X*(i-1) + 1;
%         if (sign(acraft_dsgn.con_surf_group(j)) == 1)
%             
%             R_new(jj:jj+2*size_X-1, ii:ii+2*size_X-1) = [eye(size_X)       zeros(size_X,size_X)
%                                     zeros(size_X,size_X)    eye(size_X)];
%         else
%             R_new(jj:jj+2*size_X-1, ii:ii+2*size_X-1) = [eye(size_X)       zeros(size_X,size_X)
%                                     zeros(size_X,size_X)    -eye(size_X)];
%         end
%     end
%     
%     AeroMatrices.Px_red_new = R_new'*AeroMatrices.Px_new*R_new;
%     AeroMatrices.Py_red_new = R_new'*AeroMatrices.Py_new*R_new;
%     AeroMatrices.Pz_red_new = R_new'*AeroMatrices.Pz_new*R_new;
%     AeroMatrices.Qx_red_new = R_new'*AeroMatrices.Qx_new*R_new;
%     AeroMatrices.Qy_red_new = R_new'*AeroMatrices.Qy_new*R_new;
%     AeroMatrices.Qz_red_new = R_new'*AeroMatrices.Qz_new*R_new;
%     
%     AeroMatrices.J_red_new  = J_new*R_new;
%     AeroMatrices.ReduceXtilde_new = R_new;
% else
%     if isempty(N_control_groups)
%         R_new = eye(size_X);
%     end
% end
% %%%%%%%%%

%% Velocity calculation at bound vor% tices
VelMat = zeros(size(AA(:,:,1)));
for j = 1:n_panels
    VelMat(:,:,j) = RotGeoStd*(AA(:,:,j) + HH(:,:,j)*J);
    VelMat_Red(:,:,j) = VelMat(:,:,j)*R;
end

AeroMatrices.VelMat = VelMat;
AeroMatrices.VelMat_Red = VelMat_Red;

% %%%%%%%
% %% QAIR 1.4 Construct Velocity matrices K_j
% for j = 1:n_panels
%     K_j(:,:,j) = kappa(:,:,j) + h_j(:,:,j)*J_new;
%     K_j_red(:,:,j) = K_j(:,:,j)*R_new;
% end
% AeroMatrices.VelMat_new = K_j;
% AeroMatrices.VelMat_Red_new = K_j_red;
% 
% %%%%%%%%


%% Stuff for stalled model
% tic
N = size(xyzG.col,1);
Strips = zeros(1,4);
s = 0;
% wing_strip is used to map strips to wings, by listing the start and stop
% of strip numbers for each wing
wing_strip = sum(spn_div(:,:,1),2).*(symmetry+1);
for i = 2:size(wing_strip,1)
    wing_strip(i) = wing_strip(i) + wing_strip(i-1); 
end
aux0 = spn_div(:,:,1).*repmat(symmetry+1,1,size(spn,2));
aux1 = reshape(aux0',size(spn_div(:,:,1),1)*size(spn_div(:,:,1),2),1); 
aux2 = cumsum(aux1);
part_strip = reshape(aux2, size(spn_div(:,:,1),2), size(spn_div(:,:,1),1))';
% part_strip is used to maps strips to partitions

for i = 1:N
    
    % Read y,z coord
    x = xyzG.col(i,1);
    y = xyzG.col(i,2);
    z = xyzG.col(i,3);
    
    if (isempty(find(Strips(:,2) == y & Strips(:,3) == z & x <= Strips(:,4) & x >= Strips(:,1))))
        % create new strip with y,z
        % Strips = xyz of first col point, x of 1st col point + chord, xyz
        % of quarter chord point, chord length
        s = s+1;
        Strips(s,1:3) = [x y z];
        Strips(s,4) = x + xyzG.chord(i);
        Strips(s,5:7) = (xyzG.crnr1(i,:) + xyzG.crnr2(i,:))*0.5 + [1 0 0]*xyzG.chord(i)*0.25;
        Strips(s,8)  = xyzG.chord(i);
        Strip.chord(s) = xyzG.chord(i);
        Strip.width(s) = sqrt( (xyzG.crnr2(i,2) - xyzG.crnr1(i,2))^2 + (xyzG.crnr2(i,3) - xyzG.crnr1(i,3))^2);
        
        Strip.quarterChordPoint(s,:) = Strips(s,5:7)*RotGeoStd;
        Strip.bound(s,:) = (xyzG.crnr2(i,:) - xyzG.crnr1(i,:))/norm((xyzG.crnr2(i,:) - xyzG.crnr1(i,:)))*RotGeoStd;
        
        % find the wing that corresponds to this strip s, in order to find
        % the twist of the strip (NOTE: this only works for wing with
        % constant twist, i.e. a constant offset angle
        wing = find(s <= wing_strip, 1);
        part = find(s <= part_strip(wing,:),1);
        Strip.twist(s) = twst_ang(wing,part);
        Strip.wing(s) = wing;
        Strip.part(s) = part;
        
        xyzG.Panel_strip(i) = s;
        xyzG.Strip_panel{s} = [i];
    else 
        s = find(Strips(:,2) == y & Strips(:,3) == z & x <= Strips(:,4) & x >= Strips(:,1));% & x > Strips(:,1) & x < Strips(:,4));
        xyzG.Panel_strip(i) = s;
        xyzG.Strip_panel{s} = [xyzG.Strip_panel{s} i];
    end
end

xyzG.Strips = Strips;

% Calculate for each strip the accoring P, Q matrices, and the
% corrseponding control surface
m = length(xyzG.Strip_panel);
for i = 1:m
    not_strip_panels = (xyzG.Panel_strip ~= i);
    strip_panels = (xyzG.Panel_strip == i);
    
    % Matrices for force and moment calculation
    GxV_strip = GxV;
    GyV_strip = GyV;
    GzV_strip = GzV;
    
    GxW_strip = GxW;
    GyW_strip = GyW;
    GzW_strip = GzW;
    
    Gx_gam_strip = Gx_gam;
    Gy_gam_strip = Gy_gam;
    Gz_gam_strip = Gz_gam;
    
    % Setting corresponging rows to zero
    GxV_strip(not_strip_panels,:) = 0;
    GyV_strip(not_strip_panels,:) = 0;
    GzV_strip(not_strip_panels,:) = 0;
    
    GxW_strip(not_strip_panels,:) = 0;
    GyW_strip(not_strip_panels,:) = 0;
    GzW_strip(not_strip_panels,:) = 0;
    
    Gx_gam_strip(not_strip_panels,:) = 0;
    Gy_gam_strip(not_strip_panels,:) = 0;
    Gz_gam_strip(not_strip_panels,:) = 0;
    
    Px_strip(:,:,i) = -J'*( [GxV_strip -GxW_strip zeros(n_panels,12*N_con_surf)] + Gx_gam_strip*J ); 
    Py_strip(:,:,i) = J'*( [GyV_strip -GyW_strip zeros(n_panels,12*N_con_surf)] + Gy_gam_strip*J ); 
    Pz_strip(:,:,i) = -J'*( [GzV_strip -GzW_strip zeros(n_panels,12*N_con_surf)] + Gz_gam_strip*J );
    
    Px_strip_red(:,:,i) = R'*Px_strip(:,:,i)*R;
    Py_strip_red(:,:,i) = R'*Py_strip(:,:,i)*R;
    Pz_strip_red(:,:,i) = R'*Pz_strip(:,:,i)*R;
    
    Gamma_total_strip(i,:) = strip_panels*J; 
    Gamma_total_strip_red(i,:) = Gamma_total_strip(i,:)*R; 
    
    
%     %% QAIR 1.4
%     T = zeros(1,n_panels);
%     T(strip_panels) = 1;
%     T = diag(T);
%     
%     Px_strip_new(:,:,i) =  J_new'*T*(Gx_Xtilde + Gx_gam_new*J_new);
%     Py_strip_new(:,:,i) =  J_new'*T*(Gy_Xtilde + Gy_gam_new*J_new);
%     Pz_strip_new(:,:,i) =  J_new'*T*(Gz_Xtilde + Gz_gam_new*J_new);
%     
%     Qx_strip_new(:,:,i) =  J_new'*T*(Hx_Xtilde + Hx_gam_new*J_new);
%     Qy_strip_new(:,:,i) =  J_new'*T*(Hy_Xtilde + Hy_gam_new*J_new);
%     Qz_strip_new(:,:,i) =  J_new'*T*(Hz_Xtilde + Hz_gam_new*J_new);
%     
%     Px_strip_red_new(:,:,i) = R_new'*Px_strip_new(:,:,i)*R_new;
%     Py_strip_red_new(:,:,i) = R_new'*Py_strip_new(:,:,i)*R_new;
%     Pz_strip_red_new(:,:,i) = R_new'*Pz_strip_new(:,:,i)*R_new;
%     
%     Qx_strip_red_new(:,:,i) = R_new'*Qx_strip_new(:,:,i)*R_new;
%     Qy_strip_red_new(:,:,i) = R_new'*Qy_strip_new(:,:,i)*R_new;
%     Qz_strip_red_new(:,:,i) = R_new'*Qz_strip_new(:,:,i)*R_new;
%     
%     Gamma_total_strip_new(i,:) = strip_panels*J_new; 
%     Gamma_total_strip_red_new(i,:) = Gamma_total_strip_new(i,:)*R_new;
%     
    
    % find the panel whose bound vortex is closest to the quarter chord of
    % the strip
    s = i;
    x_quarterChord = Strips(s,5);
    [x_err, index] = min(abs(x_quarterChord - xyzG.bound(xyzG.Strip_panel{s},1)));
    panel = xyzG.Strip_panel{s}(index);
    VelMat_Strips(:,:,s) = VelMat(:,:,panel);
    VelMat_Strips_red(:,:,s) = VelMat_Red(:,:,panel);
    
%     % QAIR 1.4
%     VelMat_Strips_new(:,:,s) = AeroMatrices.VelMat_new(:,:,panel);
%     VelMat_Strips_red_new(:,:,s) = AeroMatrices.VelMat_Red_new(:,:,panel);
    
%     'Strip'
%     s
%     'quarter chord'
%     100*norm(Strips(s,5:7) - xyzG.bound(panel,:))/Strip.chord(s)


    % find the control surface for each flap.  If no control surface, then
    % eqaual to 0
    Strip.con_surface(s) = 0;
    a = xyzG.Strip_panel{s}; 
    for j = 1:size(con_surf_DOF_pointer,1)
       b = con_surf_DOF_pointer(j,1):con_surf_DOF_pointer(j,2);
       c = intersect(a, b);
       if ~isempty(c)
           Strip.con_surface(s) = j;
           Strip.flap_axis_rot(:,s) = flap_axis_rot(j,:)*RotGeoStd;
           break
       end
    end
end
% toc


AeroMatrices.Gamma_total_strip = Gamma_total_strip;
AeroMatrices.Gamma_total_strip_red = Gamma_total_strip_red;
AeroMatrices.VelMat_Strips = VelMat_Strips;
AeroMatrices.VelMat_Strips_red = VelMat_Strips_red;

AeroMatrices.Px_strip = Px_strip;
AeroMatrices.Py_strip = Py_strip;
AeroMatrices.Pz_strip = Pz_strip;
AeroMatrices.Px_strip_red = Px_strip_red;
AeroMatrices.Py_strip_red = Py_strip_red;
AeroMatrices.Pz_strip_red = Pz_strip_red;
AeroMatrices.Strips = Strips;
AeroMatrices.Strip = Strip;

% % QAIR 1.4
% AeroMatrices.Gamma_total_strip_new = Gamma_total_strip_new;
% AeroMatrices.Gamma_total_strip_red_new = Gamma_total_strip_red_new;
% AeroMatrices.VelMat_Strips_new = VelMat_Strips_new;
% AeroMatrices.VelMat_Strips_red_new = VelMat_Strips_red_new;
% 
% AeroMatrices.Px_strip_new = Px_strip_new;
% AeroMatrices.Py_strip_new = Py_strip_new;
% AeroMatrices.Pz_strip_new = Pz_strip_new;
% AeroMatrices.Px_strip_red_new = Px_strip_red_new;
% AeroMatrices.Py_strip_red_new = Py_strip_red_new;
% AeroMatrices.Pz_strip_red_new = Pz_strip_red_new;
% 
% AeroMatrices.Qx_strip_new = Qx_strip_new;
% AeroMatrices.Qy_strip_new = Qy_strip_new;
% AeroMatrices.Qz_strip_new = Qz_strip_new;
% AeroMatrices.Qx_strip_red_new = Qx_strip_red_new;
% AeroMatrices.Qy_strip_red_new = Qy_strip_red_new;
% AeroMatrices.Qz_strip_red_new = Qz_strip_red_new;
% 
% AeroMatrices.Strips = Strips;
% AeroMatrices.Strip = Strip;

%% QAIR 1.2 - Construct force and moment matrices for control srufaces
J = AeroMatrices.J;
for i = 1:N_con_surf
    DOFs = con_surf_DOF_pointer(i,1):con_surf_DOF_pointer(i,2);  % DOFs corresponding to i-th control surface
    
    T = zeros(1,n_panels);
    T(DOFs) = 1;
    T = diag(T);
    
    AeroMatrices.Px_con_surf(:,:,i) =  -J'*T*([GxV -GxW zeros(n_panels,12*N_con_surf)] + Gx_gam*J);
    AeroMatrices.Py_con_surf(:,:,i) =   J'*T*([GyV -GyW zeros(n_panels,12*N_con_surf)] + Gy_gam*J);
    AeroMatrices.Pz_con_surf(:,:,i) =  -J'*T*([GzV -GzW zeros(n_panels,12*N_con_surf)] + Gz_gam*J);
    
    AeroMatrices.Qx_con_surf(:,:,i) =  -J'*T*([HxV -HxW zeros(n_panels,12*N_con_surf)] + Hx_gam*J);
    AeroMatrices.Qy_con_surf(:,:,i) =   J'*T*([HyV -HyW zeros(n_panels,12*N_con_surf)] + Hy_gam*J);
    AeroMatrices.Qz_con_surf(:,:,i) =  -J'*T*([HzV -HzW zeros(n_panels,12*N_con_surf)] + Hz_gam*J);
    
    r_bo_h1 =  -xyzG_Std.con_surf_hinge1(i,:);
    
    rx = r_bo_h1(1);
    ry = r_bo_h1(2);
    rz = r_bo_h1(3);
    
    AeroMatrices.Qx_con_surf(:,:,i) = AeroMatrices.Qx_con_surf(:,:,i) + ry*AeroMatrices.Pz_con_surf(:,:,i) - rz*AeroMatrices.Py_con_surf(:,:,i);
    AeroMatrices.Qy_con_surf(:,:,i) = AeroMatrices.Qy_con_surf(:,:,i) + rz*AeroMatrices.Px_con_surf(:,:,i) - rx*AeroMatrices.Pz_con_surf(:,:,i);
    AeroMatrices.Qz_con_surf(:,:,i) = AeroMatrices.Qz_con_surf(:,:,i) + rx*AeroMatrices.Py_con_surf(:,:,i) - ry*AeroMatrices.Px_con_surf(:,:,i);
    
    AeroMatrices.Px_con_surf_red(:,:,i) = R'*AeroMatrices.Px_con_surf(:,:,i)*R;
    AeroMatrices.Py_con_surf_red(:,:,i) = R'*AeroMatrices.Py_con_surf(:,:,i)*R;
    AeroMatrices.Pz_con_surf_red(:,:,i) = R'*AeroMatrices.Pz_con_surf(:,:,i)*R;
    
    AeroMatrices.Qx_con_surf_red(:,:,i) = R'*AeroMatrices.Qx_con_surf(:,:,i)*R;
    AeroMatrices.Qy_con_surf_red(:,:,i) = R'*AeroMatrices.Qy_con_surf(:,:,i)*R;
    AeroMatrices.Qz_con_surf_red(:,:,i) = R'*AeroMatrices.Qz_con_surf(:,:,i)*R;
    
    
    axis = acraft_dsgn.con_surf_axis_rot_Std(i,:);
    axis_x = axis(1);
    axis_y = axis(2);
    axis_z = axis(3);
    
    AeroMatrices.Torque_con_surf(:,:,i) = axis_x*AeroMatrices.Qx_con_surf(:,:,i) + ...
                                              axis_y*AeroMatrices.Qy_con_surf(:,:,i) + ...
                                              axis_z*AeroMatrices.Qz_con_surf(:,:,i);
                                            
    AeroMatrices.Torque_consurf_red(:,:,i) = R'*AeroMatrices.Torque_con_surf(:,:,i)*R;
end

% %% QAIR 1.4 - Construct force and moment matrices for control srufaces
% for i = 1:N_con_surf
%     DOFs = con_surf_DOF_pointer(i,1):con_surf_DOF_pointer(i,2);  % DOFs corresponding to i-th control surface
%     
%     T = zeros(1,n_panels);
%     T(DOFs) = 1;
%     T = diag(T);
%     
%     AeroMatrices.Px_con_surf_new(:,:,i) =  J_new'*T*(Gx_Xtilde + Gx_gam_new*J_new);
%     AeroMatrices.Py_con_surf_new(:,:,i) =  J_new'*T*(Gy_Xtilde + Gy_gam_new*J_new);
%     AeroMatrices.Pz_con_surf_new(:,:,i) =  J_new'*T*(Gz_Xtilde + Gz_gam_new*J_new);
%     
%     AeroMatrices.Qx_con_surf_new(:,:,i) =  J_new'*T*(Hx_Xtilde + Hx_gam_new*J_new);
%     AeroMatrices.Qy_con_surf_new(:,:,i) =  J_new'*T*(Hy_Xtilde + Hy_gam_new*J_new);
%     AeroMatrices.Qz_con_surf_new(:,:,i) =  J_new'*T*(Hz_Xtilde + Hz_gam_new*J_new);
%     
%     r_bo_h1 = -xyzG_Std.con_surf_hinge1(i,:);
%     
%     rx = r_bo_h1(1);
%     ry = r_bo_h1(2);
%     rz = r_bo_h1(3);
%     
%     AeroMatrices.Qx_con_surf_new(:,:,i) = AeroMatrices.Qx_con_surf_new(:,:,i) + ry*AeroMatrices.Pz_con_surf_new(:,:,i) - rz*AeroMatrices.Py_con_surf_new(:,:,i);
%     AeroMatrices.Qy_con_surf_new(:,:,i) = AeroMatrices.Qy_con_surf_new(:,:,i) + rz*AeroMatrices.Px_con_surf_new(:,:,i) - rx*AeroMatrices.Pz_con_surf_new(:,:,i);
%     AeroMatrices.Qz_con_surf_new(:,:,i) = AeroMatrices.Qz_con_surf_new(:,:,i) + rx*AeroMatrices.Py_con_surf_new(:,:,i) - ry*AeroMatrices.Px_con_surf_new(:,:,i);
%     
%     AeroMatrices.Px_con_surf_red_new(:,:,i) = R_new'*AeroMatrices.Px_con_surf_new(:,:,i)*R_new;
%     AeroMatrices.Py_con_surf_red_new(:,:,i) = R_new'*AeroMatrices.Py_con_surf_new(:,:,i)*R_new;
%     AeroMatrices.Pz_con_surf_red_new(:,:,i) = R_new'*AeroMatrices.Pz_con_surf_new(:,:,i)*R_new;
%     
%     AeroMatrices.Qx_con_surf_red_new(:,:,i) = R_new'*AeroMatrices.Qx_con_surf_new(:,:,i)*R_new;
%     AeroMatrices.Qy_con_surf_red_new(:,:,i) = R_new'*AeroMatrices.Qy_con_surf_new(:,:,i)*R_new;
%     AeroMatrices.Qz_con_surf_red_new(:,:,i) = R_new'*AeroMatrices.Qz_con_surf_new(:,:,i)*R_new;
%     
%     
%     axis = acraft_dsgn.con_surf_axis_rot_Std(i,:);
%     axis_x = axis(1);
%     axis_y = axis(2);
%     axis_z = axis(3);
%     
%     AeroMatrices.Torque_con_surf_new(:,:,i) = axis_x*AeroMatrices.Qx_con_surf_new(:,:,i) + ...
%                                               axis_y*AeroMatrices.Qy_con_surf_new(:,:,i) + ...
%                                               axis_z*AeroMatrices.Qz_con_surf_new(:,:,i);
%                                             
%     AeroMatrices.Torque_consurf_red_new(:,:,i) = R_new'*AeroMatrices.Torque_con_surf_new(:,:,i)*R_new;
% end


% % 
% % close all
% % figure
% % hold;
% % plot3(xyzG.crnr1(:,1),xyzG.crnr1(:,2),xyzG.crnr1(:,3),'o'); axis equal
% % plot3(xyzG.crnr2(:,1),xyzG.crnr2(:,2),xyzG.crnr2(:,3),'o');
% % plot3(xyzG.crnr3(:,1),xyzG.crnr3(:,2),xyzG.crnr3(:,3),'o');
% % plot3(xyzG.crnr4(:,1),xyzG.crnr4(:,2),xyzG.crnr4(:,3),'o');
% % plot3(xyzG.col(:,1),xyzG.col(:,2),xyzG.Panel_strip(:),'x');
% % grid

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xyzL = trapezoid_discretization(xyzL_pt1, xyzL_pt2, xyzL_pt3, xyzL_pt4, N, M, trpz, Flap_frac)

e1 = (xyzL_pt2 - xyzL_pt1)/N;   % vector going from point 1 to point 2
e2 = (xyzL_pt4 - xyzL_pt3)/N;   % vector going from point 3 to point 4

V = xyzL_pt3 - xyzL_pt1;        % vector going from point 1 to point 3

for i = 1:N*M
    Chrd_counter  = ceil(i/(N));           % calculate the chord-wise discrete "coordinate" of panel i, using the upward rounding function "ceil"
    Spn_counter  = i - N*(Chrd_counter-1); % calculate the span-wise discrete "coordinate" of panel i;
    
    xyzL.crnr1(i,:) = xyzL_ij(Chrd_counter-1, Spn_counter-1, xyzL_pt1, e1, e2, V, N, M) ;
    xyzL.crnr2(i,:) = xyzL_ij(Chrd_counter-1, Spn_counter, xyzL_pt1, e1, e2, V, N, M) ;
    xyzL.crnr3(i,:) = xyzL_ij(Chrd_counter, Spn_counter-1, xyzL_pt1, e1, e2, V, N, M) ;
    xyzL.crnr4(i,:) = xyzL_ij(Chrd_counter, Spn_counter, xyzL_pt1, e1, e2, V, N, M) ;

    xyzL.HS1(i,:) = xyzL.crnr1(i,:) + 0.25*(xyzL.crnr3(i,:) - xyzL.crnr1(i,:));
    xyzL.HS2(i,:) = xyzL.crnr2(i,:) + 0.25*(xyzL.crnr4(i,:) - xyzL.crnr2(i,:));
    
    xyzL_aux1  = xyzL.crnr1(i,:) + 0.5*(xyzL.crnr2(i,:) - xyzL.crnr1(i,:));
    xyzL_aux2  = xyzL.crnr3(i,:) + 0.5*(xyzL.crnr4(i,:) - xyzL.crnr3(i,:));
    
    xyzL.col(i,:)   = xyzL_aux1 + 0.75*(xyzL_aux2 - xyzL_aux1);
    
    xyzL.midpt(i,:) = xyzL_aux1 + 0.5*(xyzL_aux2 - xyzL_aux1);
    if trpz == 1
        xyzL.chord(i,:) = (norm(V + (e2 - e1)*(Spn_counter-1)) + norm(V + (e2 - e1)*(Spn_counter)))*0.5/(1-Flap_frac);
    else
        xyzL.chord(i,:) = (norm(V + (e2 - e1)*(Spn_counter-1)) + norm(V + (e2 - e1)*(Spn_counter)))*0.5/(Flap_frac);
    end
        
end



end


function xyzL_ij = xyzL_ij(i,j, xyzL_it1, e1, e2, V, N, M)
u = (V + (e2 - e1)*j)/M;
xyzL_ij = xyzL_it1 + u*i + e1*j;
end

function Normal_col = normal_col_discretization(Airfoil, X, Y, N, M, Flap_frac, trpz, Root_twst_ang, Tip_twst_ang, Dih_ang)

for i = 1:N*M
    Chrd_counter  = ceil(i/(N));           % calculate the chord-wise discrete "coordinate" of panel i, using the upward rounding function "ceil"
    Spn_counter  = i - N*(Chrd_counter-1); % calculate the span-wise discrete "coordinate" of panel i;
    
    % Calculating normalized x-coordinate at the panel's collocation
    % point
    if trpz == 1
        x_c = (1 - Flap_frac)*(Chrd_counter-0.25)/M;
    else 
        x_c = (1 - Flap_frac) + Flap_frac*(Chrd_counter-0.25)/M;
    end
    
    % Calculating the 0 twist 0 dihedral normal vector
    if Airfoil == 0
        normal_vect = [0 0 1];
    else
        normal_vect = normal_generator(X, Y, x_c);
    end
    
    % Calculating the twist angle at the span-wise station
    Twist_ang = Root_twst_ang + (Tip_twst_ang - Root_twst_ang)*(Spn_counter-0.5)/N;
    Rot_twist = [cos(Twist_ang)   0       -sin(Twist_ang)
                 0                1       0
                 sin(Twist_ang)   0       cos(Twist_ang)];

    Rot_dih   = [1   0                  0
                 0   cos(Dih_ang)       sin(Dih_ang)
                 0   -sin(Dih_ang)      cos(Dih_ang)];
    
    % Calculating the twisted with dihedral normal vector
    Normal_col(i,:) = normal_vect*Rot_twist*Rot_dih;
        
end
end

function normal_vect = normal_generator(X, Y, x_c)
% x_c is the normalized x-coordinate of the point 

for i = 1:length(X)
    if X(i) == x_c
        vect = [X(i+1)-X(i-1) 0 Y(i+1)-Y(i-1)];
        unit_vect = vect/norm(vect);    % This is the unit vector going from i-1 to i+1
        normal_vect = [-unit_vect(1,3) 0 unit_vect(1,1)];
    else if X(i) > x_c
            vect = [X(i)-X(i-1) 0 Y(i)-Y(i-1)];
            unit_vect = vect/norm(vect);    % This is the unit vector going from i-1 to i+1
            normal_vect = [-unit_vect(1,3) 0 unit_vect(1,1)];
            break;
        end
    end
end

end

function [X, Y] = arfoil_coord_parser(filename)

fid = fopen(filename);

% Disregard one line
fgetl(fid);

% Read the data
A = fscanf(fid,'%f');

n    = 2;
m    = length(A)/2;

A = reshape(A,n,m);

x1 = A(1,:);
y1 = A(2,:);

% Calculating the camberline coordinates X Y of the airfoil
i = 1;
while x1(i+1) < x1(i)
    [x, I] = min(abs(x1(i+1:length(x1))-x1(i)));
    I = I + i ;
    if x1(I) < x1(i)
        if I == length(x1)
            x1(I+1)=x1(I);
            y1(I+1)=y1(I);
        end
        y_aux = (y1(I+1)*(abs(x1(I)-x1(i))) + y1(I)*(abs(x1(I+1)-x1(i))))/abs(x1(I)-x1(I+1));
    else
        y_aux = (y1(I-1)*(abs(x1(I)-x1(i))) + y1(I)*(abs(x1(I-1)-x1(i))))/abs(x1(I)-x1(I-1));
    end
    x2(i) = x1(i);
    y2(i) = (y1(i)+y_aux)*0.5;
    
    i = i + 1;
end

% Ordering the coordinates from LE to TE
j = length(x2);
for i = 1:length(x2)
    X(j) = x2(i);
    Y(j) = y2(i);
    j = j - 1;
end
% 
% close all
% figure; 
% plot(X,Y,x1,y1,'o')
% axis equal; 
% grid

end

function [u, v, w] = HS_vortex(xp, yp, zp, x1, y1, z1, x2, y2, z2, gamma, x_inf, Index_p, Index_HS, chord)
% Calculating finite core radius
mag_ro = sqrt((y1 - y2)^2 + (z1 - z2)^2 );
if Index_p == Index_HS
    core_radius = mag_ro*1e-4;
else
    core_radius = max(chord/4,0.5*mag_ro);
end

%  Velocity induced by Y-minus seminfinite vortex filament on point P
%  (placing the x-infinity point at 20 times the chord length)
[u1, v1, w1] = vortex_segment(xp, yp, zp, x_inf, y1, z1, x1, y1, z1, gamma, core_radius);

%  Velocity induced by the bounded vortex segment on point P
[u2, v2, w2] = vortex_segment(xp, yp, zp, x1, y1, z1, x2, y2, z2, gamma, core_radius);

%  Velocity induced by Y-plus seminfinite vortex filament on point P
%  (placing the x-infinity point at 20 times the chord length)
[u3, v3, w3] = vortex_segment(xp, yp, zp, x2, y2, z2, x_inf, y2, z2, gamma, core_radius);

u = u1 + u2 + u3;
v = v1 + v2 + v3;
w = w1 + w2 + w3;
end

function [u, v, w] = vortex_segment(xp, yp, zp, x1, y1, z1, x2, y2, z2, gamma, RCORE)

% r1_cross_r2_x =  (yp - y1)*(zp - z2) - (zp - z1)*(yp - y2);
% r1_cross_r2_y = -(xp - x1)*(zp - z2) + (zp - z1)*(xp - x2);
% r1_cross_r2_z =  (xp - x1)*(yp - y2) - (yp - y1)*(xp - x2);
% mag_r1_cross_r2_square = r1_cross_r2_x^2 + r1_cross_r2_y^2 + r1_cross_r2_z^2;
% 
% ro_dot_r1 = (x2 - x1)*(xp - x1) + (y2 - y1)*(yp - y1) + (z2 - z1)*(zp - z1);
% ro_dot_r2 = (x2 - x1)*(xp - x2) + (y2 - y1)*(yp - y2) + (z2 - z1)*(zp - z2);
% mag_r1 = sqrt( (xp - x1)^2 + (yp - y1)^2 + (zp - z1)^2 );
% mag_r2 = sqrt( (xp - x2)^2 + (yp - y2)^2 + (zp - z2)^2 );
% mag_ro = sqrt( (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2 );
% 
% K = ( ro_dot_r1/mag_r1 - ro_dot_r2/mag_r2 )*gamma/(4*pi*mag_r1_cross_r2_square);
% 
% %lambda = ro_dot_r1/mag_ro;
% %d = sqrt( (xp - x1 - (x2-x1)*lambda/mag_ro)^2 + (yp - y1 - (y2-y1)*lambda/mag_ro)^2 + (zp - z1 - (z2-z1)*lambda/mag_ro)^2 );
% d = sqrt(mag_r1_cross_r2_square)/mag_ro;
% if d > RCORE
%     u = K*r1_cross_r2_x;
%     v = K*r1_cross_r2_y;
%     w = K*r1_cross_r2_z;
% else 
%     u = 0;
%     v = 0;
%     w = 0;
% end
% 
% 
% %% AVL Method for induced velocities, including finite core radius, with
% %% modifications for more efficiency
% Ax =  x1 - xp;
% Ay =  y1 - yp;
% Az =  z1 - zp;
% 
% Bx =  x2 - xp;
% By =  y2 - yp;
% Bz =  z2 - zp;
% 
% % ASQ = Ax*Ax + Ay*Ay + Az*Az;
% % BSQ = Bx*Bx + By*By + Bz*Bz;
% 
% % AMAG = sqrt(Ax*Ax + Ay*Ay + Az*Az);
% % BMAG = sqrt(Bx*Bx + By*By + Bz*Bz);
% 
% 
% 
% if (sqrt(Ax*Ax + Ay*Ay + Az*Az))*(sqrt(Bx*Bx + By*By + Bz*Bz)) ~= 0
% %     AXBx = Ay*Bz - Az*By;
% %     AXBy = Az*Bx - Ax*Bz;
% %     AXBz = Ax*By - Ay*Bx;
% %     AXBSQ = (Ay*Bz - Az*By)^2 + (Az*Bx - Ax*Bz)^2 + (Ax*By - Ay*Bx)^2;
% % 
% %     ADB = Ax*Bx + Ay*By + Az*Bz;
% %     ALSQ = Ax*Ax + Ay*Ay + Az*Az + Bx*Bx + By*By + Bz*Bz - 2*Ax*Bx + Ay*By + Az*Bz;
% %     AB = AMAG*BMAG;
%     
%     T = ( (Bx*Bx + By*By + Bz*Bz-Ax*Bx+Ay*By+Az*Bz)/sqrt(Bx*Bx + By*By + Bz*Bz+RCORE^2) + (Ax*Ax+Ay*Ay+Az*Az-Ax*Bx + Ay*By + Az*Bz)/sqrt(Ax*Ax+Ay*Ay+Az*Az+RCORE^2) ) / ((Ay*Bz - Az*By)^2 + (Az*Bx - Ax*Bz)^2 + (Ax*By - Ay*Bx)^2 + (Ax*Ax + Ay*Ay + Az*Az + Bx*Bx + By*By + Bz*Bz - 2*Ax*Bx + Ay*By + Az*Bz)*RCORE^2);
%     
%     u = (Ay*Bz - Az*By)*T*gamma/(4*pi);
%     v = (Az*Bx - Ax*Bz)*T*gamma/(4*pi);
%     w = (Ax*By - Ay*Bx)*T*gamma/(4*pi);
% else
%     u = 0;
%     v = 0;
%     w = 0;
% end


%% AVL Method for induced velocities, including finite core radius
Ax =  x1 - xp;
Ay =  y1 - yp;
Az =  z1 - zp;

Bx =  x2 - xp;
By =  y2 - yp;
Bz =  z2 - zp;

ASQ = Ax^2 + Ay^2 + Az^2;
BSQ = Bx^2 + By^2 + Bz^2;

AMAG = sqrt(ASQ);
BMAG = sqrt(BSQ);

u = 0;
v = 0;
w = 0;

if AMAG*BMAG ~= 0
    AXBx = Ay*Bz - Az*By;
    AXBy = Az*Bx - Ax*Bz;
    AXBz = Ax*By - Ay*Bx;
    AXBSQ = AXBx^2 + AXBy^2 + AXBz^2;

    ADB = Ax*Bx + Ay*By + Az*Bz;
    ALSQ = ASQ + BSQ - 2*ADB;
%     AB = AMAG*BMAG;
    
    T = ( (BSQ-ADB)/sqrt(BSQ+RCORE^2) + (ASQ-ADB)/sqrt(ASQ+RCORE^2) ) / (AXBSQ + ALSQ*RCORE^2);
    
    u = AXBx*T*gamma/(4*pi);
    v = AXBy*T*gamma/(4*pi);
    w = AXBz*T*gamma/(4*pi);
end
end

function Strip_index = assign_strip(xyz, DOFs)
% y is a vector with the y position if the "corner 1" of each HS vortex, of
% each wing
% DOFs is a vector tha goes along with y, and indicates what is the global
% DOF of each of the corresponding HS vor% tices
% The idea is, if two vor% tices have the same y coordinate, then they belong
% to the same strip

% First, we will sort vector y, in an ascending manner
y = xyz(:,2);
[y_sorted, DOF_sorted] = sort(y);
Strip = 0;
for i = 1:length(DOFs)
    if i == 1
        Strip = Strip +1;
        
        Strip_index(DOF_sorted(i)) = 1;
        
    else
        if y_sorted(i) == y_sorted(i-1)
            Strip_index(DOF_sorted(i)) = Strip;
        else
            Strip = Strip +1;
            Strip_index(DOF_sorted(i)) = Strip;
        end
    end
end
end

