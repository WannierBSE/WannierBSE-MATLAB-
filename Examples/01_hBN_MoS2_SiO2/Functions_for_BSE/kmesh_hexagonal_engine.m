function [final_kmesh, total_coarse, total_dense, orientation] = kmesh_hexagonal_engine(Nk_coarse, Nk_dense, Dense_range, a1, a2, a_const)
%KMESH_HEXAGONAL_ENGINE
%   Core geometric processing engine for dual-resolution hexagonal Brillouin 
%   zones. Generates dimensionless k-space coordinates based on symmetry.
%
%   Logic: Bit-perfect geometric reflections (Dimensionless processing).

    orientation = get_hex_orientation(a1, a2);
    
    % Define Basis Vectors U, V, and W (Dimensionless/Normalized)
    if strcmp(orientation, 'Horizontal')
        U = [2/sqrt(3), 2/3];
        V = [2/sqrt(3), -2/3];
        W = [0, -4/3];
    elseif strcmp(orientation, 'Vertical')
        U = [2/3, 2/sqrt(3)];
        V = [-2/3, 2/sqrt(3)]; 
        W = [-4/3, 0];
    else
        error('WannierBSE:LatticeError', 'Unsupported orientation: %s.', orientation);
    end

    % Geometric Increments (Coarse)
    u = U / Nk_coarse; 
    v = V / Nk_coarse; 
    w = W / Nk_coarse;

    % Geometric Increments (Dense)
    u_d = (U * (Dense_range / Nk_coarse)) / Nk_dense;
    v_d = (V * (Dense_range / Nk_coarse)) / Nk_dense;
    w_d = (W * (Dense_range / Nk_coarse)) / Nk_dense;

    % Local Origin Shifts for Dense Clusters
    Shift_U = U * (1 - Dense_range / Nk_coarse);
    Shift_V = V * (1 - Dense_range / Nk_coarse);
    Shift_W = W * (1 - Dense_range / Nk_coarse);

    limit_lower = Nk_coarse - Dense_range;

    % --- Slice 1: Segment Gamma-U-V ---
    [k_low1, ~] = triangular_kmesh(0, Nk_coarse-1, u, v, [0,0], limit_lower, 'lower');
    [k_up1,  ~] = triangular_kmesh(0, Nk_coarse-1, u, v, [0,0], limit_lower, 'upper');
    [k_dU1,  ~] = triangular_kmesh(0, Nk_dense-1, u_d, v_d, Shift_U, [], 'all');
    [k_dV1,  ~] = triangular_kmesh(0, Nk_dense-1, u_d, v_d, Shift_V, [], 'all');
    master_UV = [k_low1; k_up1; k_dU1; k_dV1];

    % --- Slice 2: Segment Gamma-V-W ---
    [k_low2, ~] = triangular_kmesh(0, Nk_coarse-1, v, w, [0,0], limit_lower, 'lower');
    [k_up2,  ~] = triangular_kmesh(0, Nk_coarse-1, v, w, [0,0], limit_lower, 'upper');
    [k_dV2,  ~] = triangular_kmesh(0, Nk_dense-1, v_d, w_d, Shift_V, [], 'all');
    [k_dW2,  ~] = triangular_kmesh(0, Nk_dense-1, v_d, w_d, Shift_W, [], 'all');
    master_VW = [k_low2; k_up2; k_dV2; k_dW2];

    % Symmetry Operators (Reflection Sign Flips)
    if strcmp(orientation, 'Vertical')
        M = { [1, 1], [1, 1], [1, -1], [1, -1], [-1, -1], [-1, 1] };
    else
        M = { [1, 1], [1, 1], [-1, 1], [-1, 1], [-1, -1], [1, -1] };
    end
    
    % Hexagonal Assembly (6-Slice Reconstruction)
    S = { master_UV, master_VW, master_VW, master_UV, master_VW, master_VW };
    final_kmesh = [];
    for i = 1:6
        final_kmesh = [final_kmesh; S{i} .* M{i}];
    end

    % Statistical Point Counts
    num_coarse_UV = size([k_low1; k_up1], 1);
    num_dense_UV  = size([k_dU1; k_dV1], 1);
    num_coarse_VW = size([k_low2; k_up2], 1);
    num_dense_VW  = size([k_dV2; k_dW2], 1);
    
    total_coarse = (num_coarse_UV * 2) + (num_coarse_VW * 4);
    total_dense  = (num_dense_UV * 2) + (num_dense_VW * 4);
end

function orientation = get_hex_orientation(a1, a2)
    % Validates hexagonal symmetry and primary axis orientation.
    a = norm(a1); tol = 1e-2;
    if abs(a - norm(a2)) > tol || ~any(abs(rad2deg(acos(dot(a1,a2)/(a*a))) - [60, 120]) < tol)
        orientation = 'Not Hexagonal'; return;
    end
    test_list = {a1, a2, a1+a2, a1-a2};
    for i = 1:length(test_list)
        v = test_list{i};
        if abs(v(1)) < tol && abs(abs(v(2)) - a) < tol
            orientation = 'Horizontal'; return;
        elseif abs(v(2)) < tol && abs(abs(v(1)) - a) < tol
            orientation = 'Vertical'; return;
        end
    end
    orientation = 'Neither';
end

function [k_list, mesh_lines] = triangular_kmesh(N_start, N_end, u, v, shift, limit, mode)
    % Core mesh generator for triangular sub-zones using centroid placement.
    % Logic uses integer-based offsets to ensure compatibility with package standards.
    k_list = []; mesh_lines = [];
    for i = N_start : N_end
        for j = 0:i
            % Downward Triangle Centroid Calculation
            P1 = (i-j)*u + j*v + shift;
            P2 = (i-j+1)*u + j*v + shift;
            P3 = (i-j)*u + (j+1)*v + shift;
            c_up = (P1 + P2 + P3) / 3.0;
            
            keep_up = false;
            if strcmp(mode, 'all'), keep_up = true;
            elseif strcmp(mode, 'lower'), if i < limit, keep_up = true; end
            elseif strcmp(mode, 'upper'), if i >= limit && ~(j > limit - 1 || (i-j) > limit - 1), keep_up = true; end
            end
            if keep_up, k_list = [k_list; c_up]; end
            
            % Upward Triangle Centroid Calculation
            if (i + 1 <= N_end)
                Pd1 = (i-j+1)*u + j*v + shift;
                Pd2 = (i-j)*u + (j+1)*v + shift;
                Pd3 = (i-j+1)*u + (j+1)*v + shift;
                c_down = (Pd1 + Pd2 + Pd3) / 3.0;
                
                keep_down = false;
                if strcmp(mode, 'all'), keep_down = true;
                elseif strcmp(mode, 'lower'), if (i + 1 < limit), keep_down = true; end
                elseif strcmp(mode, 'upper'), if (i + 1 >= limit) && (j <= limit - 1) && ((i-j) <= limit - 1), keep_down = true; end
                end
                if keep_down, k_list = [k_list; c_down]; end
            end
        end
    end
end