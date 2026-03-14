function [HBSE] = build_HBSE(vc_band_BSE, Nv, Nc, N_kpt, Number_par, ...
                             Ec_kex, Ev, k, b1, b2, d, q_data, epsQ2D_data, ...
                             Wq0, wi, tau, Basis_tau_index_cell, WF_center, N_atoms, ...
                             Cc_kex, Cv, e, eps0)
%
% build_HBSE: Encapsulated function to build the HBSE matrix.
%
% *** FIX v1.1: Modified to accept Basis_tau_index as a CELL array ***
%               (Changed Basis_tau_index(Ni).data to Basis_tau_index_cell{Ni})
%

% Ensure vc_band_BSE is oriented correctly (MATLAB engine can pass vectors as rows)
if size(vc_band_BSE, 2) > size(vc_band_BSE, 1)
    vc_band_BSE = vc_band_BSE';
end

HBSE_size = Nc * Nv * N_kpt;
HBSE = zeros(HBSE_size, HBSE_size);

for vvcci = 1:size(vc_band_BSE,1)
    vv = vc_band_BSE(vvcci,1); 
    cc = vc_band_BSE(vvcci,2);
    
    % Display progress
    if vvcci == 1
        disp(['vv = ', num2str(vv), ', cc = ', num2str(cc)]);
        loop_timer = tic;
    elseif mod(vvcci, 10) == 0 || vvcci == size(vc_band_BSE,1)
        elapsed = toc(loop_timer);
        disp(['vv = ', num2str(vv), ', cc = ', num2str(cc), ...
              ' (Processed ', num2str(vvcci), '/', num2str(size(vc_band_BSE,1)), ...
              ' blocks in ', num2str(elapsed, '%.1f'), 's)']);
        loop_timer = tic;
    end
    
    for vci = vvcci:size(vc_band_BSE,1)
        v = vc_band_BSE(vci,1); 
        c = vc_band_BSE(vci,2);

        HBSE_element = zeros(N_kpt, N_kpt);
        
        parfor (j = 1:N_kpt, Number_par)
            KE = zeros(N_kpt,1);
            KE(j) = (Ec_kex(j,c) - Ev(j,v)) * ( (c==cc) * (v==vv) );
            q = bsxfun(@minus, k, k(j,:));
            
            [Vq, epsilon_q] = Screened_Coulomb_interaction(q, b1, b2, d, q_data, epsQ2D_data);
            
            % --- THIS IS THE MODIFIED SECTION ---
            % We must call the sub-function Direct_Coulomb_interaction
            % with the new cell array name 'Basis_tau_index_cell'
            Vd = Direct_Coulomb_interaction(Vq, epsilon_q, k, b1, b2, Wq0(j), wi(j), q, ...
                tau, Basis_tau_index_cell, WF_center, N_atoms, Cc_kex, Cv, c, cc, v, vv, j, e, eps0); % Modified by Oscar
            % --- END MODIFIED SECTION ---

            HBSE_element(:, j) = KE + (wi(j) / (2*pi)^2) * Vd;
        end

        idx1 = (1:N_kpt) + (((c-1) + (v-1)*Nc) * N_kpt);
        idx2 = (1:N_kpt) + (((cc-1) + (vv-1)*Nc) * N_kpt);
        
        HBSE(idx1, idx2) = HBSE_element;

        if v ~= vv || c ~= cc
            parfor (j = 1:N_kpt, Number_par)
                HBSE_element(:,j) = HBSE_element(:,j) ./ (wi(j) / (2*pi)^2);
            end
            HBSE(idx2, idx1) = (transpose(wi) / (2*pi)^2) .* ctranspose(HBSE_element);
        end
    end
end

end

