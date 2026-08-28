function Vd = direct_coulomb_interaction(Vq, epslnq, k, b1, b2, Wq0_j, wi_j, q, tau, Basis_tau_index_cell, WF_center, N_atoms, Cc_kex, Cv, c, cc, v, vv, j, e, eps0)
%==========================================================================
% DIRECT_COULOMB_INTERACTION
%
% DESCRIPTION:
%   Computes the long-range screened direct electron-hole Coulomb interaction vector for one BSE band-pair channel.
%
% INPUT ARGUMENTS:
%   Vq, epslnq : folded bare Coulomb kernel and dielectric screening values.
%   k : k-point coordinates.
%   b1, b2 : reciprocal lattice vectors.
%   Wq0_j, wi_j : singular-cell correction and k-weight for k-point j.
%   q : momentum-transfer vectors.
%   tau : atomic positions.
%   Basis_tau_index_cell, WF_center : Wannier-center to atom mappings.
%   N_atoms : number of atoms.
%   Cc_kex, Cv : conduction and valence TB coefficients.
%   c, cc, v, vv : band indices.
%   j : reference k-point index.
%   e, eps0 : WBSE unit constants.
%
% OUTPUT:
%   Vd : direct Coulomb contribution for all k-points.
%==========================================================================

DEBUG_VD = false; % Diagnostic flag

% --- Wave Vector Folding ---
[~, MinNdkG] = wave_vector_folding(q, b1, b2);
[Gi] = folding_G_vector(b1, b2);

% --- Phase Correction (Tau) ---
% Transform coefficients to include atomic phase factors
[Cc_kex_tilde, Cv_tilde] = phase_tau_in_Vd(k, Cc_kex, Cv, tau, WF_center);

% --- G dot tau Phase factors ---
Gkkp = Gi(MinNdkG, :);
expiGtau = zeros(size(Vq, 1), size(Cc_kex_tilde, 2));
for Ni = 1:N_atoms
    % Apply atomic phases to corresponding Wannier functions
    expiGtau(:, Basis_tau_index_cell{Ni}) = repmat(exp(1i * sum(bsxfun(@times, [Gkkp(:, 1), Gkkp(:, 2), zeros(size(k, 1), 1)], tau(Ni, :)), 2)), ...
                                              1, length(Basis_tau_index_cell{Ni}));
end

% --- Matrix Element Construction (Overlap Integrals) ---
Mc = sum((conj(Cc_kex_tilde(:, :, c)) .* expiGtau) .* Cc_kex_tilde(j, :, cc), 2);
Mv = sum((conj(Cv_tilde(:, :, v)) .* expiGtau) .* Cv_tilde(j, :, vv), 2);

% Bare Direct interaction
Vd = -((e^2) ./ (4 * pi * eps0)) .* (Vq ./ epslnq) .* Mc(:) .* conj(Mv(:));

% --- Singularity Handling (k=k' term) ---
Vd(j) = Vd_singularity(e, eps0, Wq0_j, wi_j, Cc_kex_tilde, Cv_tilde, c, cc, v, vv, j);

end
