function [q_mag, MinNdkG] = wave_vector_folding(q, b1, b2)
%==========================================================================
% WAVE_VECTOR_FOLDING
%
% DESCRIPTION:
%   Folds momentum-transfer vectors by selecting the reciprocal-lattice shift with minimum magnitude.
%
% INPUT ARGUMENTS:
%   q : momentum-transfer vectors.
%   b1, b2 : reciprocal lattice vectors.
%
% OUTPUT:
%   q_mag : folded momentum magnitudes.
%   MinNdkG : index of the selected reciprocal shift.
%==========================================================================

[Gi] = folding_G_vector(b1, b2);

% --- Folding Logic ---
% Calculate q + G for all candidate G-vectors
qG = zeros(size(q, 1), 2, size(Gi, 1));
for gi = 1:size(Gi, 1)
    qG(:, :, gi) = q + Gi(gi, :);
end

% Compute magnitudes and find the minimum (closest to Gamma)
NqG = squeeze(sqrt(sum(qG.^2, 2)));
[q_mag, MinNdkG] = min(NqG, [], 2);

end
