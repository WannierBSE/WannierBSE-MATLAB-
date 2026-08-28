function [Vq, epsilon_q] = screened_coulomb_interaction(q, b1, b2, d, q_data, epsQ2D_data)
%==========================================================================
% SCREENED_COULOMB_INTERACTION
%
% DESCRIPTION:
%   Evaluates the folded z-averaged Coulomb kernel and effective dielectric screening for all momentum transfers.
%
% INPUT ARGUMENTS:
%   q : momentum-transfer vectors.
%   b1, b2 : reciprocal lattice vectors.
%   d : effective material thickness.
%   q_data, epsQ2D_data : tabulated dielectric screening model.
%
% OUTPUT:
%   Vq : folded bare quasi-2D Coulomb kernel.
%   epsilon_q : effective dielectric values at folded q.
%==========================================================================

[q_mag, ~] = wave_vector_folding(q, b1, b2);

% --- Dielectric Screening ---
% Interpolate the dielectric function for the given q-magnitudes
epsilon_q = epsilon_Q2D(q_data, epsQ2D_data, q_mag);

% --- 2D Screened Coulomb Potential ---
% Analytical form for the quasi-2D Coulomb interaction
Vq = (1 - ((1 - exp(-(d .* q_mag))) ./ (d .* q_mag))) .* ((4 * pi) ./ (d .* (q_mag.^2)));

end
