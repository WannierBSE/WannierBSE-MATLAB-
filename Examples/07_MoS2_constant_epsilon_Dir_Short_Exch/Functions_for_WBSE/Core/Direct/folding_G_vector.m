function [Gi] = folding_G_vector(b1, b2)
%==========================================================================
% FOLDING_G_VECTOR
%
% DESCRIPTION:
%   Returns the nearest-cell reciprocal vectors used to fold momentum transfers into the first Brillouin zone.
%
% INPUT ARGUMENTS:
%   b1, b2 : reciprocal lattice vectors.
%
% OUTPUT:
%   Gi : candidate reciprocal lattice shifts.
%==========================================================================

b1 = reshape(b1, [1, 2]);
b2 = reshape(b2, [1, 2]);

G0 = [0, 0];
G1 = b1;
G2 = b2;
G3 = b1 + b2;
G4 = -b1;
G5 = -b2;
G6 = -b1 - b2;
G7 = b1 - b2;
G8 = -b1 + b2;

% Aggregate into a single matrix (9 x 2)
Gi = [G0; G1; G2; G3; G4; G5; G6; G7; G8];

end
