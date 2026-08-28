function [hf, hf0] = SPint_z
%==========================================================================
% SPINT_Z
%
% DESCRIPTION:
%   Returns analytic spline integration kernels for z-direction polynomial and exponential factors in short-range exchange integrals.
%
% OUTPUT:
%   hf : exponential z-integral kernel.
%   hf0 : non-oscillatory polynomial integral kernel.
%==========================================================================

syms f(c1,c2,c3,c4,G,y0,h,yc) f0(c1,c2,c3,c4,y0,h) y

% --- Analytical Integration ---
% Define the cubic polynomial: p(y) = c1(y-y0)^3 + c2(y-y0)^2 + c3(y-y0) + c4

% 1. Integration with exponential factor
f(c1,c2,c3,c4,G,y0,h,yc) = int(((c1.*((y-y0).^3) + c2.*((y-y0).^2) + c3.*(y-y0) + c4) .* exp(-G.*(y-yc))), y, [y0, y0+h]);

% 2. Standard polynomial integration (Normalization)
f0(c1,c2,c3,c4,y0,h) = int(((c1.*((y-y0).^3) + c2.*((y-y0).^2) + c3.*(y-y0) + c4)), y, [y0, y0+h]);

% --- Optimization & Conversion ---
% Simplify the symbolic expressions for faster execution during the main loops
f  = simplify(f);
f0 = simplify(f0);

% Convert to anonymous function handles for numerical use in Step 4/BSE solvers
hf  = matlabFunction(f);
hf0 = matlabFunction(f0);

end
