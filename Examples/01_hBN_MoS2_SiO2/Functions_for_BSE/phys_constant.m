function [e,eps0] = phys_constant
load('physical_constant.mat');
e_macro = e;
eps0_macro = eps0;

% change unit
clear e eps0 hbar m0
e = e_macro./e_macro;
eps0 = eps0_macro*(1/(e_macro*1e10));  % (e^2)/(eV*angstrom)
end