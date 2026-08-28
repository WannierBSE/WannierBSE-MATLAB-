function [chi] = chi_calculation(rho,R_mNN,G_all,kex,tau,N_atoms,d_max,Basis_tau_index,Mmax_tolerance,Npar_Xi,e,eps0,z_point_t,b1)
%==========================================================================
% CHI_CALCULATION
%
% DESCRIPTION:
%   Accumulates the short-range exchange tensor chi from pair-density tensors and spline-evaluated Xi integrals.
%
% OUTPUT:
%   chi : exchange tensor over Wannier basis and retained lattice translations.
%==========================================================================


[hfz,hfz0] = spline_integral_z();

% About the Xi, please refer to section 2.5 of " BSE_documentation.ppt " file for details.

%========= Xi and chi ==========
chi = zeros(size(rho,1),size(rho,1),size(rho,1),size(rho,1),size(R_mNN,1),size(R_mNN,1));
parfor (g=1:size(G_all,1),Npar_Xi)
    G = G_all(g,:);
    Xi_g = zeros(size(rho,1),size(rho,1),size(rho,1),size(rho,1),size(R_mNN,1),size(R_mNN,1));

    % STRICTLY SHORT RANGE: G=0 is excluded in generate_g_vectors.m
    % We only process the non-zero momentum cases.
    V_kex_G = ((e^2)/(2*eps0))*(1./sqrt(sum((kex+[G,0]).^2)));

    for R_ind1=1:size(R_mNN,1)
        jj_index = 1:size(rho,1);
        R1 = R_mNN(R_ind1,:);
        j_index_tot = find_j_index_in_rho_integral(R1,tau,N_atoms,d_max,Basis_tau_index,Mmax_tolerance);
        
        for R_ind2=1:size(R_mNN,1)
            ll_index = 1:size(rho,1);
            R2 = R_mNN(R_ind2,:);
            l_index_tot = find_j_index_in_rho_integral(R2,tau,N_atoms,d_max,Basis_tau_index,Mmax_tolerance);

            % Xi of kex+G unequal to 0 
            Xi_g(:,:,:,:,R_ind1,R_ind2) = xi_integral(hfz0,hfz,rho(:,:,:,R_ind1,g),rho(:,:,:,R_ind2,g),G,kex,z_point_t,size(rho,1),b1,jj_index,j_index_tot,ll_index,l_index_tot);
        end
    end

    chi = chi + V_kex_G*Xi_g;
end
%========= Xi and chi ==========

end
