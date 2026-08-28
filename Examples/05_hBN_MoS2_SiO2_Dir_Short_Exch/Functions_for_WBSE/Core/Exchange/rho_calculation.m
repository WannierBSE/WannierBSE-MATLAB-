function [rho] = rho_calculation(W_up,W_down,W_up_t,W_down_t,G_all,R_mNN,n1_n2_mNN,...
                                        zcut_min,zcut_max,Na1,Na2,Na3,shift_a1,shift_a2,shift_a3,tau,...
                                        N_atoms,d_max,Basis_tau_index,Mmax_tolerance,tin_all,Nz_t,Na1_t,Na2_t,...
                                        kex,a1,a2,u_point_t,v_point_t,u_t,hu,hv,J,kexTol,Npar_rho)
%==========================================================================
% RHO_CALCULATION
%
% DESCRIPTION:
%   Computes Fourier-transformed Wannier pair densities needed for the short-range exchange tensor.
%
% OUTPUT:
%   rho : pair-density tensor indexed by Wannier basis, R shell, and G vector.
%==========================================================================



[hfuv,hfuv0] = spline_integral_inplane(); % functions for the spline interpolation

%========= rho ==========
rho = zeros(size(W_up_t,2),size(W_up_t,2),length(zcut_min:zcut_max),size(R_mNN,1),size(G_all,1));
for Ri=1:size(R_mNN,1)
R = R_mNN(Ri,:);

%------------W(r-R)----------
n1 = n1_n2_mNN(Ri,1);
n2 = n1_n2_mNN(Ri,2);

Wl_up_t = zeros(size(W_up_t));
Wl_down_t = zeros(size(W_up_t));
for wj=1:size(W_up,2)
Reshape_W_up = reshape(W_up(:,wj),[Na1,Na2,Na3]);
Reshape_W_up = circshift(Reshape_W_up,[n1*shift_a1,n2*shift_a2,shift_a3]);
Wl_up = reshape(Reshape_W_up,[Na1*Na2*Na3,1]);
Wl_up_t(:,wj) = Wl_up(tin_all);

Reshape_W_down = reshape(W_down(:,wj),[Na1,Na2,Na3]);
Reshape_W_down = circshift(Reshape_W_down,[n1*shift_a1,n2*shift_a2,shift_a3]);
Wl_down = reshape(Reshape_W_down,[Na1*Na2*Na3,1]);
Wl_down_t(:,wj) = Wl_down(tin_all);
end
%------------W(r-R)----------

[j_index_tot] = find_j_index_in_rho_integral(R,tau,N_atoms,d_max,Basis_tau_index,Mmax_tolerance);

%-------- inplane integration -------
for jj=1:size(W_up_t,2)
W0i_up_t = W_up_t(:,jj);
W0i_down_t = W_down_t(:,jj);

j_index = j_index_tot{jj};
for j=j_index
WRj_up_t = Wl_up_t(:,j);
WRj_down_t = Wl_down_t(:,j);
Wij = (conj(W0i_up_t).*WRj_up_t)+(conj(W0i_down_t).*WRj_down_t);

%--------integral with given z and G-------
rho_temp = zeros(Nz_t,size(G_all,1));
parfor (zi=1:Nz_t,Npar_rho)
ind_z = (1:(Na1_t*Na2_t))+(zi-1)*(Na1_t*Na2_t);
u_2D = u_t(ind_z);
Wij_z = Wij(ind_z);
[c3,c2,c1,c0] = wf_spline_interpolation(Wij_z,kex,a1,a2,u_2D,u_point_t,v_point_t,kexTol);

rho_temp_G = zeros(1,size(G_all,1));
for g=1:size(G_all,1)
G = G_all(g,:);
rho_temp_G(g) = rho_integral_v2(hfuv,hfuv0,c3,c2,c1,c0,kex,a1,a2,u_2D,u_point_t,v_point_t,kexTol,hu,hv,J,G);
end
rho_temp(zi,:) = rho_temp_G;
end
%--------integral with given z and G-------
%--------assign values-----
    for Bz = 1:Nz_t
    for Bg = 1:size(G_all,1)
    rho(jj,j,Bz,Ri,Bg) = rho_temp(Bz,Bg);
    end
    end
%--------assign values-----

end
end
%-------- inplane integration -------

clear Wl_up_t Wl_down_t
end

%========= rho ==========



end
