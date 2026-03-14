function Vd = Direct_Coulomb_interaction(Vq,epslnq,k,b1,b2,Wq0_j,wi_j,q,tau,Basis_tau_index_cell,WF_center,N_atoms,Cc_kex,Cv,c,cc,v,vv,j,e,eps0) % Modified by Oscar
%
% *** FIX v1.1: Modified to accept Basis_tau_index as a CELL array ***
%               (Changed Basis_tau_index(Ni).data to Basis_tau_index_cell{Ni})
%
[q_mag,MinNdkG] = wave_vector_folding(q,b1,b2); % k-k'+Gkk' in 1st BZ ; q_mag within 1st BZ
[Gi] = Folding_Gvector(b1,b2);
[Cc_kex_tilde,Cv_tilde] = phase_tau_in_Vd(k,Cc_kex,Cv,tau,WF_center);
%-------G dot tau-------------
Gkkp = Gi(MinNdkG,:);
expiGtau = zeros(size(Vq,1),size(Cc_kex_tilde,2));
for Ni=1:N_atoms
    % --- THIS IS THE FIX ---
    % Use curly braces {} to access cell array elements
    expiGtau(:,Basis_tau_index_cell{Ni}) = repmat(exp(1i*sum(bsxfun(@times,[Gkkp(:,1),Gkkp(:,2),zeros(size(k,1),1)],tau(Ni,:)),2)),... % Modified by Oscar
                                              1,length(Basis_tau_index_cell{Ni})); % Modified by Oscar
    % --- END FIX ---
end
%-------G dot tau-------------
Mc = sum((conj(Cc_kex_tilde(:,:,c)).*expiGtau).*Cc_kex_tilde(j,:,cc),2);
Mv = sum((conj(Cv_tilde(:,:,v)).*expiGtau).*Cv_tilde(j,:,vv),2);
Vd = -((e^2)./(4*pi*eps0)).*(Vq./epslnq).*Mc(:).*conj(Mv(:));
Vd(j) = Vd_singularity(e,eps0,Wq0_j,wi_j,Cc_kex_tilde,Cv_tilde,c,cc,v,vv,j);  % Vd(vck,v'c'k)  (the term of k=k')
end
