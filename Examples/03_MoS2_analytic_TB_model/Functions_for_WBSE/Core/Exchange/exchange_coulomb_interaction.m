function Vx = exchange_coulomb_interaction(chi_p,R_mNN,k,N_kpt,Cv,Cc,c,cc,v,vv,j)
%==========================================================================
% EXCHANGE_COULOMB_INTERACTION
%
% DESCRIPTION:
%   Contracts the precomputed short-range exchange tensor with TB coefficients to form the exchange BSE matrix contribution.
%
% OUTPUT:
%   Vx : exchange Coulomb contribution over k-points.
%==========================================================================


kex = [0, 0, 0];
Vx = zeros(N_kpt,1);

for R_ind1=1:size(R_mNN,1)
for R_ind2=1:size(R_mNN,1)
    chi_R1R2 = chi_p(:,:,:,:,R_ind1,R_ind2);
    
    S1 = sum((transpose(conj(Cv(j,:,vv))).*chi_R1R2),1);
    Cc_RR = transpose(Cc(j,:,cc).*exp(1i*sum((k(j,:)+kex(1:2)).*(R_mNN(R_ind2,1:2)),2)));
    S2 = sum(Cc_RR.*(squeeze(S1)),1);
    S2 = squeeze(S2);
    
    Cc_R = Cc(:,:,c).*exp(1i*sum((k+kex(1:2)).*(R_mNN(R_ind1,1:2)),2));
    Vx = Vx+sum((Cv(:,:,v)*S2).*conj(Cc_R),2);
        
end
end

end
