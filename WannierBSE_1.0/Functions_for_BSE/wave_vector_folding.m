function [q_mag,MinNdkG] = wave_vector_folding(q,b1,b2)

    [Gi] = Folding_Gvector(b1,b2);

    qG = zeros(size(q,1),2,size(Gi,1));
    for gi=1:size(Gi,1)
    qG(:,:,gi) = q+Gi(gi,:);
    end
    NqG = squeeze(sqrt(sum(qG.^2,2)));
    [q_mag,MinNdkG] = min(NqG,[],2);

end