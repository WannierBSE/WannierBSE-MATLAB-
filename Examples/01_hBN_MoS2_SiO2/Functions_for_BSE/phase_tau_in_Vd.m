function [Cc_tilde,Cv_tilde] = phase_tau_in_Vd(k,Cc,Cv,tau,WF_center)

NB = size(WF_center,1);
expiktau = zeros(size(k,1),NB);
for i=1:NB
    tau_i = tau(WF_center(i,2),:);
    expiktau(:,i) = exp(-1i*sum(bsxfun(@times,[k(:,1),k(:,2),zeros(size(k,1),1)],tau_i),2));
end


Cv_tilde = zeros(size(Cv));
for n=1:size(Cv,3)
for i=1:size(Cv,2)
    Cv_tilde(:,i,n) = Cv(:,i,n).*expiktau(:,i);
end
end

Cc_tilde = zeros(size(Cc));
for n=1:size(Cc,3)
for i=1:size(Cc,2)
    Cc_tilde(:,i,n) = Cc(:,i,n).*expiktau(:,i);
end
end

end