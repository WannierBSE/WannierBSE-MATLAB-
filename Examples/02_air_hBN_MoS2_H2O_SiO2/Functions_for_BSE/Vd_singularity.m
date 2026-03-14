function [Vd_kk] = Vd_singularity(e,eps0,Wq0,wi,Cc_kex_tilde,Cv_tilde,c,cc,v,vv,i)

Vd_kk = -(((e^2)./(4*pi*eps0)).*(Wq0)...
        .*sum(conj(Cc_kex_tilde(i,:,c)).*Cc_kex_tilde(i,:,cc))...
        .*sum(Cv_tilde(i,:,v).*conj(Cv_tilde(i,:,vv)))); 

Vd_kk = Vd_kk/wi;

end