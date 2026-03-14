function [Wq0,wi] = Singularity_and_kweight(N_de,N_co,Nv_coarse,Nv_dense,d,...
                                       q_data,epsQ2D_data,a,f)
                                   
[Wq0_co,wi_co] = CoulombIntegral_D3h(Nv_coarse,a,d,q_data,epsQ2D_data);
[Wq0_de,wi_de] = CoulombIntegral_D3h_dense(Nv_dense,a,d,q_data,epsQ2D_data,f);

wi = [wi_de*ones(N_de,1);wi_co*ones(N_co,1)];
Wq0 = [Wq0_de*ones(N_de,1);Wq0_co*ones(N_co,1)];
end