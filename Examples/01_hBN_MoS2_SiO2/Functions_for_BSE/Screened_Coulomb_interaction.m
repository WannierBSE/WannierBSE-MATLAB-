function [Vq,epsilon_q] = Screened_Coulomb_interaction(q,b1,b2,d,q_data,epsQ2D_data)

[q_mag,MinNdkG] = wave_vector_folding(q,b1,b2); % k-k'+Gkk' in 1st BZ ; q_mag within 1st BZ

epsilon_q = epslion_Q2D(q_data,epsQ2D_data,q_mag);
Vq = (1-((1-exp(-(d.*q_mag)))./(d.*q_mag))).*((4*pi)./(d.*(q_mag.^2)));
end