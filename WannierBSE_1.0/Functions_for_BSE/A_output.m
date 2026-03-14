function A_output = A_output(A,Nc,Nv,wi,N_kpt)

%A_OUTPUT:Converts A from unitless to angstroms (A)


scale = sqrt((4*pi^2)./wi);                 % N_kpt x 1
scale_all = kron(ones(Nc*Nv,1), scale);     % (Nc*Nv*N_kpt) x 1
A_output = A .* scale_all;                 

end    



