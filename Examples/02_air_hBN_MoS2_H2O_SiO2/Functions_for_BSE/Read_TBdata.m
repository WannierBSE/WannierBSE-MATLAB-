function [En,Cn] = Read_TBdata(Nk_all,nB_data,e,NB)


En = zeros(Nk_all,1);
Cn = zeros(NB,Nk_all);
for i = 1:Nk_all
    En(i) = nB_data(2+(i-1)*(NB+2),1).*e;
    Cn(1:NB,i) = nB_data((3+(i-1)*(NB+2)):((3+NB-1)+(i-1)*(NB+2)),1)...
                +1i*nB_data((3+(i-1)*(NB+2)):((3+NB-1)+(i-1)*(NB+2)),2);
end
Cn = transpose(Cn);


end