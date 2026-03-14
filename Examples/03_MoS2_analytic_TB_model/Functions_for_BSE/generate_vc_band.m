function [vc_band_BSE] = generate_vc_band(Nv,Nc)

vc_band_BSE = zeros(Nc*Nv,2);
for vi=1:Nv
for ci=1:Nc
    vc_band_BSE(ci+(vi-1)*Nc,:) = [vi,ci];
end
end

end