function [HSP] = High_Symmetry_Points(b1,b2)


for i=1:6

if i==1
bn = -b2;
bm = b1;
elseif i==2
bn = b1;
bm = b1+b2;
elseif i==3
bn = b1+b2;
bm = b2;
elseif i==4
bn = b2;
bm = -b1;
elseif i==5
bn = -b1;
bm = -(b1+b2);
elseif i==6
bn = -(b1+b2);
bm = -b2;
end

bnx = bn(1);
bny = bn(2);
bmx = bm(1);
bmy = bm(2);

kx = (((norm(bn)^2)*bmy)-((norm(bm)^2)*bny))/(2*(bnx*bmy-bmx*bny));
ky = ((norm(bn)^2)/(2*bny))-((bnx/bny)*kx);
Mx = bmx./2;
My = bmy./2;


HSP((1:2)+(i-1)*2,:) = [kx,ky,0;Mx,My,0];
end

HSP = [0,0,0;HSP];

end