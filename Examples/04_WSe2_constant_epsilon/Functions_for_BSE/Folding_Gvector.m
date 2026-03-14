function [Gi] = Folding_Gvector(b1,b2)

G0 = [0;0];
G1 = b1;
G2 = b2;
G3 = b1+b2;
G4 = -b1;
G5 = -b2;
G6 = -(b1+b2);
G7 = b1-b2;
G8 = -b1+b2;
Gi = [G0';G1';G2';G3';G4';G5';G6';G7';G8'];

end