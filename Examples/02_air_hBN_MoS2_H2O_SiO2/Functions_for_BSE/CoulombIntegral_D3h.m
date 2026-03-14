function [Vkkp,delta_kV] = CoulombIntegral_D3h(Nt,a,d,q_data,epsQ2D_data)

Vq = @(kx,ky) ((1-((1-exp(-(d.*sqrt(kx.^2+ky.^2))))./(d.*sqrt(kx.^2+ky.^2)))).*((4*pi)./(1.*(kx.^2+ky.^2))));
epsilon_q = @(kx,ky) (epslion_Q2D(q_data, epsQ2D_data, sqrt(kx.^2 + ky.^2)));

v1 = (pi/a).*[2/3,2/sqrt(3)];
v1_unit = v1./Nt;
Nv1_unit = norm(v1_unit);
h = Nv1_unit*(sqrt(3)/2);
l1 = Nv1_unit/sqrt(3);
l2 = Nv1_unit/(2*sqrt(3));

% AbsTol = 1e-7;
% RelTol = 1e-10; % SI unit used
AbsTol = 1e-7;
RelTol = 5*1e-4;
m = -sqrt(3);

f = @(kx, ky) Vq(kx, ky) ./ epsilon_q(kx, ky);

%--------------rectangular---------------
rec_kxmin = -l2*(sqrt(3)/2);
rec_kxmax = l2*(sqrt(3)/2);
rec_kymin = -l2;
rec_kymax = l2/2;

Vkkp_rec_iterated1 = integral2(f,rec_kxmin,0,0,rec_kymax,'method','iterated','AbsTol',AbsTol,'RelTol',RelTol);
Vkkp_rec_iterated2 = integral2(f,0,rec_kxmax,0,rec_kymax,'method','iterated','AbsTol',AbsTol,'RelTol',RelTol);
Vkkp_rec_iterated3 = integral2(f,rec_kxmin,0,rec_kymin,0,'method','iterated','AbsTol',AbsTol,'RelTol',RelTol);
Vkkp_rec_iterated4 = integral2(f,0,rec_kxmax,rec_kymin,0,'method','iterated','AbsTol',AbsTol,'RelTol',RelTol);

Vkkp_rec_iterated = Vkkp_rec_iterated1+Vkkp_rec_iterated2+Vkkp_rec_iterated3+Vkkp_rec_iterated4;
%--------------rectangular---------------

%--------------left trangle---------------
LT_kxmin = -Nv1_unit*cos(pi/3);
LT_kxmax = -l2*(sqrt(3)/2);
LT_kymin = @(kx) -l2+0*kx;
LT_kymax = @(kx) l1+m.*abs(kx);

Vkkp_LT_iterated = integral2(f,LT_kxmin,LT_kxmax,LT_kymin,LT_kymax,'method','iterated','AbsTol',1e-13,'RelTol',1e-10);
%--------------left trangle---------------
%--------------right trangle---------------
RT_kxmin = l2*(sqrt(3)/2);
RT_kxmax = Nv1_unit*cos(pi/3);
RT_kymin = @(kx) -l2+0*kx;
RT_kymax = @(kx) l1+m.*abs(kx);

Vkkp_RT_iterated = integral2(f,RT_kxmin,RT_kxmax,RT_kymin,RT_kymax,'method','iterated','AbsTol',1e-13,'RelTol',1e-10);
%--------------right trangle---------------
%--------------top trangle---------------
TT_kxmin = -l2*(sqrt(3)/2);
TT_kxmax = l2*(sqrt(3)/2);
TT_kymin = @(kx) (l2/2)+0*kx;
TT_kymax = @(kx) l1+m.*abs(kx);

Vkkp_TT_iterated = integral2(f,TT_kxmin,TT_kxmax,TT_kymin,TT_kymax,'method','iterated','AbsTol',1e-13,'RelTol',1e-10);
%--------------top trangle---------------

Vkkp = Vkkp_rec_iterated...
      +Vkkp_LT_iterated+Vkkp_RT_iterated+Vkkp_TT_iterated;

Vkkp = Vkkp/d;

delta_kV = h*Nv1_unit/2; % area of one unit parallelogram


end