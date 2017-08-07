% function [V,CoM,I] = volumeintegrate(vertices, triangle_faces)
%
% polyhedron volume, center of mass and inertia matrix
%
% polyhedron given in standard MATLAB face-vertex form
% must have triangular faces
%
% single argument is assumed to be handle of a patch object
%
% this function calculates the volume from an assumed closed surface mesh defined by
% Triangles and Positions (Nodepositions). The Triangles list must have a counterclockwise numbering
% to ensure the outward normal.
%
% /*******************************************************
% from volInt.c by Brian Mirtich
% Ref "Fast and Accurate Computation of Polyhedral Mass Properties"
% Journal of Graphics Tools 1(1) 1996.
%
% Ported to Matlab by Andrew Gosline; www.ece.ubc.ca/~andrewg/ *
%
% Simplified to give mass, centre of mass and inertia matrix
% from face-vertex data
% By Kim McKelvey, University of Otago, 2005
%
% * Copyright 1995 by Brian Mirtich *
% * mirtich@cs.berkeley.edu
% * http://www.cs.berkeley.edu/~mirtich
function [V,CoM,I] = VolumeIntegrate(Xn,Triangles);
if nargin==1
% Normals = get(Xn, 'vertexnormals');
Triangles = get(Xn, 'faces');
Xn = get(Xn, 'vertices');
if size(Triangles,2)~=3
error('faces must be triangles')
return
end
end
% Compute the Outer Normals of the Surface
norm = cross(Xn(Triangles(:,2),:)-Xn(Triangles(:,1),:),Xn(Triangles(:,3),:)-Xn(Triangles(:,1),:));
% len = 1./sqrt(diag(norm*norm'));
warning off
len = 1./sqrt(sum((norm.*norm)'))';
warning on
Normals = norm.*[len,len,len];
T0 = 0; Tx = 0; Ty = 0; Tz = 0; Txx = 0; Tyy = 0; Tzz = 0; Txy = 0; Tyz = 0; Tzx = 0;
T1 = zeros(3,1); T2 = zeros(3,1); TP = zeros(3,1);
for i = 1:length(Triangles)
nx = abs(norm(i,1));
ny = abs(norm(i,2));
nz = abs(norm(i,3));
if (nx==0)&(ny==0)&(nz==0)
else
if (nx > ny) & (nx > nz)
C = 0;
elseif (ny > nz)
C = 1;
else
C = 2;
end
% Cycle the indices so always in numerical order
A = mod(C+1,3);
B = mod(A+1,3);
A = A+1; B = B+1;, C = C+1;
% calculate the offset
w = -Normals(i,1)*Xn(Triangles(i,1),1) - Normals(i,2)*Xn(Triangles(i,1),2) - Normals(i,3)*Xn(Triangles(i,1),3);
% This part is essentially just the ComputeProjectionIntegrals part
Pa = 0; Pb = 0; P1 = 0; Paa = 0; Pab = 0; Pbb = 0; Paaa = 0; Paab = 0; Pabb = 0; Pbbb = 0;
for j = 1:3
a0 = Xn(Triangles(i,j),A);
b0 = Xn(Triangles(i,j),B);
a1 = Xn(Triangles(i, mod(j,3)+1),A);
b1 = Xn(Triangles(i, mod(j,3)+1),B);
da = a1 - a0;
db = b1 - b0;
a0_2 = a0 * a0;
a0_3 = a0_2 * a0;
a0_4 = a0_3 * a0;
b0_2 = b0 * b0;
b0_3 = b0_2 * b0;
b0_4 = b0_3 * b0;
a1_2 = a1 * a1;
a1_3 = a1_2 * a1;
b1_2 = b1 * b1;
b1_3 = b1_2 * b1;
C1 = a1 + a0;
Ca = a1*C1 + a0_2;
Caa = a1 * Ca + a0_3;
Caaa = a1*Caa + a0_4;
Cb = b1 * (b1 + b0) + b0_2;
Cbb = b1*Cb + b0_3;
Cbbb = b1*Cbb + b0_4;
Cab = 3*a1_2 + 2*a1*a0 + a0_2;
Kab = a1_2 +2*a1*a0 + 3*a0_2;
Caab = a0*Cab + 4*a1_3;
Kaab = a1*Kab + 4*a0_3;
Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;
P1 = P1 + (db * C1);
Pa = Pa + (db * Ca);
Paa = Paa + db*Caa;
Paaa = Paaa + db*Caaa;
Pb = Pb + (da * Cb);
Pbb = Pbb + da*Cbb;
Pbbb = Pbbb + da*Cbbb;
Pab = Pab + db*(b1*Cab + b0*Kab);
Paab = Paab + db*(b1*Caab + b0*Kaab);
Pabb = Pabb + da*(a1*Cabb + a0*Kabb);
end
P1 = P1/2.0;
Pa = Pa/6.0;
Paa = Paa/12.0;
Paaa = Paaa/20.0;
Pb = Pb/-6.0;
Pbb = Pbb/-12.0;
Pbbb = Pbbb/-20.0;
Pab = Pab/24.0;
Paab = Paab/60.0;
Pabb = Pabb/-60.0;
% disp(sprintf('%e %e %e %e %e %e %e %e %e %e \n',P1,Pa,Paa,Paaa,Pb,Pbb,Pbbb,Pab,Paab,Pabb));
% This Part is just the Compute Face Integrals part
k1 = 1 / Normals(i,C);
k2 = k1 * k1;
k3 = k2 * k1;
k4 = k3 * k1;
Fa = k1 * Pa;
Fb = k1 * Pb;
Fc = -k2 * (Normals(i,A)*Pa + Normals(i,B)*Pb + w*P1);
Faa = k1 * Paa;
Fbb = k1 * Pbb;
Fcc = k3 * ((Normals(i,A))^2*Paa + 2*Normals(i,A)*Normals(i,B)*Pab + (Normals(i,B))^2*Pbb + ...
w*(2*(Normals(i,A)*Pa + Normals(i,B)*Pb) + w*P1));
Faaa = k1 * Paaa;
Fbbb = k1 * Pbbb;
Fccc = -k4 * (Normals(i,A)^3*Paaa + 3*Normals(i,A)^2 * Normals(i,B)*Paab ...
+ 3*Normals(i,A)*(Normals(i,B))^2*Pabb + (Normals(i,B))^3*Pbbb ...
+ 3*w*((Normals(i,A))^2*Paa + 2*Normals(i,A)*Normals(i,B)*Pab + (Normals(i,B))^2*Pbb) ...
+ w*w*(3*(Normals(i,A)*Pa + Normals(i,B)*Pb) + w*P1));
Faab = k1 * Paab;
Fbbc = -k2 * (Normals(i,A)*Pabb + Normals(i,B)*Pbbb + w*Pbb);
Fcca = k3 * ((Normals(i,A))^2*Paaa + 2*Normals(i,A)*Normals(i,B)*Paab + (Normals(i,B))^2*Pabb ...
+ w*(2*(Normals(i,A)*Paa + Normals(i,B)*Pab) + w*Pa));
if A == 1
Part = Fa;
elseif B == 1
Part = Fb;
else
Part = Fc;
end
T0 = T0 + Normals(i,1)*Part;
T1(A) = T1(A) + Normals(i,A) * Faa;
T1(B) = T1(B) + Normals(i,B) * Fbb;
T1(C) = T1(C) + Normals(i,C) * Fcc;
T2(A) = T2(A) + Normals(i,A) * Faaa;
T2(B) = T2(B) + Normals(i,B) * Fbbb;
T2(C) = T2(C) + Normals(i,C) * Fccc;
TP(A) = TP(A) + Normals(i,A) * Faab;
TP(B) = TP(B) + Normals(i,B) * Fbbc;
TP(C) = TP(C) + Normals(i,C) * Fcca;
end
end
T1 = T1./2;
T2 = T2./3;
TP = TP./2;
% compute the desired quantities
V = T0;
CoM = T1/T0;
Iqq = [T2(2)+T2(3)-V*(CoM(2)^2+CoM(3)^2),T2(1)+T2(3)-V*(CoM(1)^2+CoM(3)^2),T2(1)+T2(2)-V*(CoM(2)^2+CoM(1)^2)];
Iqz = [TP(1)-V*CoM(1)*CoM(2);TP(2)-V*CoM(2)*CoM(3);TP(3)-V*CoM(1)*CoM(3)];
I = [Iqq(1),-Iqz(1),-Iqz(3);-Iqz(1),Iqq(2),-Iqz(2);-Iqz(3),-Iqz(2),Iqq(3)];
% volume and inertia are negative
% if faces were labelled the wrong way round
if V<0
V=-V;
I=-I;
end