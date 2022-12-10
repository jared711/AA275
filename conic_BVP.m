function [V1,V2,V1_R2,V2_R2]=conic_BVP(R1,R2,a,sx,gm,sq)
%conic_BVP solves orbital Boundary Value Problem for given semi-major axis
%[V1, V2, V1_R2] = conic_BVP(R1, R2, a, s_x, gm, s_q)
%R1, R2 are initial and final positions, 3 x n Cartesian, length units LU
%a is semi-major axis of orbit, 1 x n array, LU
%s_x is +1 or -1 to specify solution with more or less angular momentum
%gm is gravitational parameter LU^3/TU^2, time units TU
%s_q is +1 or -1 for transfer angle <180 or >180 deg
%
%V1, V2 are initial and final velocities, 3 x n Cartesian, LU/TU
%V1_R2 is derivative of V1 wrt R2, 3 x n x 3, 1/TU along dimension 3
%
%there is no solution if a < a_min = (|R1| + |R2| + |R2-R1|)/4
%otherwise there are 4 solutions through permutations of s_x and s_q
%two of these solutions are time-reverse of the other pair: (V1,V2) = -(V1,V2)

%Gooding, A Procedure for the Solution of Lambert's Orbital Boundary-Value Problem
%Battin, An Introduction to the Mathematics and Methods of Astrodynamics pp. 179, 186, 242

%could improve roundoff with input DR instead of R2 for very small DR

if nargin < 6;sq=+1;end
%geometry
r1=sqrt(sum(R1.^2));r2=sqrt(sum(R2.^2));r12=sum(R1.*R2);
DR=R2-R1;c=sqrt(sum(DR.^2));%c = sqrt(r1.^2+r2.^2-2*r12);
s=(r1+r2+c)/2;
cs=c./s;%cs=1-q.^2
q=r1.*r2+r12;q=sq.*sqrt(q/2)./s;%q < 0 for transfer angle > 180 deg 

%energy and angular momentum from Gooding
x=1-1/2*s./a; x(x<0)=nan;
x=sx.*sqrt(x);
qx=q.*x; z=sqrt(cs+qx.^2);
zqx=z+qx;
%h = sqrt(gm*s.*(r1.*r2-r12))./c.*(z+qx);%angular momentum
%velocities from Battin
cv=c.*cs./zqx.^2;
vv=sqrt(gm/2./s)./q.*zqx./c;%h = vv.*vecnorm(cross(R1,R2));
%cv*vv=sqrt(gm/2./s)./q.*cs./zqx;
V1 = DR + cv./r1.*R1; V1 = vv.*V1;
V2 = DR - cv./r2.*R2; V2 = vv.*V2;

if nargout>2%derivatives
%derivativs wrt [r1; r2; r12]
dc=-ones(3,numel(r12)); dc(1,:)=r1; dc(2,:)=r2; dc = dc./c;
ds=[1;1;0]+dc;ds=ds/2;
dcs=(dc-cs.*ds)./s;
dq=-1/2./q.*dcs;

dx=-1/4./x./a.*ds;
dqx=x.*dq+q.*dx;
dz=1/2./z.*dcs+qx./z.*dqx;
dzqx=dz+dqx;
dc_=dc./c-dzqx./zqx; ds_=1/2./s.*ds;
dcv=2*cv.*(dc_-ds_);
dvv=-dc_-ds_-dq./q;
%write derivative wrt Ri along dimension 3
pR1=permute(R1,[3 2 1]); pR2=permute(R2,[3 2 1]); 
vvI = zeros(3,numel(vv),3);for ii=1:3;vvI(ii,:,ii)=vv;end
% %dV1_dR1
% V1_R1 = (dvv(1,:)./r1.*pR1+dvv(3,:).*pR2).*V1 ...
%       + vvI.*(-1+cv./r1) ...
%       + vv.*((dcv(1,:)-cv./r1)./r1./r1.*pR1+dcv(3,:)./r1.*pR2).*R1;
%dV2_dR2
V2_R2 = (dvv(2,:)./r2.*pR2+dvv(3,:).*pR1).*V2 ...
      + vvI.*(1-cv./r2) ...
      - vv.*((dcv(2,:)-cv./r2)./r2./r2.*pR2+dcv(3,:)./r2.*pR1).*R2;
%dV1_dR2
V1_R2 = (dvv(2,:)./r2.*pR2+dvv(3,:).*pR1).*V1 ...
      + vvI ...
      + vv.*(dcv(2,:)./r1./r2.*pR2+dcv(3,:)./r1.*pR1).*R1;
end
