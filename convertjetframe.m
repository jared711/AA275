function [Rjet,Ujet] = convertjetframe
%Rjet is Enceladus-fixed position of jets (Cartesian meters)
%Ujet is directon of jets (Cartesian unit vector)

D=importdata('Porco_jetsites_3shift.txt'); D=D.data;
D=D';%rows

%calculations
Az=D(5,:);
Ze=D(6,:);%Ze(:)=0;
lat=D(2,:);
wlon=D(3,:);
U_east = sind(Az).*sind(Ze);
U_north = cosd(Az).*sind(Ze);
U_up = cosd(Ze);

Pole = [0; 0; 1];
rb=[256600 251400 248300];%Enceladus triaxial ellipsoid

Rjet = Tri2Cart(lat,wlon,0,rb);
Up = [cosd(wlon).*cosd(lat); sind(wlon).*cosd(lat); sind(lat)];
East = cross(repmat(Pole,size(lat)), Up); East = East./vecnorm(East);
North = cross(Up, East);
Ujet = U_north.*North + U_east.*East + U_up.*Up;

Rjet(2,:)=-Rjet(2,:);Ujet(2,:)=-Ujet(2,:);%East latitude

function R=Tri2Cart(lat,lon,h,rb)
%Panou, Cartesian to geodetic coordinates conversion on a triaxial ellipsoid using the bisection method
ax=rb(1);ay=rb(2);b=rb(3);
slat=sind(lat);clat=cosd(lat);slon=sind(lon);clon=cosd(lon);
ex2=1-(b/ax).^2;ee2=1-(ay./ax).^2;
N=ax./sqrt(1-ex2.*slat.^2-ee2.*clat.^2.*slon.^2);
R=[(N+h).*clon.*clat ; (N.*(1-ee2)+h).*slon.*clat ; (N.*(1-ex2)+h).*slat];
return