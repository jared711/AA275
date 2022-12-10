function [flux, density, Vave, diam] = Enceladus_jet(Rsc,Vsc,Usc,Rjet,Ujet,up_down,param_set,exit_distribution)
%calculates flux for input s/c trajectory, collector pointing, and jet pointing
%[flux, density, Vave] = Enceladus_jet(Rsc, Vsc, Usc, Rjet, Ujet, up_down)
%flux = g/s/m^2 into collector
%density = g/m^3 of jet at Rsc
%Vave = m/s of jet at Rsc averaged over density
%
%Rsc = position of s/c with respect to jet, m Cartesian 3 x n array for n positions
%Vsc = velocity of s/c wrt jet, m/s Cartesian, default = [0; 0; 0]
%Usc = direction of collector opening, unit vector, default = [0; 0; -1]
%Rjet = position of jet wrt Enceladus center, m Cartesian, default = [0; 0; 248300]
%Ujet = direction of jet opening, unit vector, default = [0; 0; 1]
%up_down = +1 to intercept particles on the way up or -1 on the way down, default = +1.
%
%assumes ballistic particle motion subject to Enceladus gravity
%does not account for Saturn gravity, electromagnetic forces, motion of jet source
%
%EXAMPLES
%hover at 400 m altitude, 50 m off jet center
% Rsc = [0; 50; 400];
% flux = Enceladus_jet(Rsc)
%
%s/c flyby at 3 km alt & 200 m/s, s/c pointed 40 deg off nadir, jet pointed 5 deg off zenith
% time = -6 : 0.1 : +6; % seconds
% Vsc = [200; 0; 0]; %m/s, s/c velocity
% Rsc = [0; 0; 3000] + Vsc*time; %m, s/c position
% Usc = [sind(40); 0; -cosd(40)]; %unit vector
% Ujet = [0; sind(5); cosd(5)]; %unit vector
% flux = Enceladus_jet(Rsc, Vsc, Usc, [], Ujet);
% plot(time, flux)

%brought to you in part by:
%Southworth et al. Surface deposition of the Enceladus plume and the zenith angle of emissions. https://doi.org/10.1016/j.icarus.2018.08.024
%Guzman et al. Collecting amino acids in the Enceladus plume. https://doi.org/10.1017/S1473550417000544

%Damon Landau (damon.landau@jpl.nasa.gov)
%4/27/22 v0
%6/23/22 v1 added "int_diam" input; allow up_down to be input as [+1 -1] to calc both in single call
%10/7/22 v3 updated parameters based on Cassini data


%% default inputs
if nargin<2||isempty(Vsc); Vsc = [0;0;0]; end
if nargin<3||isempty(Usc); Usc = [0;0;-1]; end
if nargin<4||isempty(Rjet); Rjet = [0; 0; 248.3e3]; end%m, radius at pole of Enceladus
if nargin<5||isempty(Ujet); Ujet = [0; 0; 1]; end
if nargin<6||isempty(up_down); up_down = +1;end
Rsc = Rsc + Rjet;

%% super-user inputs
if nargin<7||isempty(param_set); param_set='HRD2';end%model parameters based on Cassini data
if nargin<8||isempty(exit_distribution); exit_distribution='cos2_angle';end%options for amount of jet collimation.
calc_mass = true;%false to output particle count instead of mass
int_diam = true;%false to output array at individual particle diameters, true to integrate output across diameters
pow_law=true;%except for HRD2
d_min = 0.2; d_max = 60;%um, range of particle diameters. Southworth uses [1 30], Guzman uses [0.02 20]
switch param_set
case 'custom'
    m_total=1e3;%g/s, single jet production. Southworth uses 25e3/80, Guzman uses 16e3/100 (total plume/# of jets)
    %particle size distribution
    alpha = 3.1;%power law probability, p(d) = c_alpha*d^-alpha. Southworth uses 3.1, Guzman uses 4.0
    %exit speed distribution
    vgas = 700;%m/s, max exit speed, Southworth uses 700, Guzman uses 350
    dc = 1.6;%um, critical diameter, see speed distribution below 
case 'HRD1'
    m_total = 400e3/100;%g/s
    alpha = 4.4;
    vgas = 1500;%m/s
    dc = 3.2;%um
case 'HRD2'
    m_total = 90e3/100;%g/s
    mu = -1.2; sigma = 1;
    pow_law=false;        
    vgas = 1500;%m/s
    dc = 3.2;%um
case 'VIMS'
    m_total = 141e3/100;%g/s
    alpha = 2.3;
    vgas = 1000;%m/s
    dc = 0.3;%um
end
%exit angle distribution
max_exit_angle=15;%deg, half-angle of exit cone. Southworth uses 15, Guzman uses 15

Delta_d = 0.1;%um, resolution of diameter for numerical integration
Delta_v = 2;%m/s, resolution of exit speeds for numerical integration
d_range = [d_min d_max];%size range to calculate flux
v_range = [0 vgas];%speed range to calculate flux
Delta_d = min(Delta_d, diff(d_range));
Delta_v = min(Delta_v, diff(v_range));

%% average particle mass, particle sizes, exit speeds, exit angle distribution
if pow_law%power law distribution of particle size
    c_alpha = (alpha-1)/(d_min^(1-alpha)-d_max^(1-alpha));%normalize size distribution to unity
%     if alpha~=4; m_ave=(d_min^(4-alpha)-d_max^(4-alpha))/(alpha-4); else; m_ave=log(d_max)-log(d_min); end
    p_diam=@(d) c_alpha*d.^-alpha;
else%lognormal 
    c_alpha=2/( erf((log(d_max/2)-mu)/sigma/sqrt(2)) - erf((log(d_min/2)-mu)/sigma/sqrt(2)) );
    p_diam=@(d) c_alpha/sigma/sqrt(2*pi)./d.*exp(-(log(d/2)-mu).^2/(2*sigma^2));
end
c_mass = 0.916e-12*(4/3*pi/8);%g/um^3, particle mass = c_mass*d^3
diam = d_range(1) : Delta_d : d_range(2); diam= diam(:);%sample particle diameters
m_ave=c_mass*trapz(diam,p_diam(diam).*diam.^3);%g, average particle mass

%exit angle distribution
emax=max_exit_angle*pi/180;
switch exit_distribution%distribution of jet exit direction, fraction of total per steradian (integrate over 2-D exit directions to get unity)
    case 'cos2_angle'; mass_exit=@(q) (q<emax)./(pi*emax*sin(q)).*cos(pi/2/emax*q).^2;%cos^2 exit radian, most collimated
    case 'uni_angle'; mass_exit=@(q) (q<emax)./(emax*2*pi*sin(q));%uniform exit radian
    case 'cos2_area'; n=2;q_=linspace(0,emax,1e4);c_=trapz(q_,2*pi*sin(q_).*cos(pi/2/emax*q_).^n);
                      mass_exit=@(q) 1/c_*(q<emax).*cos(pi/2/emax*q).^n;%cos^n exit steradian
    case 'uni_area'; mass_exit=@(q) (q<emax)./((1-cos(emax))*2*pi);%uniform exit steradian, least collimated
end

%exit speed sampling
gm = 7.2105e9;%m^3/s^2, gravitational parameter of Enceladus
rjet = vecnorm(Rjet);
vesc = sqrt(2*gm./rjet);%m/s, conic escape speed
sma=(vecnorm(Rsc)+rjet+vecnorm(Rsc-Rjet))/4; vmin = sqrt(gm*(2./rjet-1./sma));%minimum speed to reach Rsc (see conic_BVP.m)
vesc = min(vesc, v_range(2)); v_range(1) = max(v_range(1), min(vmin) );%adjust bounds if necessary
if up_down==-1; v_range(2) = vesc; end%omit degenerate (negative time) cases for hyperbolas
speed_bins = [v_range(1) : Delta_v : vesc, v_range(2)];%sample different exit speeds, lump all escaped in single bin
speeds = (speed_bins(1:end-1)+speed_bins(2:end))/2;%approximate speed in middle of bins

%% calculate flux across all speed bins
flux = 0; density = 0; Vave = 0;
for sbi=1:numel(speeds)
    %% exit speed and mass distributions
    % p=(1+diam/dc).*(diam/dc).*speeds/vgas^2.*(1-speeds/vgas).^(diam/dc-1);%speed distribution
    v1=speed_bins(sbi); v2=speed_bins(sbi+1); if v1==v2;continue;end
    p = (1+diam/dc*v1/vgas).*(1 - v1/vgas).^(diam/dc) - (1+diam/dc*v2/vgas).*(1 - v2/vgas).^(diam/dc);%analytic integration across speed bin
    p = p_diam(diam).*p;%particle size distribution

    if calc_mass
        pm = c_mass/m_ave*diam.^3.*p;%particle mass distribution
        pm = .5*diff(diam).*(pm(1:end-1,:)+pm(2:end,:));%trapezoidal weighting
        if int_diam; pm = sum( pm, 1); end%integration across masses
        flux_bin = m_total*pm;%total mass within individual size/speed bin
    else
        p = .5*diff(diam).*(p(1:end-1,:)+p(2:end,:));%trapezoidal weighting
        if int_diam; p=sum( p, 1); end% integration across sizes
        flux_bin = m_total/m_ave*p;%total particles within individual size/speed bin
    end

    %% particle mass and velocity at different s/c positions
    vjet = speeds(sbi);%jet exit speed
    imin = v1<vmin & vmin<v2;%check if min possible speed is within bin, omit portion < vmin
    if any(imin)%reset speed bin with v1 = vmin
        vjet = repmat(vjet,size(vmin)); flux_bin = repmat(flux_bin,size(vmin));%make array
        vjet(imin) = (vmin(imin)+v2)/2;%vjet in middle of smaller bin 
        flux_bin(:,imin) = (v2-vmin(imin))/(v2-v1) .* flux_bin(:,imin);%decrement flux for smaller bin
    end
    sma=2*gm./rjet-vjet.^2; sma=gm./sma;%semi-major axis of particle orbit
    for ud_=up_down
        if ud_==-1 && speeds(sbi)>vesc; continue; end%omit degenerate (negative time) cases for hyperbolas
        %Vjet is jet exit velocity, Vint is particle velocity at s/c intercept (both inertial), Vj_Rs is partial of Vjet wrt Rsc
        [Vjet, Vint, Vj_Rs] = conic_BVP(Rjet, Rsc, sma, ud_, gm);%calc conic orbit connecting Rjet to Rsc
        ejet=sum(Vjet.*Ujet)./vjet; ejet(ejet>1)=1; ejet=acos(ejet);%exit angle from jet center (check roundoff cos(q)>1)
        mass_Vj = flux_bin.*mass_exit(ejet);%g/s/steradian in exit direction
        if ~any(mass_Vj);continue;end
        %% calculate mass in 1 m^3 box (density) for flux at Rsc
        %orbital coordinate frame along velocity and angular momentum at intercept
        vint = sqrt(sum(Vint.^2));
        H=cross(Rsc,Vint); h=sqrt(sum(H.^2)); H=H./h;%unit vector along angular momentum (normal to orbit plane)
        T=cross(H,Vint)./vint;%unit vector to complete orthogonal frame with Vint and H
        %correct for rectilinear orbit when Rsc and Rjet are aligned and plane is undefined
        for rect = find(h==0); HT = null(Vint(:,rect)'); H(:,rect)=HT(:,1); T(:,rect)=HT(:,2); end
        Vj_H = sum(Vj_Rs.*permute(H, [3 2 1]), 3);%change in Vjet for 1-m change in Rsc along H (matrix multiply along dim 3)
        Vj_T = sum(Vj_Rs.*permute(T, [3 2 1]), 3);%change in Vjet for 1-m change in Rsc along T
        %for conics out- and in-plane motion are decoupled so Vj_H & Vj_T remain orthogonal because H & T are orhogonal
        %vjet = |Vjet| is fixed, so Vj_H & Vj_T remain orthogonal to Vjet
        %A_Vj = cross(Vj_H./vjet, Vj_T./vjet); A_Vj = dot(A_Vj, Vjet./vjet);
        A_Vj = sqrt(sum(Vj_H.^2).*sum(Vj_T.^2))./vjet.^2;%steradian of jet exit direction to cover 1 m^2 normal to Vint

        time_V = 1./vint;%change in time for 1-m change in Rsc along Vint, duration jet produced mass to fill box
        density_v = mass_Vj.*time_V.*A_Vj;%g/m^3 at Rsc
        density_v = max( density_v, 0);
        if int_diam; Vave = Vave + density_v.*Vint; end

        Vrel = Vsc - Vint;%relative velocity at Rsc
        flux_v = density_v.*sum(Vrel.*Usc);%g/s/m^2, flux into Usc 
        flux = flux + max( flux_v, 0); density = density + density_v;%correct for negative or NaN (when vjet insufficient to reach Rsc)
    end
end%sbi
Vave = Vave./density;
diam = (diam(1:end-1)+diam(2:end))/2;%midpoint of diameter bins

return
