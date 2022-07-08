function [rho, density, vllb, vperpb, pll, pperp, pcross] = beams3d_calc_moments(beam_data)
%BEAMS3D_CALC_MOMENTS Calculates moments of the distribution function
%   The BEAMS3D_CALC_MOMENTS routine calculates the zeroth (density), first
%   (bulk velocity), and second (pressure tensor) moments of the
%   distribution function.  It is asumed the distribution function is in
%   units of s^3/m^6.  
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [rho, density, vllb, vperpb, pll, pperp, pcross] = beams3d_calc_moments(beam_data);
%
%   Note:
%       For old versions with s^3/m^3 try the following and check density
%       against dense_prof:
%           [rhoV,V,dVds] = beams3d_volume(beam_data);
%           dVdrho=dVds.*2.*rho./(beam_data.ns_prof2.*beam_data.ns_prof3.*beam_data.ns_prof1);
%           for i=1:beam_data.ns_prof1, beam_data.dist_prof(:,i,:,:,:,:) = beam_data.dist_prof(:,i,:,:,:,:)./dVdrho(i); end
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


% Get volume density
[rhoV,~,dVds] = beams3d_volume(beam_data);
dVdrho = 2.*dVds.*rhoV;
dVdrho = repmat(dVdrho,[beam_data.nbeams,1]);

% Get mass arrary
mass = ones(beam_data.nbeams,beam_data.ns_prof1);
for i=1:beam_data.nbeams
    j=find(beam_data.Beam==i,1,'first');
    mass(i,:) = beam_data.mass(j);
end

% Rho
x=0:1./beam_data.ns_prof1:1;
rho = 0.5.*(x(1:end-1)+x(2:end));

% Distribution
dist = beam_data.dist_prof;

% axes
x=-1:2./beam_data.ns_prof4:1;
vll = 0.5.*(x(1:end-1)+x(2:end)).*beam_data.partvmax;
x=0:1./beam_data.ns_prof5:1;
vperp = 0.5.*(x(1:end-1)+x(2:end)).*beam_data.partvmax;
dV = diff(vll(1:2)).*diff(vperp(1:2)).*pi.*2;
dthdphi = 1./(beam_data.ns_prof2.*beam_data.ns_prof3);

d3V = ones([beam_data.nbeams beam_data.ns_prof1 beam_data.ns_prof2 beam_data.ns_prof3 beam_data.ns_prof4 beam_data.ns_prof5]);

for i=1:beam_data.ns_prof5
    d3V(:,:,:,:,:,i) = d3V(:,:,:,:,:,i).*dV.*vperp(i);
end

density = squeeze(sum(dist.*d3V.*dthdphi,[3 4 5 6])); % Close to actual value [part/m^3]

% vllb
temp = ones([beam_data.nbeams beam_data.ns_prof1 beam_data.ns_prof2 beam_data.ns_prof3 beam_data.ns_prof4 beam_data.ns_prof5]);
for i=1:beam_data.ns_prof4
    temp(:,:,:,:,i,:) = dist(:,:,:,:,i,:).*vll(i);
end
vllb = squeeze(sum(temp.*d3V.*dthdphi,[3 4 5 6]))./density;

% vperpb
temp = ones([beam_data.nbeams beam_data.ns_prof1 beam_data.ns_prof2 beam_data.ns_prof3 beam_data.ns_prof4 beam_data.ns_prof5]);
for i=1:beam_data.ns_prof5
    temp(:,:,:,:,:,i) = dist(:,:,:,:,:,i).*vperp(i);
end
vperpb = squeeze(sum(temp.*d3V.*dthdphi,[3 4 5 6]))./density;

% pll
wll = ones([beam_data.nbeams beam_data.ns_prof1 beam_data.ns_prof2 beam_data.ns_prof3 beam_data.ns_prof4 beam_data.ns_prof5]);
for i=1:beam_data.nbeams
    for j=1:beam_data.ns_prof1
        wll(i,j,:,:,:,:) = -vllb(i,j);
    end
end
for i=1:beam_data.ns_prof4
    wll(:,:,:,:,i,:) = wll(:,:,:,:,i,:)+vll(i);
end
w2 = wll.*wll;
pll = squeeze(sum(w2.*dist.*d3V.*dthdphi,[3 4 5 6])).*mass;

% pperp
wperp = ones([beam_data.nbeams beam_data.ns_prof1 beam_data.ns_prof2 beam_data.ns_prof3 beam_data.ns_prof4 beam_data.ns_prof5]);
for i=1:beam_data.nbeams
    for j=1:beam_data.ns_prof1
        wperp(i,j,:,:,:,:) = -vperpb(i,j);
    end
end
for i=1:beam_data.ns_prof5
    wperp(:,:,:,:,:,i) = wperp(:,:,:,:,:,i)+vperp(i);
end
w2 = wperp.*wperp;
pperp = squeeze(sum(w2.*dist.*d3V.*dthdphi,[3 4 5 6])).*mass;
pcross = squeeze(sum(wll.*wperp.*dist.*d3V.*dthdphi,[3 4 5 6])).*mass;


end