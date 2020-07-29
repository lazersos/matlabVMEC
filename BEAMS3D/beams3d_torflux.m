function flux = beams3d_torflux(beam_data)
%BEASM3D_TORFLUX Extracts the toroidal flux as a function of s.
%   The BEASM3D_TORFLUX function returns a profile of enclosed toroidal
%   flux on the toroidal flux grid.
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       torflux = beams3d_torflux(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


flux=[];

% We assume R/Z/phi equidistant points
dr = beam_data.raxis(2)-beam_data.raxis(1);
dz = beam_data.zaxis(2)-beam_data.zaxis(1);
nphi = beam_data.nphi;

% Area
area=dr.*dz;

% Equilibrium mask
dexs = beam_data.S_ARR <= 1;

% Calc local toroidal flux
flux_local = beam_data.B_PHI;
flux_local(~dexs) = 0;
flux_local = flux_local.*area;

% Setup grid
ns = double(beam_data.ns_prof1);
ds = 1./(ns-1);
edges = 0:ds:1;
s = 0.5.*(edges(1:end-1)+edges(2:end));

% Create profile
flux=zeros(nphi,ns-1);
for i=1:ns-1
    for j=1:nphi
        dex = beam_data.S_ARR(:,j,:) <= edges(i+1);
        temp = flux_local(:,j,:);
        flux(j,i) = sum(temp(dex),'all');
    end
end
flux = mean(flux);

return;

end

