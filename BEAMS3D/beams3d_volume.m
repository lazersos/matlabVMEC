function [s, plasma_vol, plasma_dvolds] = beams3d_volume(beam_data,varargin)
%BEASM3D_VOLUME Extracts plasma volume data
%   The BEASM3D_VOLUME function calculates volume and differential volume
%   (dV/ds) from the background grid data.  Values are returned on the half
%   grid in a similar fashion to the distribution functions.
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [s, V, dVds] = beams3d_volume(beam_data);
%
%       [s, V, dVds] = beams3d_voluem(beam_data,'ns',128);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0

ns = [];

% Handle varargin
if ~isempty(varargin)
    i = 1;
    while i <= length(varargin)
        switch varargin{i}
            case 'ns'
                ns = varargin{i+1};
        end
        i=i+1;
    end
end

% We assume R/Z/phi equidistant points
nfp = round(2*pi/beam_data.phiaxis(end));
dr = beam_data.raxis(2)-beam_data.raxis(1);
dz = beam_data.zaxis(2)-beam_data.zaxis(1);
dphi = beam_data.phiaxis(2) - beam_data.phiaxis(1);

% Area
area=dr.*dz;

% Volume (function of R)
vol = beam_data.raxis.*dphi.*area;
vol2d=repmat(vol,[1 beam_data.nz]);

dv = zeros(beam_data.nr,beam_data.nphi,beam_data.nz);

for i = 1:beam_data.nphi-1
    grid = squeeze(beam_data.S_ARR(:,i,:));
    grid(grid<=1)=1;
    grid(grid>1)=0;
    dv(:,i,:) = grid.*vol2d;
end

if isempty(ns)
    ns = beam_data.ns_prof1;
end
ds = 1./(ns);
edges = 0:ds:1;
s = 0.5.*(edges(1:end-1)+edges(2:end));

% Create profile
plasma_dvolds = zeros(1,ns-1);
for i=1:ns
    dex1 = beam_data.S_ARR > edges(i);
    dex2 = beam_data.S_ARR <= edges(i+1);
    dex = and(dex1,dex2);
    plasma_dvolds(i) = sum(dv(dex),'all');
end
plasma_dvolds=smooth(plasma_dvolds).*nfp./ds;
plasma_vol = cumsum(plasma_dvolds).*ds;

pp=polyfit(s,plasma_vol,3);
pp2 = polyder(pp);
plasma_dvolds = polyval(pp2,s);
plasma_vol = cumsum(plasma_dvolds).*ds;


end

