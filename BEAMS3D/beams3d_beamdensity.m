function [BEAM,raxis,phiaxis,zaxis] = beams3d_beamdensity(beam_data,varargin)
%BEAMS3D_BEAMDENSITY Calculates the beam density on a grid
%   The BEAMS3D_BEAMDENSITY routine calculates the neutral beam density on
%   a cylindrical grid. Note that the units are [m^-2] as the weights are
%   mutiplied by the particle velocities. The user can specify new axis
%   arrays (R,phi(rad),Z) if the BEAMS3D background grid is not to be used.
%   The code returns the beam density, R grid, phi grid, and Z grid
%   [m,rad,m].
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [BEAM,R,PHI,Z] = beams3d_beamdensity(beam_data);
%
%       raxis = beam_data.raxis;
%       phiaxis = deg2rad(linspace(15,29,64));
%       zaixs = beam_data.zaxis;
%       [BEAM,R,PHI,Z] = beams3d_beamdensity(beam_data,'raxis',raxis,...
%           'phiaxis',phiaxis,'zaxis',zaxis);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.5


%Defaults
nr = double(beam_data.nr);
nphi = double(beam_data.nphi);
nz = double(beam_data.nz);
nbeams = double(beam_data.nbeams);
raxis = beam_data.raxis;
zaxis = beam_data.zaxis;
paxis = beam_data.phiaxis;
pmax  = max(paxis); % Note this is just for rescaling from 0 to B3D grid
weights = beam_data.Weight';


% Handle Varargin
if nargin > 1
    j = 1;
    while j<=length(varargin)
        switch(varargin{j})
            case 'raxis'
                j=j+1;
                raxis = varargin{j};
                nr = length(raxis);
            case 'zaxis'
                j=j+1;
                zaxis = varargin{j};
                nz = length(zaxis);
            case 'phiaxis'
                j=j+1;
                paxis = varargin{j};
                nphi = length(raxis);
        end
        j = j + 1;
    end
end

% Adjust paxis in case it's out of domain
paxis = mod(paxis,pmax);
if(paxis(1)>paxis(end))
    paxis = linspace(min(paxis),max(paxis),nphi);
end

% Setup helpers
rhalf = 0.5.*(raxis(1:end-1)+raxis(2:end));
zhalf = 0.5.*(zaxis(1:end-1)+zaxis(2:end));
phalf = 0.5.*(paxis(1:end-1)+paxis(2:end));

% Calc grid volume
Vol = zeros(nbeams,nr-1,nphi-1,nz-1);
dR  = diff(raxis);
dZ  = diff(zaxis);
dP  = diff(paxis);
dA  = dR*dZ';
dR0 = repmat(rhalf,[1 nz-1]);
for i = 1:nphi-1
    Vol(1,:,i,:) = dA.*dR0.*dP(i);
end
for i = 2:nbeams
    Vol(i,:,:,:) = Vol(1,:,:,:);
end

% Extract beam data
X = beam_data.X_lines(1:2,:);
Y = beam_data.Y_lines(1:2,:);
Z = beam_data.Z_lines(1:2,:);
SOURCE = beam_data.Beam(:)';
VR = beam_data.vr_lines(1,:);
VPHI = beam_data.vphi_lines(1,:);
VZ = beam_data.vz_lines(1,:);
V  = sqrt(VR.*VR+VPHI.*VPHI+VZ.*VZ);
dX = diff(X);
dY = diff(Y);
dZ = diff(Z);
n = sqrt(dX.*dX+dY.*dY+dZ.*dZ);
dX = dX./n;
dY = dY./n;
dZ = dZ./n;
x0 = X(2,:);
y0 = Y(2,:);
z0 = Z(2,:);

% Initialize beam
BEAM = zeros(nbeams,nr-1,nphi-1,nz-1);

% Iterate over beams
for dl = 0:0.01:2
    x1 = x0 - dX.*dl;
    y1 = y0 - dY.*dl;
    z1 = z0 - dZ.*dl;
    r1 = sqrt(x1.*x1+y1.*y1);
    p1 = mod(atan2(y1,x1),pmax);
    n1 = sum(r1 >= raxis,1);
    n2 = sum(p1 >= paxis,1);
    n3 = sum(z1 >= zaxis,1);
    d1 = and(n1>0,n1<nr);
    d2 = and(n2>0,n2<nphi);
    d3 = and(n3>0,n3<nz);
    d  = and(and(d1,d2),d3);
    i1=sub2ind(size(BEAM),SOURCE(d),n1(d),n2(d),n3(d));
    BEAM(i1)=BEAM(i1)+weights(d)./V(d);
end

BEAM = BEAM./Vol;
raxis = rhalf;
phiaxis = phalf;
zaxis = zhalf;


end