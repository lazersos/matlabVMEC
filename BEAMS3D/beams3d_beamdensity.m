function BEAM = beams3d_beamdensity(beam_data)
%BEAMS3D_BEAMDENSITY Calculates the beam density on a grid
%   The BEAMS3D_BEAMDENSITY routine calculates the neutral beam density on
%   a cylindrical grid.  If no grid is provided the background grid is
%   used.
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       BEAM = beams3d_beamdensity(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


%Defaults
nr = double(beam_data.nr);
nphi = double(beam_data.nphi);
nz = double(beam_data.nz);
raxis = beam_data.raxis;
zaxis = beam_data.zaxis;
paxis = beam_data.phiaxis;
pmax  = max(paxis);
weights = beam_data.Weight';

% Initialize beam
BEAM = zeros(nr-1,nphi-1,nz-1);

% Setup helpers
rhalf = 0.5.*(raxis(1:end-1)+raxis(2:end));
zhalf = 0.5.*(zaxis(1:end-1)+zaxis(2:end));
phalf = 0.5.*(paxis(1:end-1)+paxis(2:end));

% Calc grid volume
Vol = zeros(nr-1,nphi-1,nz-1);
dR  = diff(raxis);
dZ  = diff(zaxis);
dP  = diff(paxis);
dA  = dR*dZ';
dR0 = repmat(rhalf,[1 nr-1]);
for i = 1:nphi-1
    Vol(:,i,:) = dA.*dR0.*dP(i);
end

% Extract beam data
X = beam_data.X_lines(1:2,:);
Y = beam_data.Y_lines(1:2,:);
Z = beam_data.Z_lines(1:2,:);
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
    i1=sub2ind(size(BEAM),n1(d),n2(d),n3(d));
    BEAM(i1)=BEAM(i1)+weights(d);
end

BEAM = BEAM./Vol;


end