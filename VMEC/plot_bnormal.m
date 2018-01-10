function plot_bnormal( coil_data,vmec_data,extcur,varargin)
%PLOT_BNORMAL(coil_data,vmec_data,extcur) Creates a plot of B-Normal.
%   The PLOT_BNORMAL function calculates the normal field on an equilibrium
%   given a vacuum coil set and an equilibria.  The coil_data variable is a
%   coil_data structure as returned by READ_COILS.  The vmec_data variable
%   is a vmec_data struture as read by READ_VMEC.  The extcur variable is
%   an array corresponding to the desired vacuum coil currents for the
%   normal plot.  Additional parameters may then be passed to control how
%   the code exectues.
%   
%   Optional parameters
%      'k',k_val     k_val is the surface on which to work
%      'nu',nu_val   nu_val is the number of poloidal points to use.
%      'nv',nv_val   nv_val is the number of toroidal points to use.
%      'plottype',plottype
%          '2D'      2D pixplot of the B-Normal field.
%          default   3D Isosurface plot.
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.00

% Set Some Defaults
plot_type = 'default';
k      = -1;
ntheta = 90;
nzeta  = 90;


numdefargs=3;   %Number of default arguments
if nargin >numdefargs
    i=1;
    while i<=(nargin-numdefargs)
        if isstruct(varargin{i})
            if isfield(varargin{i},'datatype')
                switch varargin{i}.datatype
                    case {'coil_data'}
                        coil_data=varargin{i};
                end
            end
        else
            switch varargin{i}
                case('k')
                    i=i+1;
                    k=varargin{i};
                case('nu')
                    i=i+1;
                    ntheta=varargin{i};
                case('nv')
                    i=i+1;
                    nzeta=varargin{i};
                case('plottype')
                    i=i+1;
                    plot_type=varargin{i};
            end
        end
        i=i+1;
    end
end

% Prep the coils file
coil_data=coil_biot_prep(coil_data);

% Extract the coordinates
% (do this so in the future we can simply pass a generic equilibrium and
% don't have to rewite the code)
rmnc  = vmec_data.rmnc;
zmns  = vmec_data.zmns;
xm    = vmec_data.xm;
xn    = vmec_data.xn;
rumns = vmec_data.rumns;
zumnc = vmec_data.zumnc;
rvmns = vmec_data.rvmns;
zvmnc = vmec_data.zvmnc;
if (k==-1), k=vmec_data.ns; end;  % If user doesn't pass a surface.


% Transform
theta  = 0:2*pi/(ntheta-1):2*pi;
zeta  = 0:2*pi/(nzeta-1):2*pi;
r    = cfunct(theta,zeta,rmnc,xm,xn);
z    = sfunct(theta,zeta,zmns,xm,xn);
ru   = sfunct(theta,zeta,rumns,xm,xn);
zu   = cfunct(theta,zeta,zumnc,xm,xn);
rv   = sfunct(theta,zeta,rvmns,xm,xn);
zv   = cfunct(theta,zeta,zvmnc,xm,xn);
if vmec_data.iasym
    rmns  = vmec_data.rmns;
    zmnc  = vmec_data.zmnc;
    rumnc = vmec_data.rumnc;
    zumns = vmec_data.zumns;
    rvmnc = vmec_data.rvmnc;
    zvmns = vmec_data.zvmns;
    r    = r  + sfunct(theta,zeta,rmns,xm,xn);
    z    = z  + cfunct(theta,zeta,zmnc,xm,xn);
    ru   = ru + cfunct(theta,zeta,rumnc,xm,xn);
    zu   = zu + sfunct(theta,zeta,zumns,xm,xn);
    rv   = rv + cfunct(theta,zeta,rvmnc,xm,xn);
    zv   = zv + sfunct(theta,zeta,zvmns,xm,xn);
end

% Now extract just the interesting surface and go to xyz coordinates
r = squeeze(r(k,:,:));
z = squeeze(z(k,:,:));
ru = (ru(k,:,:));
zu = (zu(k,:,:));
rv = (rv(k,:,:));
zv = (zv(k,:,:));
for v = 1:nzeta
    x(:,v) = r(:,v).*cos(zeta(v));
    y(:,v) = r(:,v).*sin(zeta(v));
end
x = reshape(x,[1 ntheta*nzeta]);
y = reshape(y,[1 ntheta*nzeta]);
z = reshape(z,[1 ntheta*nzeta]);

% Loop over the equilibria and calculate B on the equilibria
bx = zeros([1 ntheta*nzeta]);
by = zeros([1 ntheta*nzeta]);
bz = zeros([1 ntheta*nzeta]);
for v = 1: ntheta*nzeta
    [bx(v) by(v) bz(v)] = coil_biot(coil_data,x(v),y(v),z(v),extcur);
end

% Now go back to other format
x  = reshape(x,[1 ntheta nzeta]);
y  = reshape(y,[1 ntheta nzeta]);
z  = reshape(z,[1 ntheta nzeta]);
r  = reshape(r,[1 ntheta nzeta]);
z  = reshape(z,[1 ntheta nzeta]);
bx  = reshape(bx,[1 ntheta nzeta]);
by  = reshape(by,[1 ntheta nzeta]);
bz  = reshape(bz,[1 ntheta nzeta]);


% Now make a test plot
%isotoro(r,z,zeta,1,by);

% Calculate the normal vector
snr = r.*zu;
snp = ru.*zv-rv.*zu;
snz = -ru.*r;
for v = 1 : nzeta
   snx(1,:,v) = snr(1,:,v).*cos(zeta(v)) - snp(1,:,v).*sin(zeta(v));
   sny(1,:,v) = snr(1,:,v).*sin(zeta(v)) + snp(1,:,v).*cos(zeta(v));
end
sn = sqrt(snx.*snx+sny.*sny+snz.*snz);
snx = snx./sn;
sny = sny./sn;
snz = snz./sn;

% Calculates B-Normal
bn = bx.*snx+by.*sny+bz.*snz;

% Now make a plot
switch plot_type
    case('fft_mid')
        val  = squeeze(bn(1,1,:));
        L    = nzeta;
        NFFT = 2^nextpow2(L);
        Y    = fft(val,NFFT)/L;
        subplot(2,1,1);
        plot(zeta,val);
        ylabel('B.N midplane');
        xlabel('Toroidal Angle [rad]');
        subplot(2,1,2);
        bar(2*abs(Y(1:NFFT/2+1)));
        ylabel('Amplitude');
        xlabel('Toroidal Mode Number');
    case('test')
    case('2D')
        pixplot(theta,zeta,squeeze(bn));
        xlabel('Toroidal Angle [rad]');
        ylabel('Poloidal Angle [rad]');
        title(['B-Normal (k=' num2str(k,'%3d') ') [T]']);
    otherwise
        isotoro(r,z,zeta,1,bn);
        title(['B-Normal (k=' num2str(k,'%3d') ') [T]']);
end

end

