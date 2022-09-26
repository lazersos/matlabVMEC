function dist = beams3d_getdistrpzEpitch(data,r,phi,z,E,pitch)
%BEAMS3D_GETDISTPRZ Calculates the distribution function in energy/pitch
%   This BEAMS3D_GETDISTRPZ function calculates the distribution function
%   at a position in the cylindrical phase space given a BEAMS3D data
%   structure as returned by READ_BEAMS3D.  Either a single point or array
%   of points can be passed to the funciton as R [m], PHI (rad), Z [m],
%   Energy [eV], and pitch (vll/v).  Points which fall outside the
%   domain return 0.  The function returns a 2D array of size (NBEAMS,
%   NPTS) where NPTS is the nubmer of requested datapoints. The
%   distribution function units are the same as that of FIDASIM:
%   [1/cm^3/keV/(dOmega/4pi(?))]
%
% Example usage
%      % Return an RPZ shaped array
%      beam_data = read_beams3d('beams3d_test.h5');
%      raxis   = 4.5:0.05:6.5;
%      zaxis   = -1.0:0.05:1.0;
%      paxis   = 0:2*pi/40:2*pi;
%      Eaxis   = 0:10E3:100E3;
%      pitchaxis = -1:0.1:1;
%      [R,P,Z,E,PITCH] = ndgrid(raxis,paxis,zaxis,Eaxis,pitchaxis);
%      nsave = size(R);
%      ntotal = prod(nsave);
%      R = reshape(R,[1 ntotal]);
%      P = reshape(P,[1 ntotal]);
%      Z = reshape(Z,[1 ntotal]);
%      E = reshape(E,[1 ntotal]);
%      PITCH = reshape(PITCH,[1 ntotal]);
%      dist=beams3d_getdistrpzEpitch(beam_data,R,P,Z,E,PITCH);a
%      dist = reshape(dist,[size(dist,1) nsave]);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.20

% Developer note, as of v2.9 BEAMS3D normalizes dist5D by m^3/s^3, meaning
% there is no volume normalization.  Now we can get dV for each flux
% surface easily.  However we have problems when calculating the dV for
% each voxel.  A trivial (but oversimplified) view is to divide by the
% number of voxels (nu*nv).  Which is what we do for now.

% Helpers
ec = 1.60217662E-19;
cspeed = 299792458;
ds=1.0;
du=2*pi;
dp=2*pi;
vmax = data.partvmax;
dv=2*data.partvmax;
dw=data.partvmax;

mass = data.mass(1); % Assume same mass
V = sqrt(ec.*2.*E./mass);
Vval = V.*pitch; %Vll
Wval = sqrt(V.^2-Vval.^2); %Vperp
%Vval = Vval+vmax;



% Interpolate
S   = permute(data.S_ARR,[2 1 3]);
U   = permute(data.U_ARR,[2 1 3]);
pgrid = mod(phi,data.phiaxis(end));
sval= interp3(data.raxis,data.phiaxis,data.zaxis,...
    S,r,pgrid,z);
uval= interp3(data.raxis,data.phiaxis,data.zaxis,...
    U,r,pgrid,z);
pval = mod(phi,2*pi);
rhoval = sqrt(sval);

% Jacobian from D. Moseev paper
% https://doi.org/10.1063/1.5085429
%jac = 1.0./(mass.*sqrt(1-pitch.*pitch)); % this is wrong for the current
%beams3d normalization (f2D_Ep = 2pi v/m f3D_Car)

jac = 2* pi * V ./mass .* ec / 1000;
dist_norm = squeeze(sum(data.dist_prof,1));

uval(uval<min(data.dist_uaxis))=min(data.dist_uaxis);
uval(uval>max(data.dist_uaxis))=max(data.dist_uaxis);

pval(pval<min(data.dist_paxis))=min(data.dist_paxis);
pval(pval>max(data.dist_paxis))=max(data.dist_paxis);

dist=interpn(data.dist_rhoaxis,data.dist_uaxis,data.dist_paxis,data.dist_Vaxis,data.dist_Waxis,dist_norm,rhoval,uval,pval,Vval,Wval,'linear');

dist = dist.*jac;

end