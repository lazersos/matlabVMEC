function dist = beams3d_getdistrpz(data,r,phi,z,vll,vperp)
%BEAMS3D_GETDISTPRZ Calculates the distribution function
%   This BEAMS3D_GETDISTRPZ function calculates the distribution function
%   at a position in the cylindrical phase space given a BEAMS3D data
%   structure as returned by READ_BEAMS3D.  Either a single point or array
%   of points can be passed to the funciton as R [m], PHI (rad), Z [m],
%   V_parallel [m/s], and vperp [m/s].  Points which fall outside the
%   domain return 0.  The function returns a 2D array of size (NBEAMS,
%   NPTS) where NPTS is the nubmer of requested datapoints.
%
% Example usage
%      % Return an RPZ shaped array
%      beam_data = read_beams3d('beams3d_test.h5');
%      raxis   = 4.5:0.05:6.5;
%      zaxis   = -1.0:0.05:1.0;
%      paxis   = 0:2*pi/40:2*pi;
%      vmax    = beam_data.partvmax;
%      vllaxis = -vmax:2.*vmax./31:vmax;
%      vperpaxis = 0:vmax./15:vmax;
%      [R,P,Z,V,W] = ndgrid(raxis,paxis,zaxis,vllaxis,vperpaxis);
%      nsave = size(R);
%      ntotal = prod(nsave);
%      R = reshape(R,[1 ntotal]);
%      P = reshape(P,[1 ntotal]);
%      Z = reshape(Z,[1 ntotal]);
%      V = reshape(V,[1 ntotal]);
%      W = reshape(W,[1 ntotal]);
%      dist=beams3d_getdistrpz(beam_data,R,P,Z,V,W);
%      dist = reshape(dist,[size(dist,1) nsave]);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       2.00

% Developer note, as of v2.9 BEAMS3D normalizes dist5D by m^3/s^3, meaning
% there is no volume normalization.  Now we can get dV for each flux
% surface easily.  However we have problems when calculating the dV for
% each voxel.  A trivial (but oversimplified) view is to divide by the
% number of voxels (nu*nv).  Which is what we do for now.

% Helpers

% Check version
if data.VERSION<=2.9
    disp('WARNING: Please check the units of the distribution function.');
    disp(data.dist_prof_description);
end

% Create mass array based on beams;
mass = zeros(1,data.nbeams);
for i = 1:data.nbeams
    j = find(data.Beam==i,1,'first');
    mass(i) = data.mass(j);
end

% Interpolate
S   = permute(data.S_ARR,[2 1 3]);
U   = permute(data.U_ARR,[2 1 3]);
pgrid = mod(phi,data.phiaxis(end));
sval= interp3(data.raxis,data.phiaxis,data.zaxis,...
    S,r,pgrid,z);
uval= interp3(data.raxis,data.phiaxis,data.zaxis,...
    U,r,pgrid,z);
rhoval = sqrt(sval);

% We must pad the arrays as data is stored on half-grid mesh
n = [data.ns_prof1+2, data.ns_prof2+2, data.ns_prof3+2, data.ns_prof4+2, data.ns_prof5+2];
n1 = 2:(data.ns_prof1+1);
n2 = 2:(data.ns_prof2+1);
n3 = 2:(data.ns_prof3+1);
n4 = 2:(data.ns_prof4+1);
n5 = 2:(data.ns_prof5+1);
a = data.dist_rhoaxis(1)   - (data.dist_rhoaxis(2)-data.dist_rhoaxis(1));
b = data.dist_rhoaxis(end) + (data.dist_rhoaxis(end)-data.dist_rhoaxis(end-1));
rhoaxis2 = [a data.dist_rhoaxis b];
a = data.dist_uaxis(1)   - (data.dist_uaxis(2)-data.dist_uaxis(1));
b = data.dist_uaxis(end) + (data.dist_uaxis(end)-data.dist_uaxis(end-1));
uaxis2 = [a data.dist_uaxis b];
a = data.dist_paxis(1)   - (data.dist_paxis(2)-data.dist_paxis(1));
b = data.dist_paxis(end) + (data.dist_paxis(end)-data.dist_paxis(end-1));
paxis2 = [a data.dist_paxis b];
a = data.dist_Vaxis(1)   - (data.dist_Vaxis(2)-data.dist_Vaxis(1));
b = data.dist_Vaxis(end) + (data.dist_Vaxis(end)-data.dist_Vaxis(end-1));
Vaxis2 = [a data.dist_Vaxis b];
a = data.dist_Waxis(1)   - (data.dist_Waxis(2)-data.dist_Waxis(1));
b = data.dist_Waxis(end) + (data.dist_Waxis(end)-data.dist_Waxis(end-1));
Waxis2 = [a data.dist_Waxis b];

% Loop over beams
dist=zeros(data.nbeams,length(r));
for i = 1:data.nbeams
    % Pad dist_norm
    dist_norm = zeros(n);
    dist_norm(n1,n2,n3,n4,n5) = squeeze(data.dist_prof(i,:,:,:,:,:));
    % Lower bound
    dist_norm(1,:,:,:,:)    = dist_norm(2,:,:,:,:);
    dist_norm(:,1,:,:,:)    = dist_norm(:,n(2)-1,:,:,:);
    dist_norm(:,:,1,:,:)    = dist_norm(:,:,n(3)-1,:,:);
    %dist_norm(:,:,:,1,:)    = 0;
    dist_norm(:,:,:,:,1)    = dist_norm(:,:,:,:,2);
    % Upper bound
    %dist_norm(n(1),:,:,:,:) = 0;
    dist_norm(:,n(2),:,:,:) = dist_norm(:,2,:,:,:);
    dist_norm(:,:,n(3),:,:) = dist_norm(:,:,2,:,:);
    %dist_norm(:,:,:,n(4),:) = 0;
    %dist_norm(:,:,:,:,n(5)) = 0;
    % Interpolate
    dist(i,:)=interpn(rhoaxis2,uaxis2,paxis2,Vaxis2,Waxis2,dist_norm,rhoval,uval,pval,vll,vperp,'linear',0);
end

end

