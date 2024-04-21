function dist = beams3d_getdistEpitch(data,rho,u,phi,E,pitch)
%BEAMS3D_GETDISTEPITCH Calculates the distribution function in energy/pitch
%   This BEAMS3D_GETDISTEPITCH function calculates the distribution 
%   function at a position in the equilibrium space given a BEAMS3D data
%   structure as returned by READ_BEAMS3D.  Either a single point or array
%   of points can be passed to the funciton as rho [arb], u (rad), 
%   phi [rad], Energy [eV], and pitch (vll/v).  Points which fall outside
%   the domain return 0.  The function returns a 2D array of size (NBEAMS,
%   NPTS) where NPTS is the nubmer of requested datapoints.
%
% Example usage
%      % Return an RPZ shaped array
%      beam_data = read_beams3d('beams3d_test.h5');
%      rhoaxis   = linspace(0,1,32);
%      uaxis   = linspace(0,2.*pi,16);
%      paxis   = linspace(0,2.*pi,40);
%      Eaxis   = 0:10E3:100E3;
%      pitchaxis = -1:0.1:1;
%      [RHO,U,PHI,E,PITCH] = ndgrid(rhoaxis,uaxis,paxis,Eaxis,pitchaxis);
%      nsave = size(rhoaxis);
%      ntotal = prod(nsave);
%      RHO = reshape(RHO,[1 ntotal]);
%      U = reshape(U,[1 ntotal]);
%      PHI = reshape(PHI,[1 ntotal]);
%      E = reshape(E,[1 ntotal]);
%      PITCH = reshape(PITCH,[1 ntotal]);
%      dist=beams3d_getdistEpitch(beam_data,RHO,U,PHI,E,PITCH);
%      dist = reshape(dist,[size(dist,1) nsave]);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       2.0

% Helpers
ec = 1.60217662E-19;

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

% Loop over beams calculating distribution
dist=zeros(data.nbeams,length(r));
for i = 1:data.nbeams
    % Calc Jacobian
    V = sqrt(ec.*2.*E./mass(i));
    Vval = V.*pitch; %Vll
    Wval = sqrt(V.^2-Vval.^2); %Vperp
    % Jacobian from D. Moseev paper
    % https://doi.org/10.1063/1.5085429
    jac = 2 * pi * V .*ec ./mass(i); % eq 29, ec for J->eV
    % Setup distribution
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
    dist(i,:)=interpn(rhoaxis2,uaxis2,paxis2,Vaxis2,Waxis2,dist_norm,rho,u,phi,Vval,Wval,'linear',0).*jac;
end

end