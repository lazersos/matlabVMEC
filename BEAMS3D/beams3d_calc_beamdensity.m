function beam_density = beams3d_calc_beamdensity(beam_data)
%BEAMS3D_CALC_BEAMDENSITY calculates the spatial neutral beam density
%profiles given a BEAMS3D run with deposition. The beam density is
%calculated on the normal cylindrical background grid. For now, only one
%beam density is calculated, which is the sum over all simulated beams.
%
% Example usage
%      beam_data=read_beams3d('test.h5');
%      beam_density = beams3d_calc_beamdensity(beam_data);
%
% Maintained by: David Kulla (david.kulla@ipp.mpg.de)
% Version:       1.00

nl=128;
norm_nl=ones(nl,1)/nl;
s=linspace(0,1,nl);
rmax=beam_data.raxis(end);


dex=beam_data.end_state~=3&beam_data.end_state~=4;
r0=beam_data.R_lines(1,dex);
phi0=beam_data.PHI_lines(1,dex);
z0=beam_data.Z_lines(1,dex);


r1=beam_data.R_lines(2,dex);
phi1=beam_data.PHI_lines(2,dex);
z1=beam_data.Z_lines(2,dex);

%calculate phi and z for r0=rmax
smax=(rmax-r0)./(r1-r0);
r0=r0+smax.*(r1-r0);
z0=z0+smax.*(z1-z0);
phi0=phi0+smax.*(phi1-phi0);

%Track  the particles in real space
x0 = r0.*cos(phi0);
y0 = r0.*sin(phi0);
x1 = r1.*cos(phi1);
y1 = r1.*sin(phi1);
%m  = beam_data.Beam;

d  = sqrt((x1-x0).^2 + (y1-y0).^2 + (z1-z0).^2);
denbeam = beam_data.Weight(dex)'.*d./beam_data.vll_lines(1,dex); % This is the total number of particles
denp=norm_nl.*denbeam; %Distribute the weight over all points along the track
xl = x0 + s'.*(x1-x0);
yl = y0 + s'.*(y1-y0);
zl = z0 + s'.*(z1-z0);
rl = sqrt(xl.*xl+yl.*yl);
pl = atan2(yl,xl);
pl = mod(pl,beam_data.phiaxis(end));


%Use axis as edges and adjust for bin width
hp2=(beam_data.phiaxis(2)-beam_data.phiaxis(1))/2;
hr2=(beam_data.raxis(2)-beam_data.raxis(1))/2;
hz2=(beam_data.zaxis(2)-beam_data.zaxis(1))/2;
[discphi,~]=discretize(pl(:),beam_data.phiaxis-hp2);
[discr,~]=discretize(rl(:),beam_data.raxis-hr2);
[discz,~]=discretize(zl(:),beam_data.zaxis-hz2);
discr(isnan(discr))=1;
discphi(isnan(discphi))=1;
discz(isnan(discz))=1;

beam_density = accumarray([discr,discphi,discz],denp(:),size(beam_data.B_PHI));

end
