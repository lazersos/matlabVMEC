function [intensity, Ti, B, S] = beams3d_lineint(beam_data,r0,p0,z0,r1,p1,z1,width,varargin)
%BEASM3D_LINEINT Performs line integral across birth locations
%   The BEASM3D_LINEINT function returns beam deposition weighted line
%   integrals along a LOS.  It takes a beams3d data structure as returned
%   by READ_BEAMS3D, cylindical (r0,p0,z0) starting and ending (r1,p1,z1)
%   positions, and a LOS width parameter as input.  It returns the beam
%   deposition weighted particle number, ion temperature, magnetic field,
%   and radial position.
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [intenstiy, Ti, B, S] = beams3d_lineint(beam_data,r0,p0,z0,r1,p1,z1,width);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0

lplot=0;
res = 0.01;
intensity=[];
Ti=[];
B=[];
beamdex=1:double(max(beam_data.Beam));

% Handle varargin
if ~isempty(varargin)
    i=1;
    while i<=length(varargin)
        if isstr(varargin{i})
            switch varargin{i}
                case 'plots'
                    lplot=1;
                case 'beamdex'
                    i=i+1;
                    beamdex=varargin{i};
            end
        end
        i=i+1;
    end
end

% Handle geometry
x0=r0.*cos(p0);
y0=r0.*sin(p0);
x1=r1.*cos(p1);
y1=r1.*sin(p1);
nx=x1-x0; ny=y1-y0; nz=z1-z0;
n=sqrt(nx.*nx+ny.*ny+nz.*nz);
nstep=n./res;
nx=nx./n; ny=ny./n; nz=nz./n;

% Select beam data (and downselect)
dex = beams3d_finddex(beam_data,'orbit_birth');
r_beam = beam_data.R_lines(2,dex~=0);
phi_beam = beam_data.PHI_lines(2,dex~=0);
z_beam = beam_data.Z_lines(2,dex~=0);
b_beam = beam_data.B_lines(2,dex~=0);
s_beam = beam_data.S_lines(2,dex~=0);
w_beam = beam_data.Weight(dex~=0)';
x_beam = r_beam.*cos(phi_beam);
y_beam = r_beam.*sin(phi_beam);

% Create helpers
temp   = permute(beam_data.TI,[2 1 3]);
ti_beam= interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
    temp,r_beam,mod(phi_beam,beam_data.phiaxis(end)),z_beam);

if (lplot)
    fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
    subplot(3,4,[1 2 5 6 9 10]);
    plot3([x0,x1]',[y0,y1]',[z0,z1]');
    hold on;
    plot3(x_beam,y_beam,z_beam,'.k');
    axis equal;
end

% Loop over chords
nchords = length(x0);
intensity = zeros(1,nchords);
Ti = zeros(1,nchords);
B = zeros(1,nchords);
S = 1.5.*ones(1,nchords);
for i=1:nchords
    xt=x0(i); yt=y0(i); zt=z0(i);
    mask = zeros(1,length(x_beam));
    for j = 1:nstep(i)
        xt = xt+res.*nx(i);
        yt = yt+res.*ny(i);
        zt = zt+res.*nz(i);
        dx = xt-x_beam;
        dy = yt-y_beam;
        dz = zt-z_beam;
        d  = sqrt(dx.*dx+dy.*dy+dz.*dz);
        mask(d<=width) = 1;
    end
    n = sum(mask);
    if n==0, continue; end
    dex = mask==1;
    intensity(i) = sum(w_beam(dex));
    Ti(i) = sum(w_beam(dex).*ti_beam(dex))./intensity(i);
    B(i) = sum(w_beam(dex).*b_beam(dex))./intensity(i);
    S(i) = sum(w_beam(dex).*s_beam(dex))./intensity(i);
end

if (lplot)
    subplot(3,4,[3 4]);
    plot(sqrt(S),intensity,'ok','LineWidth',2,'MarkerSize',12);
    set(gca,'FontSize',24);
    xlabel('r/a');
    ylabel('Intensity');
    subplot(3,4,[7 8]);
    plot(sqrt(S),Ti,'ok','LineWidth',2,'MarkerSize',12);
    set(gca,'FontSize',24);
    xlabel('r/a');
    ylabel('T_i');
    subplot(3,4,[11 12]);
    plot(sqrt(S),B,'ok','LineWidth',2,'MarkerSize',12);
    set(gca,'FontSize',24);
    xlabel('r/a');
    ylabel('|B|');
    
    
end
return;

end

