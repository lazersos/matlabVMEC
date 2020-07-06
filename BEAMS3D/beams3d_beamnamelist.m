function beams3d_beamnamelist(vmec_data,energy,power,r_beam,phi_beam,z_beam,div_beam,varargin)
%BEAMS3D_BEAMNAMELIST Creates an beam input namelist from VMEC run
%   The BEASM3DINPUTNAMELIST function outputs to screen an BEASM3D input
%   namelist based on beam information and a VMEC run.  It takes as
%   input a vmec_data data structure as returned by READ_VMEC.  Energy is
%   specified in [eV], power in [W], r_beam in [m], phi_beam in [rad],
%   z_beam in [m], div_beam in [rad]. Optional arguments include species
%   type ('H' (default), 'D', 'T', or 'He'), power fracions, and whether
%   plots are requested 'plots'.
%
%   Example usage
%       vmec_data=read_vmec('wout_test.nc');
%       r_min = vmec_data.rmin_surf;
%       r_max = vmec_data.rmax_surf+vmec_data.Aminor;
%       energy = [55e3 55e3];
%       power  = [1.7e6 1.8e6];
%       div_beam = [0.05 0.05];
%       r_beam = [r_max r_max; r_min r_min]; 
%       p_beam = [0.0 2*pi; 0.0 pi]./vmec_data.nfp;
%       z_beam = [0.0 0.0; 0.0 0.0];
%       beams3d_beamnamelist(vmec_data,energy,power,r_beam,p_beam,z_beam,div_beam,'D','plots')
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


% Defaults
mass=[];
charge=[];
ashape=[];
pfrac=[];
Zatom=1;
species_type='H';
ec  = 1.60217662E-19; % electron charge [C]
amu = 1.66053906660E-27; % Dalton [kg]
lplots=0;

if (nargin > 6)
    j=1;
    while j<=numel(varargin)
        switch(varargin{j})
            case {'H','D','T','He'}
                species_type = varargin{j};
            case {'pfrac','power_frac','power_fraction'}
                j=j+1;
                pfrac=varargin{j};
            case {'Z','Zatom'}
                j=j+1;
                Zatom=varargin{j};
                charge = charge.*ec;
            case {'mass'}
                j=j+1;
                A=varargin{j};
                mass = A.*amu;
            case {'plots','plot'}
                lplots=1;
        end
        j=j+1;
    end
end

% Handle mass
if isempty(mass)
    switch species_type
        case 'H'
            A   = 1.00784; % Hydrogen mass [amu]
        case 'D'
            A   = 2.01410177811; % Deuterium mass [amu]
        case 'T'
            A   = 3.0160492; % Tritium mass [amu]
        case 'He'
            A   = 4.002602; % Helium mass [amu]
    end
    mass = A.*amu;
end
% Handle charge
if isempty(charge)
    switch species_type
        case {'H','D','T'}
            Zatom   = 1; % Hydrogen mass [amu]
        case 'He'
            Zatom   = 2; % Helium mass [amu]
    end
    charge = Zatom.*ec;
end
% Handle Power fractions
if isempty(pfrac)
    switch species_type
        case {'H','D','T'}
            pfrac=[0.5 0.4 0.1];
        case 'He'
            pfrac=1.0;
    end
end

%Setup beams
nbeams = length(energy);
npower = length(pfrac);
ntotal = nbeams*npower;

% Handle asize
if isempty(ashape)
    asize=ones(1,nbeams).*0.25;
    adist=zeros(1,nbeams);
end


% VMEC helpers
delta = 0.2.*vmec_data.Aminor;
rmin = vmec_data.rmin_surf-delta;
rmax = vmec_data.rmax_surf+delta;
zmax = vmec_data.zmax_surf+delta;
if (vmec_data.iasym==1)
    zmin = vmec_data.zmin_surf-delta;
else
    zmin = -vmec_data.zmax_surf-delta;
end

% Output values to screen
disp( '&BEAMS3D_INPUT');
disp( '  NR = 128');
disp(['  NPHI = ' num2str(360/(2*vmec_data.nfp),'%d')]);
disp( '  NZ = 128');
disp(['  RMIN = ' num2str(rmin,'%20.10E')]);
disp(['  RMAX = ' num2str(rmax,'%20.10E')]);
disp(['  ZMIN = ' num2str(zmin,'%20.10E')]);
disp(['  ZMAX = ' num2str(zmax,'%20.10E')]);
disp( '  PHIMIN = 0.0');
disp(['  PHIMAX = ' num2str(2.*pi./vmec_data.nfp,'%20.10E')]);
disp( '  INT_TYPE = ''LSODE''');
disp( '  FOLLOW_TOL = 1.0E-9');
disp( '  VC_ADAPT_TOL = 1.0E-2');
disp( '  NPOINC = 2');
disp( '!--------Dummy profiles Te=Ti=10kev ----');
disp( '  NE_AUX_S = 0.00E+00  0.20E+00  0.40E+00  0.60E+00  0.80E+00  1.00E+00');
disp( '  NE_AUX_F = 1.20E+20  1.00E+00  0.80E+20  0.60E+20  0.40E+20  0.20E+20');
disp( '  TE_AUX_S = 0.00E+00  0.20E+00  0.40E+00  0.60E+00  0.80E+00  1.00E+00');
disp( '  TE_AUX_F = 1.00E+04  0.80E+04  0.60E+04  0.40E+04  0.20E+04  0.00E+00');
disp( '  TI_AUX_S = 0.00E+00  0.20E+00  0.40E+00  0.60E+00  0.80E+00  1.00E+00');
disp( '  TI_AUX_F = 1.00E+04  0.80E+04  0.60E+04  0.40E+04  0.20E+04  0.00E+00');
disp( '  ZEFF_AUX_S = 0.00E+00  0.20E+00  0.40E+00  0.60E+00  0.80E+00  1.00E+00');
disp( '  ZEFF_AUX_F = 1.00E+00  1.00E+00  1.00E+00  1.00E+00  1.00E+00  1.00E+00');
disp(['  PLASMA_MASS = ' num2str(mass,'%20.10E')]);
disp(['  PLAMSA_ZMEAN = ' num2str(1,'%20.10E')]);
disp(['  PLAMSA_ZAVG  = ' num2str(1,'%20.10E')]);
disp( '!--------Universal Beam Parameters------');
disp( '  NPARTICLES_START = 65536');
disp(['  T_END_IN = ' num2str(ntotal,'%d') '*1.0E-3']);
disp(['  MASS_BEAMS = ' num2str(ntotal,'%d') '*' num2str(mass,'%-20.10E')]);
disp(['  ZATOM_BEAMS = ' num2str(ntotal,'%d') '*' num2str(Zatom,'%-20.10E')]);
disp(['  CHARGE_BEAMS = ' num2str(ntotal,'%d') '*' num2str(charge,'%-20.10E')]);
n=1;
for i=1:nbeams
    for j=1:npower
        disp(['!----------BEAM ' num2str(n,'%d') ' ----------']);
        disp(['  E_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(energy(i).*(1.0/j),'%20.10E')]);
        disp(['  P_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(power(i).*pfrac(j),'%20.10E')]);
        disp(['  DIV_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(div_beam(i),'%20.10E')]);
        disp(['  ADIST_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(adist(i),'%20.10E')]);
        disp(['  ASIZE_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(asize(i),'%20.10E')]);
        disp(['  R_BEAMS(' num2str(n,'%2.2d') ',1:2) = ' num2str(r_beam(1:2,i)','%20.10E')]);
        disp(['  PHI_BEAMS(' num2str(n,'%2.2d') ',1:2) = ' num2str(phi_beam(1:2,i)',' %20.10E')]);
        disp(['  Z_BEAMS(' num2str(n,'%2.2d') ',1:2) = ' num2str(z_beam(1:2,i)',' %20.10E')]);
        n=n+1;
    end
end
disp('/');

if lplots
    fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
    thv=0:2*pi./63:2*pi;
    dphi=2*pi./vmec_data.nfp;
    n1=min(min(phi_beam));
    n2=max(max(phi_beam));
    temp=max(min(ceil(0.01+(n2-n1)/dphi)*dphi,2*pi),dphi);
    n1 = (n2-n1)./2; azi=rad2deg(n1);
    zeta = (-.5:1./63:.5).*temp+n1;
    rv = cfunct(thv,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    zv = sfunct(thv,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    if vmec_data.iasym==1
        rv = rv+sfunct(thv,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
        zv = zv+cfunct(thv,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
    end
    x_beam = r_beam.*cos(phi_beam);
    y_beam = r_beam.*sin(phi_beam);
    ha=isotoro(rv,zv,zeta,vmec_data.ns);
    set(ha,'FaceAlpha',0.33); hold on;
    plot3(squeeze(rv(1,1,:)).*cos(zeta)',squeeze(rv(1,1,:)).*sin(zeta)',squeeze(zv(1,1,:)),'k');
    plot3(x_beam,y_beam,z_beam);
    nx=diff(x_beam);
    ny=diff(y_beam);
    nz=diff(z_beam);
    n=sqrt(nx.^2+ny.^2+nz.^2);
    nx=nx./n; ny=ny./n; nz=nz./n;
    ax=ny;
    ay=-nx;
    az=ax.*0;
    bx=ay.*nz-az.*ny;
    by=az.*nx-ax.*nz;
    bz=ax.*ny-ay.*nx;
    %quiver3(x_beam(1,:),y_beam(1,:),z_beam(1,:),nx,ny,nz)
    %quiver3(x_beam(1,:),y_beam(1,:),z_beam(1,:),ax,ay,az)
    %quiver3(x_beam(1,:),y_beam(1,:),z_beam(1,:),bx,by,bz)
    % Add apperature
    for i=1:nbeams
        appx=(ax(i).*cos(0:2*pi./64:2*pi)+bx(i).*sin(0:2*pi./64:2*pi)).*asize(i)+nx(i).*adist(i)+x_beam(1,i);
        appy=(ay(i).*cos(0:2*pi./64:2*pi)+by(i).*sin(0:2*pi./64:2*pi)).*asize(i)+ny(i).*adist(i)+y_beam(1,i);
        appz=(az(i).*cos(0:2*pi./64:2*pi)+bz(i).*sin(0:2*pi./64:2*pi)).*asize(i)+nz(i).*adist(i)+z_beam(1,i);
        plot3(appx,appy,appz,'k');
    end
    % Add cone (resue app variable)
    for i=1:nbeams
        r=n(i).*tan(div_beam(i));
        appx=(ax(i).*cos(0:2*pi./64:2*pi)+bx(i).*sin(0:2*pi./64:2*pi)).*r+nx(i).*n(i)+x_beam(1,i);
        appy=(ay(i).*cos(0:2*pi./64:2*pi)+by(i).*sin(0:2*pi./64:2*pi)).*r+ny(i).*n(i)+y_beam(1,i);
        appz=(az(i).*cos(0:2*pi./64:2*pi)+bz(i).*sin(0:2*pi./64:2*pi)).*r+nz(i).*n(i)+z_beam(1,i);
        x1x2x3=[]; y1y2y3=[]; z1z2z3=[];
        for j=1:length(appx)-1
            x1x2x3(:,j)=[x_beam(1,i) appx(j) appx(j+1)];
            y1y2y3(:,j)=[y_beam(1,i) appy(j) appy(j+1)];
            z1z2z3(:,j)=[z_beam(1,i) appz(j) appz(j+1)];
        end
        ha=patch('XData',x1x2x3,'YData',y1y2y3,'ZData',z1z2z3);
        set(ha,'FaceColor','blue','FaceAlpha',0.33,'LineStyle','none');
    end
    title('');
    view(azi+90,30);
    set(gca,'Clipping','off');
    axis off;
end

end

