function beams3d_beamnamelist_temp(vmec_data,energy,power,r_beam,phi_beam,z_beam,div_beam,varargin)
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
beam_dex=[];
note={};
ne_s=[]; ne_f=[];
te_s=[]; te_f=[];
ti_s=[]; ti_f=[];
zeff_s=[]; zeff_f=[];
pot_s=[]; pot_f=[];
nr=128; nz=128; nphi=[];
npoinc=2;
Zatom=1;
species_type='H';
vc_adapt_tol = 0.005;
nparticles_start = 65536;
ec  = 1.60217662E-19; % electron charge [C]
amu = 1.66053906660E-27; % Dalton [kg]
t_end = 1.0E-3;
lplots=0;
fid=[];

if (nargin > 6)
    j=1;
    while j<=numel(varargin)
        if ischar(varargin{j})
            switch(varargin{j})
                case {'H2'}
                    species_type = 'H';
                case {'D2'}
                    species_type = 'D';
                case {'H','D','T','He'}
                    species_type = varargin{j};
                case {'pfrac','power_frac','power_fraction'}
                    j=j+1;
                    pfrac=varargin{j};
                case {'file','filename'}
                    j=j+1;
                    filename=varargin{j};
                    fid=fopen(filename,'a');
                case {'Z','Zatom'}
                    j=j+1;
                    Zatom=varargin{j};
                    charge = charge.*ec;
                case {'mass'}
                    j=j+1;
                    A=varargin{j};
                    mass = A.*amu;
                case {'beam_dex'}
                    j=j+1;
                    beam_dex=varargin{j};
                case {'note'}
                    j=j+1;
                    note=varargin{j};
                case {'plots','plot'}
                    lplots=1;
                case {'TE'}
                    j=j+1;
                    temp=varargin{j};
                    te_s=temp(1,:);
                    te_f=temp(2,:);
                case {'NE'}
                    j=j+1;
                    temp=varargin{j};
                    ne_s=temp(1,:);
                    ne_f=temp(2,:);
                case {'TI'}
                    j=j+1;
                    temp=varargin{j};
                    ti_s=temp(1,:);
                    ti_f=temp(2,:);
                case {'ZEFF'}
                    j=j+1;
                    temp=varargin{j};
                    zeff_s=temp(1,:);
                    zeff_f=temp(2,:);
                case {'POT'}
                    j=j+1;
                    temp=varargin{j};
                    pot_s=temp(1,:);
                    pot_f=temp(2,:);
                case {'t_end'}
                    j=j+1;
                    t_end=varargin{j};
                case {'nr'}
                    j=j+1;
                    nr=varargin{j};
                case {'nphi'}
                    j=j+1;
                    nphi=varargin{j};
                case {'nz'}
                    j=j+1;
                    nz=varargin{j};
                case {'NPOINC'}
                    j=j+1;
                    npoinc=varargin{j};
                case {'VC_ADAPT_TOL'}
                    j=j+1;
                    vc_adapt_tol=varargin{j};
                case {'NPARTICLES_START'}
                    j=j+1;
                    nparticles_start=varargin{j};
            end
        end
        j=j+1;
    end
end

% Hanled fileid
if isempty(fid)
    fid=1;
end
if fid == -1
    fprintf(fid,[' ERROR- Problem opening file: ' filename]);
    return;
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

% Handle bad energy
if (max(energy) > 1)
    energy = energy*ec;
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
if isempty(nphi)
    nphi = max(round(360/(2*vmec_data.nfp)),1);
end

% Handle no profiles passed
if isempty(ne_f)
    ne_s = 0:0.2:1;
    ne_f = polyval([-1.0E20 1.2E20],ne_s);
end
if isempty(te_f)
    te_s = 0:0.2:1;
    te_f = polyval([-1.0E4 1.0E4],te_s);
end
if isempty(ti_f)
    ti_s = 0:0.2:1;
    ti_f = polyval([-1.0E4 1.0E4],ti_s);
end
if isempty(zeff_f)
    zeff_s = 0:0.2:1;
    zeff_f = polyval([1],zeff_s);
end

% Output values to screen
fprintf(fid, '&BEAMS3D_INPUT\n');
fprintf(fid,['  NR = ' num2str(nr,'%d') '\n']);
fprintf(fid,['  NPHI = ' num2str(nphi,'%d') '\n']);
fprintf(fid,['  NZ = ' num2str(nz,'%d') '\n']);
fprintf(fid,['  RMIN = ' num2str(rmin,'%20.10E') '\n']);
fprintf(fid,['  RMAX = ' num2str(rmax,'%20.10E') '\n']);
fprintf(fid,['  ZMIN = ' num2str(zmin,'%20.10E') '\n']);
fprintf(fid,['  ZMAX = ' num2str(zmax,'%20.10E') '\n']);
fprintf(fid, '  PHIMIN = 0.0\n');
fprintf(fid,['  PHIMAX = ' num2str(2.*pi./vmec_data.nfp,'%20.10E') '\n']);
fprintf(fid, '  INT_TYPE = ''LSODE''\n');
fprintf(fid, '  FOLLOW_TOL = 1.0E-9\n');
fprintf(fid,['  VC_ADAPT_TOL = ' num2str(vc_adapt_tol,'%20.10E') '\n']);
fprintf(fid,['  NPOINC = ' num2str(npoinc,'%d') '\n']);
fprintf(fid, '!--------PROFILES ----\n');
fprintf(fid,[ '  NE_AUX_S = ' num2str(ne_s,'%12.6E  ') '\n']);
fprintf(fid,[ '  NE_AUX_F = ' num2str(ne_f,'%12.6E  ') '\n']);
fprintf(fid,[ '  TE_AUX_S = ' num2str(te_s,'%12.6E  ') '\n']);
fprintf(fid,[ '  TE_AUX_F = ' num2str(te_f,'%12.6E  ') '\n']);
fprintf(fid,[ '  TI_AUX_S = ' num2str(ti_s,'%12.6E  ') '\n']);
fprintf(fid,[ '  TI_AUX_F = ' num2str(ti_f,'%12.6E  ') '\n']);
fprintf(fid,[ '  ZEFF_AUX_S = ' num2str(zeff_s,'%12.6E  ') '\n']);
fprintf(fid,[ '  ZEFF_AUX_F = ' num2str(zeff_f,'%12.6E  ') '\n']);
if ~isempty(pot_f)
fprintf(fid,[ '  POT_AUX_S = ' num2str(pot_s,'%12.6E  ') '\n']);
fprintf(fid,[ '  POT_AUX_F = ' num2str(pot_f,'%12.6E  ') '\n']);
end
fprintf(fid,['  PLASMA_MASS = ' num2str(mass,'%20.10E') '\n']);
fprintf(fid,['  PLASMA_ZMEAN = ' num2str(1,'%20.10E') '\n']);
fprintf(fid,['  PLASMA_ZAVG  = ' num2str(1,'%20.10E') '\n']);
fprintf(fid, '!--------Universal Beam Parameters------\n');
fprintf(fid,['  NPARTICLES_START = ' num2str(nparticles_start,'%d') '\n']);
fprintf(fid,['  T_END_IN = ' num2str(ntotal,'%d') '*' num2str(t_end,'%-8.2E') '\n']);
fprintf(fid,['  MASS_BEAMS = ' num2str(ntotal,'%d') '*' num2str(mass,'%-20.10E') '\n']);
fprintf(fid,['  ZATOM_BEAMS = ' num2str(ntotal,'%d') '*' num2str(Zatom,'%-20.10E') '\n']);
fprintf(fid,['  CHARGE_BEAMS = ' num2str(ntotal,'%d') '*' num2str(charge,'%-20.10E') '\n']);
n=1;
for i=1:nbeams
    for j=1:npower
        fprintf(fid,['!----------BEAM ' num2str(n,'%d') ' ----------\n']);
        if ~isempty(note)
            fprintf(fid,['!-- ' note{i} '  Energy=' num2str(energy(i).*(1.0/j)./(ec.*1E3),'%5.2f') '\n']);
        end
        fprintf(fid,['  E_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(energy(i).*(1.0/j),'%20.10E') '\n']);
        fprintf(fid,['  P_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(power(i).*pfrac(j),'%20.10E') '\n']);
        fprintf(fid,['  DIV_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(div_beam(i),'%20.10E') '\n']);
        fprintf(fid,['  ADIST_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(adist(i),'%20.10E') '\n']);
        fprintf(fid,['  ASIZE_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(asize(i),'%20.10E') '\n']);
        fprintf(fid,['  R_BEAMS(' num2str(n,'%2.2d') ',1:2) = ' num2str(r_beam(1:2,i)','%20.10E') '\n']);
        fprintf(fid,['  PHI_BEAMS(' num2str(n,'%2.2d') ',1:2) = ' num2str(phi_beam(1:2,i)',' %20.10E') '\n']);
        fprintf(fid,['  Z_BEAMS(' num2str(n,'%2.2d') ',1:2) = ' num2str(z_beam(1:2,i)',' %20.10E') '\n']);
        if ~isempty(beam_dex)
            fprintf(fid,['  DEX_BEAMS(' num2str(n,'%2.2d') ') = ' num2str(beam_dex(i),'%2i') '\n']);
        end
        n=n+1;
    end
end
fprintf(fid,'/\n');

if lplots
    fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
    thv=0:2*pi./63:2*pi;
    dphi=2*pi./vmec_data.nfp;
    n1=min(min(phi_beam));
    n2=max(max(phi_beam));
    temp=max(min(ceil(0.01+(n2-n1)/dphi)*dphi,2*pi),dphi);
    n1 = (n2-n1)./2; 
    azi=rad2deg(n1);
    zeta = (0.55:1./63:1.0).*temp+n1;%(-1.03:1./63:-0.2).*temp+n1;
    zeta2= (1.75:1./63:3.26).*temp+n1;
    zeta_axis = (0:1./63:5.0).*temp+n1;
    rv = cfunct(thv,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    zv = sfunct(thv,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    rv2 = cfunct(thv,zeta2,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    zv2 = sfunct(thv,zeta2,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    rv_axis = cfunct(thv,zeta_axis,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    zv_axis = sfunct(thv,zeta_axis,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    if vmec_data.iasym==1
        rv = rv+sfunct(thv,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
        zv = zv+cfunct(thv,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
        rv2 = rv+sfunct(thv,zeta2,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
        zv2 = zv+cfunct(thv,zeta2,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
        rv_axis = cfunct(thv,zeta_axis,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
        zv_axis = sfunct(thv,zeta_axis,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
    end
    x_beam = r_beam.*cos(phi_beam);
    y_beam = r_beam.*sin(phi_beam);
%     rv = repmat(rv(:,:,1:end-1),[1 1 vmec_data.nfp]);
%     zv = repmat(zv(:,:,1:end-1),[1 1 vmec_data.nfp]);
%     rv(:,:,end+1) = rv(:,:,1);
%     zv(:,:,end+1) = zv(:,:,1);   
%     phi = (-.5:1./63:4.5).*temp+n1; %0:2*pi/(size(rv,3)-1):2*pi;
    ha=isotoro(rv,zv,zeta,vmec_data.ns);
    hold on
    ha2=isotoro(rv2,zv2,zeta2,vmec_data.ns);
    set(ha,'FaceAlpha',0.33);
    set(ha2,'FaceAlpha',0.33); 
    gcf.GraphicsSmoothing = 'on';
    %plot3(squeeze(rv(1,1,:)).*cos(zeta)',squeeze(rv(1,1,:)).*sin(zeta)',squeeze(zv(1,1,:)),'k');
    plot3(squeeze(rv_axis(1,1,:)).*cos(zeta_axis)',squeeze(rv_axis(1,1,:)).*sin(zeta_axis)',squeeze(zv_axis(1,1,:)),'k');
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
    % Add aperture
    for i=1:nbeams
        appx=(ax(i).*cos(0:2*pi./64:2*pi)+bx(i).*sin(0:2*pi./64:2*pi)).*asize(i)+nx(i).*adist(i)+x_beam(1,i);
        appy=(ay(i).*cos(0:2*pi./64:2*pi)+by(i).*sin(0:2*pi./64:2*pi)).*asize(i)+ny(i).*adist(i)+y_beam(1,i);
        appz=(az(i).*cos(0:2*pi./64:2*pi)+bz(i).*sin(0:2*pi./64:2*pi)).*asize(i)+nz(i).*adist(i)+z_beam(1,i);
        plot3(appx,appy,appz,'k');
    end
    % Add cone (resue app variable)
    for i=1:nbeams
        rv=n(i).*tan(div_beam(i));
        appx=(ax(i).*cos(0:2*pi./64:2*pi)+bx(i).*sin(0:2*pi./64:2*pi)).*rv+nx(i).*n(i)+x_beam(1,i);
        appy=(ay(i).*cos(0:2*pi./64:2*pi)+by(i).*sin(0:2*pi./64:2*pi)).*rv+ny(i).*n(i)+y_beam(1,i);
        appz=(az(i).*cos(0:2*pi./64:2*pi)+bz(i).*sin(0:2*pi./64:2*pi)).*rv+nz(i).*n(i)+z_beam(1,i);
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
fclose(fid);

end

