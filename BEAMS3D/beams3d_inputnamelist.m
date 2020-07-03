function beams3d_inputnamelist(vmec_data,energy,pitch,rho,theta,zeta,varargin)
%BEAMS3D_INPUTNAMELIST Creates an input namelist from VMEC run
%   The BEASM3DINPUTNAMELIST function outputs to screen an BEASM3D input
%   namelist based on particle information and a VMEC run.  It takes as
%   input a vmec_data data structure as returned by READ_VMEC.  Energy is
%   specified in [eV], pitch in degrees, rho in r/a, theta in radian, and
%   zeta in radians.  Optional arguments include species type ('H'
%   (default), 'D', 'T', or 'He') and whether plots are requested 'plots'.
%
%   Example usage
%       vmec_data=read_vmec('wout_test.nc');
%       energy = 55e3;
%       pitch  = -85:5:85;
%       rho    = [0.25 0.5 0.75];
%       theta  = [0 pi];
%       zeta   = [0 pi];
%       beams3d_inputnamelist(vmec_data,energy,pitch,rho,theta,zeta,'D','plots');
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


% Defaults
mass=[];
charge=[];
Zatom=1;
species_type='H';
ec  = 1.60217662E-19; % electron charge [C]
amu = 1.66053906660E-27; % Dalton [kg]
lplots=0;

if length(energy) >1
    disp('At this time only one energy may be specified');
    return;
end

if (nargin > 6)
    j=1;
    while j<=numel(varargin)
        switch(varargin{j})
            case {'H','D','T','He'}
                species_type = varargin{j};
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
if isempty(charge)
    switch species_type
        case {'H','D','T'}
            Zatom   = 1; % Hydrogen mass [amu]
        case 'He'
            Zatom   = 2; % Helium mass [amu]
    end
    charge = Zatom.*ec;
end

% Convert energy
V = sqrt(2.*energy.*ec./mass);

%Setup spatial arrays
npitch = length(pitch);
nrho   = length(rho);
ntheta = length(theta);
nzeta  = length(zeta);
[rho3, theta3, zeta3]=ndgrid(rho,theta,zeta);
ntotal=numel(rho3);
rho_lin = reshape(rho3,[1 ntotal]);
theta_lin = reshape(theta3,[1 ntotal]);
zeta_lin = reshape(zeta3,[1 ntotal])/vmec_data.nfp;

% VMEC helpers
s=0:1/(vmec_data.ns-1):1;
delta = 0.2.*vmec_data.Aminor;
rmin = vmec_data.rmin_surf-delta;
rmax = vmec_data.rmax_surf+delta;
zmax = vmec_data.zmax_surf+delta;
if (vmec_data.iasym==1)
    zmin = vmec_data.zmin_surf-delta;
else
    zmin = -vmec_data.zmax_surf-delta;
end

% Calculate particle values
r=[]; z=[]; phi=[]; vll=[]; vperp=[]; mu=[];
for i=1:numel(rho3)
    rhot=rho_lin(i); u=theta_lin(i); v=zeta_lin(i);
    rtemp=cfunct(u,v,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    ztemp=sfunct(u,v,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    btemp=cfunct(u,v,vmec_data.bmnc,vmec_data.xm_nyq,vmec_data.xn_nyq);
    r=[r ones(1,npitch).*pchip(s,rtemp,rhot.*rhot)];
    z=[z ones(1,npitch).*pchip(s,ztemp,rhot.*rhot)];
    phi=[phi ones(1,npitch).*v];
    vll=[vll V.*sind(pitch)];
    vperp=[vperp V.*cosd(pitch)];
    mu =[mu 0.5.*mass.*(V.*cosd(pitch)).^2./(pchip(s,btemp,rhot.*rhot).*ones(1,npitch))];
end

% Output values to screen
disp(['&BEAMS3D_INPUT']);
disp(['  NR = 128']);
disp(['  NPHI = ' num2str(360/(2*vmec_data.nfp),'%d')]);
disp(['  NZ = 128']);
disp(['  RMIN = ' num2str(rmin,'%20.10E')]);
disp(['  RMAX = ' num2str(rmax,'%20.10E')]);
disp(['  ZMIN = ' num2str(zmin,'%20.10E')]);
disp(['  ZMAX = ' num2str(zmax,'%20.10E')]);
disp(['  PHIMIN = 0.0']);
disp(['  PHIMAX = ' num2str(pi./vmec_data.nfp,'%20.10E')]);
disp(['  INT_TYPE = ''LSODE''']);
disp(['  FOLLOW_TOL = 1.0E-9']);
disp(['  VC_ADAPT_TOL = 1.0E-2']);
disp(['  NPOINC = 1000']);
disp(['  R_START_IN = ' num2str(r,' %20.10E')]);
disp(['  Z_START_IN = ' num2str(z,' %20.10E')]);
disp(['  PHI_START_IN = ' num2str(phi,' %20.10E')]);
disp(['  VLL_START_IN = ' num2str(vll,' %20.10E')]);
disp(['  MU_START_IN = ' num2str(mu,' %20.10E')]);
disp(['  MASS_IN = ' num2str(length(r),'%d') '*' num2str(mass,'%-20.10E')]);
disp(['  CHARGE_IN = ' num2str(length(r),'%d') '*' num2str(charge,'%-20.10E')]);
disp(['  ZATOM_IN = ' num2str(length(r),'%d') '*' num2str(Zatom,'%-3.1d')]);
disp(['  T_END_IN = ' num2str(length(r),'%d') '*1.0E-3']);
disp(['/']);

if lplots
    thv=0:2*pi./64:2*pi;
    rv = cfunct(thv,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    zv = sfunct(thv,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
    if vmec_data.iasym==1
        rv = rv+sfunct(thv,zeta,vmec_data.rmns,vmec_data.xm,vmec_data.xn);
        zv = zv+cfunct(thv,zeta,vmec_data.zmnc,vmec_data.xm,vmec_data.xn);
    end
    n1=1;
    nstep = length(r)/nzeta;
    n2 = nstep;
    for i=1:nzeta
        fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
        plot(rv(vmec_data.ns,:,i),zv(vmec_data.ns,:,i),'r','LineWidth',4); hold on;
        plot(rv(vmec_data.ns/4,:,i),zv(vmec_data.ns/4,:,i),'--k','LineWidth',2);
        plot(rv(1,1,i),zv(1,1,i),'+k','LineWidth',2,'MarkerSize',10);
        plot(r(n1:n2),z(n1:n2),'ok','Linewidth',2,'MarkerSize',10);
        xlim([rmin rmax]);
        ylim([zmin zmax]);
        set(gca,'FontSize',24);
        xlabel('R [m]');
        ylabel('Z [m]');
        title(['BEASM3D Starting Points (\zeta=' num2str(rad2deg(zeta(i)),'%d') '°)']);
        n1=n2+1;
        n2=n2+nstep;
    end
    fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
    plot(vll./1E6,vperp./1E6,'ok','LineWidth',2,'MarkerSize',10);
    vmax=max(max(abs(vll)),max(vperp))./1E6;
    xlim([-1.1 1.1].*vmax);
    ylim([0 1.1].*vmax);
    set(gca,'FontSize',24);
    xlabel('V_{parallel} x10^{6} [m/s]');
    ylabel('V_{perpendicular} x10^{6} [m/s]');
    title('BEAMS3D Starting Points Pitch');
end

end

