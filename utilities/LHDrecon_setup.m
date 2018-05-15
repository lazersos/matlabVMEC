function LHDrecon_setup(data,timedex1,timedex2,varargin)
%LHDrecon_setup(data,timedex1,timedex2) Create VMEC/STELLOPT namelist file
%
%   LHDrecon_setup(data,timedex1,timedex2)  Creates a VMEC input namelist
%   file from the values found in data.  Mean values from timedex1 to
%   timedex2 are utilized.
%
%   LHDrecon_setup(data,timedex1,timedex2,'recon')  Creates a STELLOPT
%   input namelist file from the values found in data.  Mean values from
%   timedex1 to timedex2 are utilized.
%       Options:
%           'recon'         STELLOPT reconstruction.
%           'shotnum'       Will use the next input variable as a name for
%                           the file.
%           'user'          The next argument should be the initials of the
%                           user creating the shot.
%           'magdiags'      Include matching of magnetic diagnostics.
%           'eplasma'       Target stored energy.
%           'diamag'        Target diamagnetic loop signal.
%           'curtor'        Target total toroidal current.
%           'isote'         Target iso-Te constraint. (under development)
%           'iota'          Optimize rotational transform profile instead
%                           of current profile.
%           'msedata'       Include matching of MSE data.  User must also
%                           pass a coil_data structure for vacuum field
%                           calculation.  User may also pass stellopt
%                           calculated MSE data by passing a structure with
%                           datatype 'stellopt_mse', as produced by
%                           read_stellopt routine.
%           'msevac'        If passed, a single iteration vacuum run will
%                           be created in order to facilitate vacuum MSE
%                           calculations.
%           ves_data        If the user passes a vessel data structure it
%                           be used as an optimization target.
%           'pres_splines'  Use Akima Splines to define VMEC pressure.
%           'curr_splines'  Use Akima Splines to define VMEC I'.
%           'iota_splines'  Use Akima Splines to define VMEC rotational
%                           transform. (automatically sets iota
%                           optimization)
%           'plots'         Will produce plots of various quantities.
%
%   See also read_stellopt, read_vessel, read_coil, and read_netcdf.
%
%   Example:
%      shot_data=read_netcdf('LHD85384.nc');
%      coil_data=read_coil('coil.lhd_ys4');
%      ves_data=read_vessel('vessel.dat');
%      LHDrecon_setup(shot_data,1.95,2.05,'recon','shotnum','LHD85384','user','SAL','pres_splines','curr_splines',ves_data,coil_data,'msedata','magdiags','plots')
%
%   Written by: Samuel Lazerson (lazerson@pppl.gov)
%   Version:    1.5
%   Date:       1/5/2012
 
warning off
% Defaults
shotnum='LHDshot';
filename=['input.' shotnum '_'];
outputtype='VMEC';
percf='%+7.4f';
perct='%04.4i';
percd='%6.6d';
user='XXX';
bigno=1e30;
vesdata=[];
vmec_data=[];
coil_data=[];
mse_vac=[];
makeplots=0;
pcurr_type = 'power_series';
piota_type = 'power_series';
pmass_type = 'power_series';
signfactor=-1;    % FOR LHD SHOTS (Applied to curtor and flux loop signals)
sigma_diagno_flx_default=[  1.0 1.0 1.0 bigno 1.0 bigno...
                    1.0 1.0 1.0 1.0 1.0 1.0...
                    1.0 1.0 1.0 1.0 1.0 1.0...
                    1.0 1.0 1.0 1.0 1.0 1.0 ];
ndt=10;                 % Number of slices used in averageing between t1 and t2
tunits_press=0.001;     % Thomson Time units (ms)
tunits_mse=1.0;         % MSE Time units (s)
ne_units=1e16;          % Electron Number Density Units
ec=1.60217653e-19;      % NIST 2008 pg 34
phi_press=-18.0;        % Phi Location of Pressure Profile (degrees)
Runits_press=0.001;     % units/meter for Radial Pressure Profile. [mm/m]
z_press=0.0;            % Z Location of Pressure Profile
scale_curtor=1000.;     % Toroidal Current is in [kA]
scale_extcur=1000.;     % Coil Currents are in [kA]
scale_eplasma=1000.;    % Scaling Factor for total stored energy [kJ]
scale_fdia=1.0e-4;      % Diamagnetic loop in (10^-4) [Wb]
mu0=4*pi*10^-7;         % Mu0
use_mags=0;             % Switch to use DIAGNO
luse_mse=0;             % Switch to use MSE data
lmse_vac=0;             % Switch to run a single iteration to calculate VMEC MSE vacuum response
luse_pres_spline=0;     % Switch to use pressure splines (v8.47)
luse_curr_spline=0;     % Switch to use current splines (v8.47)
luse_iota_spline=0;     % Switch to use iota splines (v8.47)
lcurtor=0;              % Switch to target total toroidal current
leplasma=0;             % Switch to target total stored energy
ldiamag=0;              % Switch to use diagno for diamagnetic loop
liota_opt=0;            % Switch to optimize iota instead of I'
lisote=0;
lcurr_opt=0;
npspline=15;            % Number of datapoints in Pressure spline
njspline=15;            % Number of datapoints in Current spline
nnodes=0;                % Node counter
fluxloop_legtext={'SL0011' 'SL0012' 'SL0013'...
                  'SL0021' 'SL0022' 'SL0023'...
                  'SL0031' 'SL0032' 'SL0033'...
                  'SL0041' 'SL0042' 'SL0043'...
                  'SL0051' 'SL0052' 'SL0053'...
                  'SL0061' 'SL0062' 'SL0063'...
                  'SL0071' 'SL0072' 'SL0073'...
                  'SL0081' 'SL0082' 'SL0083'};
% Handle varargin
numdefargs=3;   %Number of default arguments
if nargin >numdefargs
    i=1;
    while i<=(nargin-numdefargs)
        if isstruct(varargin{i})
            if isfield(varargin{i},'datatype')
                switch varargin{i}.datatype
                    case {'coil_data'}
                        coil_data=varargin{i};
                    case {'vessel'}
                        vesdata=varargin{i};
                    case {'wout'}
                        vmec_data=varargin{i};
                    case {'stellopt_mse'}
                        mse_vac=varargin{i};
                end
            end
        elseif strcmp(varargin{i},'shotnum')
            i=i+1;
            shotnum=varargin{i};
            filename=['input.' shotnum '_'];
        elseif strcmp(varargin{i},'user')
            i=i+1;
            user=varargin{i};
        else
            switch varargin{i}
                case {'recon','VMEC'}
                    outputtype=varargin{i};
                case {'plots'}
                    makeplots=1;
                case {'magdiags'}
                    use_mags=1;
                case {'msedata'}
                    luse_mse=1;
                case {'msevac'}
                    lmse_vac=1;
                    luse_mse=1;
                case {'curtor'}
                    lcurtor=1;
                case {'pres_splines'}
                    luse_pres_spline=1;
                case {'curr_splines'}
                    luse_curr_spline=1;
                    lcurr_opt=1;
                case {'iota_splines'}
                    luse_iota_spline=1;
                    liota_opt=1;
                case {'eplasma'}
                    leplasma=1;
                case {'iota'}
                    liota_opt=1;
                case {'isote'}
                    lisote=1;
                case {'diamag'}
                    ldiamag=1;
            end
        end
        i=i+1;
    end
end
% Check for diagno.control file
diagno_file=dir('diagno.control');
if isempty(diagno_file)
    lwrite_diagno=1;
else
    lwrite_diagno=0;
end
% Begin Program
t=(0:data.tdim-1).*data.dt+data.t0;     % Time Vector
% -----  INDATA NAMELIST VARS  ----- 
%mgrid_file='/u/slazerso/Sims/LHD/coils/lhd_ys4/mgrid.lhd_ys4.msize';
mgrid_file='/p/lhd/coils/lhd_ys4/mgrid_lhd_ys4_847_nzeta32.nc';
nzeta=32;
delt=0.9;
tcon0=1.0;
niter=2500;
nstep=200;
nvacskip=6;
gamma=0.0;
nfp=10;
mpol=8;
ntor=6;
ncurr=1;
bloat=1.0;
spres_ped=1.0;
ns_array=[9 29 49 99];
ftol_array=[1.0e-6 1.0e-8 1.0e-10 1.0e-12];
am=     [0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00];
ai=     [0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00];
ac=     [0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00];
raxis=  [3.80 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00];
zaxis=  [0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00];
am_aux_s=[];
am_aux_f=[];
ac_aux_s=[];
ac_aux_f=[];
ai_aux_s=[];
ai_aux_f=[];
pmass_type=[];
pcurr_type=[];
piota_type=[];
niter=20000;
rbc00=   3.70;
rbc01=   0.5;
rbc11=  -0.13;
rbcn11=  0.00;
zbs00=   0.0;
zbs01=  -0.5;
zbsn11=  0.00;
zbs11=  -0.14;
phiedge= 1.35;
ac(1)=1;
ac(2)=-1;
ac(3)=0;
ac(4)=0;
ac(5)=-1;
ac(6)=1;
% Make some default parameters
ns=max(ns_array);
npoly=11;
% -----  OPTIMUM NAMELIST VARS  -----
%        Run Control
epsfcn = 1e-6;
niter_opt = 5000;
num_processors = 16;
num_levmar_params = 16;
nopt_alg = 0;
nopt_boundary = 0;
lreset_opt = 1;
ldiag_opt = 0;
lkeep_mins = 1;         %Keep all minimum outputs
%        Physics Modules
lcur_prof_opt=0;
lcur_opt_edge0=0;
liota_prof_opt=0;
ldiagno_opt=0;
lpres_prof_opt=1;
lpres_opt_edge=1;
lpres_opt_edgegr0=0;
%        Scalar Optimizations
target_beta=0.0;                sigma_beta=bigno;
target_curtor=0.0;              sigma_curtor=bigno;
target_eplasma=1.0;             sigma_eplasma=bigno;
target_iota_max = 2.0;         sigma_iota_max = bigno;
target_iota_min = 0.0;         sigma_iota_min =bigno;
target_iota_p = zeros(1,11);
sigma_iota_pmax = 0.5.*ones(1,max(ns_array));
%        Vacuum Vessel Matching Optimizations
target_vv=1.0;                  sigma_vv=bigno;
lvv_tgt_min=0;
target_vv_rms=1.0;              sigma_vv_rms=bigno;
                                sigma_vv_max=bigno;
mpol_vv=0;
ntor_vv=0;
rbc_vv=zeros(2*ntor_vv+1,mpol_vv+1);
zbs_vv=zeros(2*ntor_vv+1,mpol_vv+1);
%        Boundary Matching Optimizations
target_bd=1.0;                  sigma_bd=bigno;
target_bd_rms=1.0;              sigma_bd_rms=bigno;
                                sigma_bd_max=bigno;
mpol_bd=2;
ntor_bd=2;
rbc_bd=zeros(2*ntor_vv+1,mpol_vv+1);
zbs_bd=zeros(2*ntor_vv+1,mpol_vv+1);
%        Limiter Specification
phi_lim=zeros(1,4); %Max 60 elements
r_lim=zeros(2,4);   %max 40,60 elements
z_lim=zeros(2,4);   %max 40,60 elements
%        Boundary Targets
shapeweight=1;
target_rbtor=1.0;               sigma_rbtor=bigno;
planes_bw=[0 pi/2. pi];
amplw_bw=[1 1 1];
theta0_bw=[0 pi/2. pi];
phi0_bw=pi;
wtheta_bw=pi/8.;
wphi_bw=pi/8.;
% Pressure profile optimizations
np_prof=100;
factor_p_prof=1.000;
lp_prof_incl_edge = 1;
ne_prof=zeros(1,np_prof);
te_prof=zeros(1,np_prof);
sigma_p_prof=zeros(1,np_prof);
r_p_prof=zeros(1,np_prof);
z_p_prof=zeros(1,np_prof);
phi_p_prof=zeros(1,np_prof);
%        DIAGNO Optimizations
diagno_control='diagno.control';
ndiagno_seg=0;
if (ndiagno_seg > 0)
    target_diagno_seq=zeros(1,ndiagno_seg);
    sigma_diagno_seq=zeros(1,ndiagno_seg);
end
ndiagno_flx=24;
if (ndiagno_flx > 0)
    target_diagno_flx=zeros(1,ndiagno_flx);
end

% Handle Vessel data
target_vv=0.0;
if ~isempty(vesdata)
    target_vv=1.0;
    mpol_vv=0;
    ntor_vv=0;
    target_vv=0;
    sigma_vv=1e-2;
    % Now setup limiter;
    nphi_lim=max(vesdata.coords(:,4));
    phi_lim=1:nphi_lim;
    nr_lim=max(find(vesdata.coords(:,4)==1,'first'));
    r_lim=zeros(nr_lim-1,nphi_lim);
    z_lim=zeros(nr_lim-1,nphi_lim);
    for i=1:nphi_lim
        phi_lim(i)=vesdata.coords(find(vesdata.coords(:,4)==i,1,'first'),3);
        for j=1:nr_lim-1
            k=j+(i-1)*nr_lim;
            r_lim(j,i)=vesdata.coords(k,1);
            z_lim(j,i)=vesdata.coords(k,2);
        end
    end
    phi_lim=phi_lim.*180./pi;
end

% For total run
ldiagno_opt=0;
lpres_prof_opt=1;
lpres_opt_edge=1;
lcur_prof_opt=0;
lcur_opt_edge0=0;
pres_opt_nmax=11;

% So Dave's script reads in Wp then calculates fdia by
%  fdia = 0.39456*wp/1.5
%  Wp is the reported value so we go back to Wp
wp=1.5.*data.fdia./0.39456;
% For 85384
%Rax = 3.6;  Bt=1.5;
% For 82716
Rax = 3.75;  Bt=1.3;
% OK so this is for shot 85384 but
% Rax = 3.6
% Bt  = -1.5
% So Wp=3*pi*Bt*R*dFdia/mu0
% Then fdia = mu0*wp / (3*pi*Bt*R)
%fdia2= mu0*wp*scale_eplasma / (3*pi*1.5*3.6);
fdia2= mu0*wp*scale_eplasma / (3*pi*Bt*Rax);  % May not need to scale eplasma
% The paramagnetic flux is fpara = (mu0*I)^2/(8*pi*Bt)
fpara= (mu0.*data.Ip.*scale_curtor).*(mu0.*data.Ip.*scale_curtor)./(8*pi*Bt);
% The helical flux is fhelic = k*mu0*I*phi0/(pi*R*Bt)
fdia3= fdia2+fpara;

% Get the rogowski signal


% Now output some shot stuff
% Plot the coil currents
if makeplots
    LHDplot_current(data,'shot',shotnum);
    saveas(gca,['coil_currents_' shotnum '.fig']);
    pause(.01);
end

% Plot Thomson
if isfield(data,'T_e') && makeplots
    % Plot the Thompson Data
    pixplot(double(data.R_TS).*Runits_press,double(data.t_TS).*tunits_press,double(data.T_e')./1000);
    xlabel('R [m]');
    ylabel('Time [s]');
    title(['Electron Temperature (TS) (shot: ' shotnum ')']);
    hc=colorbar;
    ylabel(hc,'[keV]')
    saveas(gca,['thomson_te_' shotnum '.fig']);
    pause(.01);
    pixplot(double(data.R_TS).*Runits_press,double(data.t_TS).*tunits_press,double(data.n_e')./1000);
    xlabel('R [m]');
    ylabel('Time [s]');
    title(['Electron Density (TS) (shot: ' shotnum ')']);
    hc=colorbar;
    ylabel(hc,'x10^{19} [m^{-3}]');
    saveas(gca,['thomson_ne_' shotnum '.fig']);
    pause(.01);
    pixplot(double(data.R_TS).*Runits_press,double(data.t_TS).*tunits_press,double(data.n_e').*double(data.T_e').*ec.*ne_units.*.001);
    xlabel('R [m]');
    ylabel('Time [s]');
    title(['Electron Pressure (TS) (shot: ' shotnum ')']);
    hc=colorbar;
    ylabel(hc,'[kPa]');
    saveas(gca,['thomson_pe_' shotnum '.fig']);
    pause(.01);
end
% Plot MSE
if isfield(data,'Pol_Angle') && makeplots
    % Plot the MSE Data
    pixplot(double(data.R_MSE),double(data.t_MSE).*tunits_mse,double(data.Pol_Angle'));
    xlabel('R [m]');
    ylabel('Time [s]');
    title(['Polarization Angle (MSE) (shot: ' shotnum ')']);
    hc=colorbar;
    ylabel(hc,'[degrees]')
    saveas(gca,['mse_' shotnum '.fig']);
    pause(.01);
end

% Plot the Fluxloop Data
if isfield(data,'FluxLoop0') && makeplots
    set(gcf,'Position',[1 1 1024 768]);
    plot(t,data.('FluxLoop0'),'Color','k','LineStyle','-')
    hold on
    for i=2:ndiagno_flx
        plot(t,data.(['FluxLoop' num2str(i-1)]),'Color','k','LineStyle','-')
    end
    xlabel('Time [s]');
    ylabel('\delta Flux [Wb]');
    title(['Flux Loops (shot: ' shotnum ')']);
    legend(fluxloop_legtext);
    saveas(gca,['delta_fluxloops_' shotnum '.fig']);
    pause(.01);
end

% Plot the toroidal current Data
if isfield(data,'Ip') && makeplots
    clf;
    set(gcf,'Position',[1 1 1024 768]);
    plot(t,data.Ip,'Color','k','Linestyle','-');
    xlabel('Time [s]');
    ylabel('Current [kA]');
    title(['Toroidal Current (shot: ' shotnum ')']);
    saveas(gca,['curtor_' shotnum '.fig']);
    pause(.01);
end

% Plot the diamagnetic flux
if isfield(data,'fdia') && makeplots
    clf;
    set(gcf,'Position',[1 1 1024 768]);
    plot(t,data.fdia,'Color','k','Linestyle','-');
    xlabel('Time [s]');
    ylabel(' Flux 10^{-4} [Wb]');
    title(['Diamagnetic Flux (shot: ' shotnum ')']);
    saveas(gca,['diamag_' shotnum '.fig']);
    pause(.01);
end

% Plot total stored energy
if isfield(data,'fdia') && makeplots
    clf;
    set(gcf,'Position',[1 1 1024 768]);
    plot(t,wp,'Color','k','Linestyle','-');
    xlabel('Time [s]');
    ylabel('Energy [J]');
    title(['Stored Energy (shot: ' shotnum ')']);
    saveas(gca,['eplasma_' shotnum '.fig']);
    pause(.01);
end

% First we convert dex to timslice data for file names
% Output Table Headings
disp('  TIME     I(HCI)     I(HCM)      I(HCO)       I(POV)       I(PIS)     I(PIV)   ITOR [kA]  EPLASMA [kJ]');
n=numel(timedex1);
for i=1:n
    % Find 
    % Create input.val filename
    timestamp=num2str(int32(timedex1(i)*1000),perct);
    % Create EXTCUR ARRAY
    t=(0:data.tdim-1).*data.dt+data.t0;     % Time Vector
    t_temp=timedex1(i):(timedex2(i)-timedex1(i))/(ndt-1):timedex2(i);
    hci=pchip(t,data.CurrentHCI2);
    hcm=pchip(t,data.CurrentHCM0);
    hco=pchip(t,data.CurrentHCO1);
    pov=pchip(t,data.CurrentPOV3);
    piv=pchip(t,data.CurrentPIV5);
    pis=pchip(t,data.CurrentPIS4);
    extcur=[mean(ppval(hci,t_temp)) mean(ppval(hcm,t_temp)) mean(ppval(hco,t_temp))...
        -mean(ppval(pov,t_temp)) -mean(ppval(pis,t_temp)) mean(ppval(piv,t_temp))].*scale_extcur;
    Ip_spline=pchip(t,data.Ip*scale_curtor);
    % Handle reconstructions
    if strcmp(outputtype,'recon')
        % Create Thomson Splines
        ne_spline=pchip(data.t_TS*tunits_press,permute(data.n_e.*ne_units,[2 1]));
        te_spline=pchip(data.t_TS*tunits_press,permute(data.T_e,[2 1]));
        % Calculate Thomson Pressure
        ne_prof=mean(permute(ppval(ne_spline,t_temp), [2 1]));
        te_prof=mean(permute(ppval(te_spline,t_temp), [2 1]));
        p_prof=ne_prof.*ec.*te_prof;
        % Calculate Standard Deviations
        % Note the we define sigma_p_prof by p=n*kB*(T_e+T_i+T_h) so
        % sigma_p=dn*kB*(T_e+T_i+T_h)+n*kB*(dT_e+dT_i+dT_h)
        std_ne=std(permute(ppval(ne_spline,t_temp), [2 1]),1);
        std_te=std(permute(ppval(te_spline,t_temp), [2 1]),1);
        % Set STELLOPT pressure
        np_prof       = max(size(ne_prof));
        r_p_prof      = data.R_TS.*Runits_press;
        z_p_prof      = z_press*ones(1,np_prof);
        phi_p_prof    = phi_press.*ones(1,np_prof);
        factor_p_prof = 3.0;
        sigma_p_prof  = (std_ne.*ec.*te_prof+ne_prof.*ec.*std_te);
        % Force sigma p_prof to be at least 15%
        percent=sigma_p_prof./p_prof;
        percent(p_prof < 1000) = 1.;
        sigma_p_prof(percent < 0.15) = p_prof(percent < 0.15).*0.15;
        %sigma_p_prof=sigma_p_prof.*sqrt(np_prof);
        % Set Pressure profile constraints
        sigma_pedge   = 1000.;  % This could be improved
        sigma_pgrad   = 0.1*ones(1,ns_array(length(ns_array)));  % This could be improved
        % Clean Pressure and set sigma.
        ne_cutoff=max(ne_prof)*.02;
        te_cutoff=max(te_prof)*.02;
        p_cutoff=max(p_prof)*.02;
        sigma_p_zero=p_cutoff;
        ne_prof(ne_prof < ne_cutoff)=0.0;
        te_prof(te_prof < te_cutoff)=0.0;
        p_prof(p_prof < p_cutoff)=0.0;
        ne_prof(p_prof < p_cutoff)=0.0;
        te_prof(p_prof < p_cutoff)=0.0;
        sigma_p_prof(ne_prof < ne_cutoff)=sigma_p_zero;
        sigma_p_prof(te_prof < te_cutoff)=sigma_p_zero;
        sigma_p_prof(p_prof < p_cutoff)=sigma_p_zero;
        % Check for erroneous points
        for j=2:np_prof
            if (p_prof(j) == 0.0) && (p_prof(j-1) > 0.0) && (p_prof(j+1)>0.0)
                sigma_p_prof(j) = 1e30;
            end
            if (te_prof(j) == 0.0) && (te_prof(j-1) > 0.0) && (te_prof(j+1)>0.0)
                sigma_p_prof(j) = 1e30;
            end
            if (ne_prof(j) == 0.0) && (ne_prof(j-1) > 0.0) && (ne_prof(j+1)>0.0)
                sigma_p_prof(j) = 1e30;
            end
        end
        % Set VMEC Pressure
        gamma=0.0;
        am(1:11)=0.0;
        am(1)=max(ne_prof)*ec*max(te_prof)*factor_p_prof;
        am(2)=-am(1);
        am(3)=0.0;
        am(4)=0.0;
        am(5)=-am(1);
        am(6)=am(1);
        nnodes = nnodes + 11; % Assume polynomial
        lphiedge=1;  nnodes = nnodes + 1;
        % Extract Segmented Rogowski values
        for j=1:ndiagno_seg
            % Spline over data
            seg_spline=pchip(t,data.(['FluxLoop' num2str(j-1)]));
            target_diagno_seg(j)=mean(ppval(seg_spline,t_temp));
            sigma_diagno_seg(i)=std(ppval(seg_spline,t_temp),1);
        end
        if (ndiagno_seg > 0)
            percent=sigma_diagno_seg./target_diagno_seg;
            sigma_diagno_seg(percent < 0.15)=0.15.*target_diagno_seg(percent < 0.15);
        end
        % Extract Flux values (use delta flux if no vacdata supplied)
        sigma_diagno_flx=sigma_diagno_flx_default;
        for j=1:ndiagno_flx
            flux_spline=pchip(t,data.(['FluxLoop' num2str(j-1)]));
            target_diagno_flx(j)=mean(ppval(flux_spline,t_temp));
            if (sigma_diagno_flx(j) < bigno)
                sigma_diagno_flx(j)=std(ppval(flux_spline,t_temp),1)*3;
            end
        end
        if (ndiagno_flx > 0)
            percent=sigma_diagno_flx./target_diagno_flx;
            sigma_diagno_flx(percent < 0.15)=0.15.*target_diagno_flx(percent < 0.15);
        end
        % Extract Toroidal Current
        if lcurtor
            target_curtor=mean(ppval(Ip_spline,t_temp));
            sigma_curtor=std(ppval(Ip_spline,t_temp),1)*3;
            if (sigma_curtor/target_curtor < 0.05)
                sigma_curtor = 0.05*target_curtor;
            end
        end
        % Extract Total Stored Energy
        if leplasma
            eplasma_spline=pchip(t,wp*scale_eplasma);
            target_eplasma=mean(ppval(eplasma_spline,t_temp));
            sigma_eplasma=std(ppval(eplasma_spline,t_temp),1)*3;
            if (sigma_eplasma/target_eplasma < 0.40)
                sigma_eplasma = 0.40*target_eplasma;
            end
        end
        if ldiamag
            lwrite_diagno = 1;
            %diamag_spline=pchip(t,data.fdia*scale_fdia);
            diamag_spline=pchip(t,fdia2);
            target_diamag=mean(ppval(diamag_spline,t_temp));
            sigma_diamag=std(ppval(diamag_spline,t_temp),1)*3;
            if (sigma_diamag/target_diamag < 0.05)
                sigma_diamag = 0.05*target_diamag;
            end
        end
        % MSE Section
        if luse_mse && isfield(data,'Pol_Angle') && ~isempty(coil_data) && ~lmse_vac && isempty(mse_vac)
            pol_spline=pchip(data.t_MSE*tunits_mse,permute(pi.*data.Pol_Angle./180,[2 1]));
            dpol_spline=pchip(data.t_MSE*tunits_mse,permute(pi.*(abs(data.dPol_Angle)+1)./180,[2 1]));  % Add 1 degree uncertainty from discrepancy between vac MSE and vac VMEC
            nmse_chords=length(data.R_MSE);
            nmse_cams=1; % LHD has one camera
            coil_data=coil_biot_prep(coil_data);
            mse_r=zeros(nmse_cams,nmse_chords);
            mse_phi=zeros(nmse_cams,nmse_chords);
            mse_z=zeros(nmse_cams,nmse_chords);
            mse_alpha=zeros(nmse_cams,nmse_chords);
            mse_beta=zeros(nmse_cams,nmse_chords);
            mse_theta=zeros(nmse_cams,nmse_chords);
            mse_pol_vac=zeros(nmse_cams,nmse_chords);
            mse_pol=zeros(nmse_cams,nmse_chords);
            sigma_mse_pol=zeros(nmse_cams,nmse_chords);
            %fid=fopen('mse_test.txt','w');
            for j=1:nmse_cams
                mse_r(j,:)=data.R_MSE(:);
                mse_phi(j,:)=pi*(data.phi_MSE(:)+108)/180; %Note all angles in degrees (LHD assumes MSE is in PHI angles = CAD - 105)
                mse_z(j,:)=0.0*data.R_MSE(:);        %In Z=0 plane for LHD
                mse_alpha(j,:)=pi*data.alpha_MSE(:)/180;
                mse_beta(j,:)= pi*data.beta_MSE(:)/180;
                mse_theta(j,:)=pi*data.THETA_MSE(:)/180;
                % Now calculate the vacuum field
                for k=1:nmse_chords
                    %fprintf(fid,'%d',k);
                    x=mse_r(j,k)*cos(mse_phi(j,k));
                    y=mse_r(j,k)*sin(mse_phi(j,k));
                    z=mse_z(j,k);
                    [bx by bz]=coil_biot(coil_data,x,y,z,extcur);
                    br=bx*cos(mse_phi(j,k))+by*sin(mse_phi(j,k));
                    bphi=by*cos(mse_phi(j,k))-bx*sin(mse_phi(j,k));
                    mse_pol_vac(j,k)=atan(...
                                          bz*cos(mse_theta(j,k)) /...
                                         ( bz*sin(mse_beta(j,k))*sin(mse_theta(j,k))...
                                    +( br*cos(mse_alpha(j,k))+bphi*sin(mse_alpha(j,k)))*cos(mse_beta(j,k)) ) );
                    %fprintf(fid,'%20.10e',mse_r(j,k),mse_phi(j,k),mse_z(j,k),...
                    %    br,bphi,bz,bx,by,bz,mse_pol_vac(j,k));
                    %fprintf(fid,'\n');
                end
                mse_pol(j,:)=mean(permute(ppval(pol_spline,t_temp), [2 1]));
                sigma_mse_pol(j,:)=mean(permute(ppval(dpol_spline,t_temp), [2 1]));
            end
            %fclose(fid);
        elseif luse_mse && isfield(data,'Pol_Angle') && ~isempty(mse_vac)
            pol_spline=pchip(data.t_MSE*tunits_mse,permute(pi.*data.Pol_Angle./180,[2 1]));
            dpol_spline=pchip(data.t_MSE*tunits_mse,permute(pi.*(abs(data.dPol_Angle)+1)./180,[2 1]));  % Add 1 degree uncertainty from discrepancy between vac MSE and vac VMEC
            nmse_chords=length(data.R_MSE);
            nmse_cams=1; % LHD has one camera
            mse_r=zeros(nmse_cams,nmse_chords);
            mse_phi=zeros(nmse_cams,nmse_chords);
            mse_z=zeros(nmse_cams,nmse_chords);
            mse_alpha=zeros(nmse_cams,nmse_chords);
            mse_beta=zeros(nmse_cams,nmse_chords);
            mse_theta=zeros(nmse_cams,nmse_chords);
            mse_pol_vac=zeros(nmse_cams,nmse_chords);
            mse_pol=zeros(nmse_cams,nmse_chords);
            sigma_mse_pol=zeros(nmse_cams,nmse_chords);
            for j=1:nmse_cams
                mse_r(j,:)=data.R_MSE(:);
                mse_phi(j,:)=pi*(data.phi_MSE(:)+108)/180; %Note all angles in degrees (LHD assumes MSE is in PHI angles = CAD - 105)
                mse_z(j,:)=0.0*data.R_MSE(:);        %In Z=0 plane for LHD
                mse_alpha(j,:)=pi*data.alpha_MSE(:)/180;
                mse_beta(j,:)= pi*data.beta_MSE(:)/180;
                mse_theta(j,:)=pi*data.THETA_MSE(:)/180;
                %fid=fopen('mse_test.txt','w');
                % Now calculate the vacuum field
                for k=1:nmse_chords
                    if (mse_vac.data(k,6) == 0.0)
                        mse_pol_vac(j,k) = mse_vac.data(k,7);
                    else
                        disp('ERROR: mse_vac not zero in supplied mse_data structure');
                        return
                    end
                end
                mse_pol(j,:)=mean(permute(ppval(pol_spline,t_temp), [2 1]));
                sigma_mse_pol(j,:)=mean(permute(ppval(dpol_spline,t_temp), [2 1]));
            end
        elseif luse_mse && isfield(data,'Pol_Angle') && lmse_vac && isempty(mse_vac)
            pol_spline=pchip(data.t_MSE*tunits_mse,permute(pi.*data.Pol_Angle./180,[2 1]));
            dpol_spline=pchip(data.t_MSE*tunits_mse,permute(pi.*(abs(data.dPol_Angle)+1)./180,[2 1]));  % Add 1 degree uncertainty from discrepancy between vac MSE and vac VMEC
            nmse_chords=length(data.R_MSE);
            nmse_cams=1; % LHD has one camera
            mse_r=zeros(nmse_cams,nmse_chords);
            mse_phi=zeros(nmse_cams,nmse_chords);
            mse_z=zeros(nmse_cams,nmse_chords);
            mse_alpha=zeros(nmse_cams,nmse_chords);
            mse_beta=zeros(nmse_cams,nmse_chords);
            mse_theta=zeros(nmse_cams,nmse_chords);
            mse_pol_vac=zeros(nmse_cams,nmse_chords);
            mse_pol=zeros(nmse_cams,nmse_chords);
            sigma_mse_pol=zeros(nmse_cams,nmse_chords);
            for j=1:nmse_cams
                mse_r(j,:)=data.R_MSE(:);
                mse_phi(j,:)=pi*(data.phi_MSE(:)+108)/180; %Note all angles in degrees (LHD assumes MSE is in PHI angles = CAD - 105)
                mse_z(j,:)=0.0*data.R_MSE(:);        %In Z=0 plane for LHD
                mse_alpha(j,:)=pi*data.alpha_MSE(:)/180;
                mse_beta(j,:)= pi*data.beta_MSE(:)/180;
                mse_theta(j,:)=pi*data.THETA_MSE(:)/180;
                %fid=fopen('mse_test.txt','w');
                % Now calculate the vacuum field
                for k=1:nmse_chords
                    mse_pol_vac(j,k) = 0.0;
                end
                mse_pol(j,:)=mean(permute(ppval(pol_spline,t_temp), [2 1]));
                sigma_mse_pol(j,:)=mean(permute(ppval(dpol_spline,t_temp), [2 1]));
            end
            luse_pres_spline = 0;
            luse_curr_spline = 0;
            lpres_prof_opt = 0;
            lpres_opt_edge = 0;
            lpres_opt_edgegr0 = 0;
            lp_prof_incl_edge = 0;
            lcurr_opt = 0;
            lcur_prof_opt = 0;
            ldiagno_opt = 0;
            ldiamag = 0;
            leplasma = 0;
            lcurtor = 0;
            curtor= 0.0;
            ac_aux_s = [];
            am_aux_s = [];
            ai_aux_s = [];
            ac = 0.0*(1:11);
            am = 0.0*(1:11);
            ai = 0.0*(1:11);
            niter_opt = 1;
        elseif luse_mse && isfield(data,'Pol_Angle') && isempty(coil_data) && ~isempty(vmec_data)
            pol_spline=pchip(data.t_MSE*tunits_mse,permute(pi.*data.Pol_Angle./180,[2 1]));
            dpol_spline=pchip(data.t_MSE*tunits_mse,permute(pi.*(abs(data.dPol_Angle)+1)./180,[2 1]));  % Add 1 degree uncertainty from discrepancy between vac MSE and vac VMEC
            nmse_chords=length(data.R_MSE);
            nmse_cams=1; % LHD has one camera
            coil_data=coil_biot_prep(coil_data);
            mse_r=zeros(nmse_cams,nmse_chords);
            mse_phi=zeros(nmse_cams,nmse_chords);
            mse_z=zeros(nmse_cams,nmse_chords);
            mse_alpha=zeros(nmse_cams,nmse_chords);
            mse_beta=zeros(nmse_cams,nmse_chords);
            mse_theta=zeros(nmse_cams,nmse_chords);
            mse_pol_vac=zeros(nmse_cams,nmse_chords);
            mse_pol=zeros(nmse_cams,nmse_chords);
            sigma_mse_pol=zeros(nmse_cams,nmse_chords);
            for j=1:nmse_cams
                mse_r(j,:)=data.R_MSE(:);
                mse_phi(j,:)=pi*(data.phi_MSE(:)+108)/180; %Note all angles in degrees (LHD assumes MSE is in PHI angles = CAD - 105)
                mse_z(j,:)=0.0*data.R_MSE(:);        %In Z=0 plane for LHD
                mse_alpha(j,:)=pi*data.alpha_MSE(:)/180;
                mse_beta(j,:)= pi*data.beta_MSE(:)/180;
                mse_theta(j,:)=pi*data.THETA_MSE(:)/180;
                % Now calculate the vacuum field
                theta_vmec = 0:2*pi/359:2*pi;
                for k=1:nmse_chords
                    % First Calculate Fields
                    zeta_vmec = mse_phi(j,k);
                    r=cfunct(theta_vmec,zeta_vmec,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
                    bu=cfunct(theta_vmec,zeta_vmec,vmec_data.bsupumnc,vmec_data.xm,vmec_data.xn);
                    bv=cfunct(theta_vmec,zeta_vmec,vmec_data.bsupvmnc,vmec_data.xm,vmec_data.xn);
                    drdu=sfunct(theta_vmec,zeta_vmec,vmec_data.rumns,vmec_data.xm,vmec_data.xn);
                    drdv=sfunct(theta_vmec,zeta_vmec,vmec_data.rvmns,vmec_data.xm,vmec_data.xn);
                    dzdu=cfunct(theta_vmec,zeta_vmec,vmec_data.zumnc,vmec_data.xm,vmec_data.xn);
                    dzdv=cfunct(theta_vmec,zeta_vmec,vmec_data.zvmnc,vmec_data.xm,vmec_data.xn);
                    br_vmec = bu.*drdu+bv.*drdv;
                    bphi_vmec = r.*bv;
                    bz_vmec = bu.*dzdu+bv.*dzdv;
                    % Now spline
                    % Now calculate MSE response
                    mse_pol_vac(j,k)=atan(...
                                          bz*cos(mse_theta(j,k)) /...
                                         ( bz*sin(mse_beta(j,k))*sin(mse_theta(j,k))...
                                    +( br*cos(mse_alpha(j,k))+bphi*sin(mse_alpha(j,k)))*cos(mse_beta(j,k)) ) );
                    %fprintf(fid,'%d %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n',k,mse_pol_vac(j,k),br,bphi,bz,bx,by,bz);
                end
                mse_pol(j,:)=mean(permute(ppval(pol_spline,t_temp), [2 1]));
                sigma_mse_pol(j,:)=mean(permute(ppval(dpol_spline,t_temp), [2 1]));
            end
        elseif luse_mse && isfield(data,'Pol_Angle') && isempty(coil_data) && isempty(vmec_data)
            disp('ERROR: Either WOUT or COIL_DATA must be supplied to calculated vacuum MSE response.');
            luse_mse = 0;
        elseif luse_mse && ~isfield(data,'Pol_Angle')
            disp('ERROR: No MSE data found in shot data structure.');
            luse_mse = 0;
        end
        % Handle Currnet Profile Optimization
        if lcurr_opt
            lcur_prof_opt=1;
            nnodes = nnodes + 11; % Assume polynomial
            nnodes = nnodes + 1;  % Assume curtor
        end
        % Handle v8.47 spline values
        if luse_pres_spline
            nnodes = nnodes - 11 + npspline;
            %am_aux_s=(0:1/(npspline-1):1).^2;
            am_aux_s=(0:1/(npspline-1):1);
            am_aux_f=polyval(fliplr(am),am_aux_s);
            am_aux_f(length(am_aux_f))=0.0;  % Default edge to zero
            pmass_type='Akima_spline';
        end
        if liota_opt && ~isempty(vmec_data)
            liota_prof_opt=1;
            ncurr = 0;
            curtor = 0; % This is important
            if luse_iota_spline
                %ai_aux_s=(0:1/(njspline-1):1).^2;
                %ai_aux_s=[0.0  0.01  0.1  0.25  0.5  0.75  0.8  0.9  1.0];
                ai_aux_s=(0:1/(njspline-1):1);
                nnodes = nnodes + njspline;
                ai_aux_f=pchip(0:1/(vmec_data.ns-1):1,vmec_data.iotaf,ai_aux_s);
                piota_type='Akima_spline';
            else
                nnodes = nnodes + 11;
                ai=fliplr(polyfit(0:1/(vmec_data.ns-1):1,vmec_data.iotaf,11));
            end
            target_iota_p = zeros(1,11);
            sigma_iota_pmax = 0.5.*ones(1,max(ns_array));
            sigma_iota_min = 0.1;
            sigma_iota_max = 0.1;
            pcurr_type = 'power_series';
        elseif liota_opt
            liota_prof_opt=1;
            ncurr = 0;
            curtor = 0; % This is important
            ai(1) = 0.4;
            ai(2) = 0.0;
            ai(3) = 1.0;
            if luse_iota_spline
                %ai_aux_s=(0:1/(njspline-1):1).^2;
                %ai_aux_s=[0.0  0.01  0.1  0.25  0.5  0.75  0.8  0.9  1.0];
                ai_aux_s=(0:1/(njspline-1):1);
                nnodes = nnodes + njspline;
                ai_aux_f = polyval(fliplr(ai),am_aux_s);
                piota_type='Akima_spline';
            else
                nnodes = nnodes + 11;
                ai(1) = 0.4;
                ai(2) = 1e-2;
                ai(3) = 1.0;
                ai(4:11) = 0.0;
                piota_type='power_series';
            end 
            target_iota_p = zeros(1,11);
            sigma_iota_pmax = 0.5.*ones(1,max(ns_array));
            sigma_iota_min = 0.1;
            sigma_iota_max = 0.1;
            pcurr_type = 'power_series';
        elseif ~liota_opt
            % AC array is defaulted
            pcurr_type = 'power_series';
            piota_type = 'power_series';
            ncurr = 1;
            curtor=mean(ppval(Ip_spline,t_temp));
            if luse_curr_spline
                nnodes = nnodes - 11 + njspline;
                %ac_aux_s=[ 0.0  0.16  0.25  0.36  0.49  0.64  0.81  1.0];
                %ac_aux_s=(0:1/(njspline-1):1).^2;
                ac_aux_s=(0:1/(njspline-1):1);
                ac_aux_f=polyval(fliplr(ac),ac_aux_s);
                ac_aux_f(length(ac_aux_f))=0.0;  % Default edge to zero
                pcurr_type='Akima_spline_Ip';
            end
        else
            disp('ERROR: VMEC output file needed for Iota Optimization');
            return
        end
            
    end
    % Creat filename extenstion
    file_ext='';
    if strcmp(outputtype,'recon')
        if lpres_prof_opt, 
            if lisote
                file_ext=[file_ext 't'];
            else
                file_ext=[file_ext 'p'];
            end
        if lcur_prof_opt, file_ext=[file_ext 'j']; end
        if liota_prof_opt, file_ext=[file_ext 'i']; end
        if ldiagno_opt, file_ext=[file_ext 'm']; end
        if ldiamag, file_ext=[file_ext 'd']; end
        if leplasma, file_ext=[file_ext 'e']; end
        if lcurtor, file_ext=[file_ext 'I']; end
        if luse_mse, file_ext=[file_ext 'M']; end
        file_ext=[file_ext '_'];
        if lmse_vac, file_ext=[file_ext 'vac']; end % Do this here because no amac tags will be printed for MSE_VAC
        if lpres_prof_opt
            if isempty(am_aux_s)
                file_ext=[file_ext 'am'];
            else
                file_ext=[file_ext 'sm'];
            end
        end
        if lcur_prof_opt
            if isempty(ac_aux_s)
                file_ext=[file_ext 'ac'];
            else
                file_ext=[file_ext 'sc'];
            end
        end
        if liota_prof_opt
            if isempty(ai_aux_s)
                file_ext=[file_ext 'ai'];
            else
                file_ext=[file_ext 'si'];
            end
        end
    else
        file_ext=[file_ext '_vac'];
        curtor=0.0;
        am(:)=0.0;
        ai(:)=0.0;
        ac(:)=0.0;
        ncurr=1;
    end
    file_ext=[file_ext '_' user];
    % Output File
    fid=fopen([filename timestamp file_ext],'w+');
    fprintf(fid,'&INDATA\n');    
    fprintf(fid,'%s\n','!----- Runtime Parameters -----');
    write_namelist_flt(fid,'DELT',delt);
    write_namelist_int(fid,'NITER',niter);
    write_namelist_int(fid,'NSTEP',nstep);
    write_namelist_flt(fid,'TCON0',tcon0);
    write_namelist_vec(fid,'NS_ARRAY',ns_array,'int');
    write_namelist_vec(fid,'FTOL_ARRAY',ftol_array);
%    fprintf(fid,'  LDIAGNO = F\n');
    fprintf(fid,'%s\n','!----- Grid Parameters -----');
    fprintf(fid,'%s\n','  LASYM = F');
    write_namelist_int(fid,'NFP',nfp);
    write_namelist_int(fid,'MPOL',mpol);
    write_namelist_int(fid,'NTOR',ntor);
    write_namelist_int(fid,'NZETA',nzeta);
    write_namelist_flt(fid,'PHIEDGE',phiedge);
    fprintf(fid,'%s\n','!----- Free Boundary Parameters -----');
    fprintf(fid,'  LFREEB = T\n');
    write_namelist_str(fid,'MGRID_FILE',mgrid_file);
    write_namelist_vec(fid,'EXTCUR',extcur,'vert');
    write_namelist_int(fid,'NVACSKIP',nvacskip);
    fprintf(fid,'%s\n','!----- Pressure Parameters -----');
    write_namelist_flt(fid,'GAMMA',gamma);
    write_namelist_flt(fid,'BLOAT',bloat);
    write_namelist_flt(fid,'SPRES_PED',spres_ped);
    if isempty(am_aux_s)
        write_namelist_str(fid,'PMASS_TYPE',pmass_type);
        write_namelist_vec(fid,'AM',am);
    else
        write_namelist_str(fid,'PMASS_TYPE',pmass_type);
        write_namelist_vec(fid,'AM_AUX_S',am_aux_s);
        write_namelist_vec(fid,'AM_AUX_F',am_aux_f);
    end
    fprintf(fid,'%s\n','!----- Current/Iota Parameters -----');
    write_namelist_flt(fid,'CURTOR',signfactor*curtor);
    write_namelist_int(fid,'NCURR',ncurr);
    if isempty(ai_aux_s)
        write_namelist_str(fid,'PIOTA_TYPE',piota_type);
        write_namelist_vec(fid,'AI',ai);
    else
        write_namelist_str(fid,'PIOTA_TYPE',piota_type);
        write_namelist_vec(fid,'AI_AUX_S',ai_aux_s);
        write_namelist_vec(fid,'AI_AUX_F',ai_aux_f);
    end
    if isempty(ac_aux_s)
        write_namelist_str(fid,'PCURR_TYPE',pcurr_type);
        write_namelist_vec(fid,'AC',ac);
    else
        write_namelist_str(fid,'PCURR_TYPE',pcurr_type);
        write_namelist_vec(fid,'AC_AUX_S',ac_aux_s);
        write_namelist_vec(fid,'AC_AUX_F',ac_aux_f);
    end
    fprintf(fid,'%s\n','!----- Axis Parameters -----');
    write_namelist_vec(fid,'RAXIS',raxis);
    write_namelist_vec(fid,'ZAXIS',zaxis);
    fprintf(fid,'%s\n','!----- Boundary Parameters -----');
    fprintf(fid,'  RBC( 0,0) =  %3.2f     ZBS( 0,0) =  %3.2f\n',rbc00,zbs00);
    fprintf(fid,'  RBC(-1,1) =  %3.2f     ZBS(-1,1) =  %3.2f\n',rbcn11,zbsn11);
    fprintf(fid,'  RBC( 0,1) =  %3.2f     ZBS( 0,1) =  %3.2f\n',rbc01,zbs01);
    fprintf(fid,'  RBC( 1,1) =  %3.2f     ZBS( 1,1) =  %3.2f\n',rbc11,zbs11);
    fprintf(fid,'/\n'); % End INDATA namelist
    % Reconstruction Namelist
    if strcmp(outputtype,'recon')
        % Begin Writing Namelist
        fprintf(fid,'&OPTIMUM\n'); 
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          OPTIMIZER RUN CONTROL PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        write_namelist_flt(fid,'EPSFCN',epsfcn);
        write_namelist_int(fid,'NITER_OPT',niter_opt);
        write_namelist_int(fid,'NUM_PROCESSORS',num_processors);
        write_namelist_int(fid,'NUM_LEVMAR_PARAMS',num_levmar_params);
        write_namelist_int(fid,'NOPT_ALG',nopt_alg);
        write_namelist_int(fid,'NOPT_BOUNDARY',nopt_boundary);
        write_namelist_boo(fid,'LRESET_OPT',lreset_opt);
        write_namelist_boo(fid,'LDIAG_OPT',ldiag_opt);
        write_namelist_boo(fid,'LKEEP_MINS',lkeep_mins);
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          PHYSICS MODULES\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        write_namelist_boo(fid,'LPHIEDGE',lphiedge);
        write_namelist_boo(fid,'LIOTA_PROF_OPT',liota_prof_opt);
        write_namelist_boo(fid,'LCUR_PROF_OPT',lcur_prof_opt);
        write_namelist_boo(fid,'LCUR_OPT_EDGE0',lcur_opt_edge0);
        if isempty(ac_aux_s)
            fprintf(fid,'  AC_MASK =  1  0  1  0  1  0  0  0  0  0  0\n');
        end
        if use_mags || ldiamag
            write_namelist_boo(fid,'LDIAGNO_OPT',1);
        else
            write_namelist_boo(fid,'LDIAGNO_OPT',0);
        end
        write_namelist_boo(fid,'LPRES_PROF_OPT',lpres_prof_opt);
        write_namelist_boo(fid,'LPRES_OPT_EDGE0',lpres_opt_edge);
        write_namelist_boo(fid,'LPRES_OPT_EDGEGR0',lpres_opt_edgegr0);
        write_namelist_boo(fid,'LP_PROF_INCL_EDGE',lp_prof_incl_edge);
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          EQUILIBRIUM AND GEOMETRY OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        write_namelist_flt(fid,'TARGET_BETA',target_beta);
        write_namelist_flt(fid,'SIGMA_BETA',abs(sigma_beta));
        write_namelist_flt(fid,'TARGET_CURTOR',signfactor*target_curtor);
        write_namelist_flt(fid,'SIGMA_CURTOR',abs(sigma_curtor));
        write_namelist_flt(fid,'TARGET_EPLASMA',target_eplasma);
        write_namelist_flt(fid,'SIGMA_EPLASMA',abs(sigma_eplasma));
        write_namelist_flt(fid,'TARGET_IOTA_MAX',target_iota_max);
        write_namelist_flt(fid,'SIGMA_IOTA_MAX',abs(sigma_iota_max));
        write_namelist_flt(fid,'TARGET_IOTA_MIN',target_iota_min);
        write_namelist_flt(fid,'SIGMA_IOTA_MIN',abs(sigma_iota_min));
        write_namelist_flt(fid,'SIGMA_PEDGE',abs(sigma_pedge));
        write_namelist_vec(fid,'SIGMA_PGRAD',abs(sigma_pgrad));
        %write_namelist_vec(fid,'TARGET_IOTA_P',target_iota_p);
        %write_namelist_vec(fid,'SIGMA_IOTA_PMAX',sigma_iota_pmax);
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          PRESSURE PROFILE OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        write_namelist_int(fid,'PRES_OPT_NMAX',pres_opt_nmax);
        if lisote, write_namelist_boo(fid,'ISOTE',lisote); end
        write_namelist_int(fid,'NP_PROF',np_prof);
        write_namelist_flt(fid,'FACTOR_P_PROF',factor_p_prof);
        write_namelist_vec(fid,'NE_PROF',ne_prof);
        write_namelist_vec(fid,'TE_PROF',te_prof);
        write_namelist_vec(fid,'SIGMA_P_PROF',abs(sigma_p_prof));
        write_namelist_vec(fid,'R_P_PROF',r_p_prof);
        write_namelist_vec(fid,'Z_P_PROF',z_p_prof);
        write_namelist_vec(fid,'PHI_P_PROF',phi_p_prof);
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          DIAGNO OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        write_namelist_str(fid,'DIAGNO_CONTROL',diagno_control);
        if (ndiagno_seg > 0)
            write_namelist_vec(fid,'TARGET_DIAGNO_SEG',target_diagno_seg);
            write_namelist_vec(fid,'SIGMA_DIAGNO_SEG',abs(sigma_diagno_seg));
        end
        if (ndiagno_flx > 0) && ~ldiamag
            write_namelist_vec(fid,'TARGET_DIAGNO_FLX',signfactor*target_diagno_flx);
            write_namelist_vec(fid,'SIGMA_DIAGNO_FLX',abs(sigma_diagno_flx));
        end
        if (ldiamag)
            write_namelist_flt(fid,'TARGET_DIAGNO_FLX',target_diamag);
            write_namelist_flt(fid,'SIGMA_DIAGNO_FLX',abs(sigma_diamag));
        end
        if luse_mse
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            fprintf(fid,'!          MSE PROFILE OPTIMIZATION PARAMETERS\n');
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            write_namelist_int(fid,'NMSE_CAMS',nmse_cams);
            write_namelist_int(fid,'NMSE_CHORDS',nmse_chords);
            for n=1:nmse_cams
                for k=1:nmse_chords
                    fprintf(fid,'  MSE_R(%3d,%3d) = %20.10e ',n,k,mse_r(n,k));
                    fprintf(fid,'  MSE_PHI(%3d,%3d) = %20.10e ',n,k,mse_phi(n,k));
                    fprintf(fid,'  MSE_Z(%3d,%3d) = %20.10e ',n,k,mse_z(n,k));
                    fprintf(fid,'  MSE_ALPHA(%3d,%3d) = %20.10e ',n,k,mse_alpha(n,k));
                    fprintf(fid,'  MSE_BETA(%3d,%3d) = %20.10e ',n,k,mse_beta(n,k));
                    fprintf(fid,'  MSE_THETA(%3d,%3d) = %20.10e ',n,k,mse_theta(n,k));
                    fprintf(fid,'  MSE_POL(%3d,%3d) = %20.10e ',n,k,mse_pol(n,k));
                    fprintf(fid,'  SIGMA_MSE_POL(%3d,%3d) = %20.10e ',n,k,sigma_mse_pol(n,k));
                    fprintf(fid,'  MSE_VAC(%3d,%3d) = %20.10e \n',n,k,mse_pol_vac(n,k));
                end
            end
        end
        fprintf(fid,'/\n'); % End Optimizer namelist
        fprintf(fid,'&END\n'); % End Namelists
    end
    time=clock;
    time_comment=sprintf('%2d-%2d-%4d %02d:%02d:%02d',time(2),time(3),time(1),time(4),time(5),round(time(6)));
    fprintf(fid,['!-----  Created by LHDrecon_setup (' time_comment ')\n']);    %Add signature
    fclose(fid); % Close File
    % Write batch script
    batchname=strtrim([timestamp file_ext '.batch']);
    fid=fopen(batchname,'w+');
    fprintf(fid,'#Batchname\n');
    fprintf(fid,'#!/bin/tcsh\n');  % NOTE TCSH
    fprintf(fid,'#PBS -N %s\n',strtrim([shotnum '_' timestamp file_ext]));
    fprintf(fid,'# Email\n');
    fprintf(fid,'#PBS -m aeb\n');
    switch user
        case 'SAL'
            fprintf(fid,'#PBS -M lazerson@pppl.gov\n');
        case 'NAP'
            fprintf(fid,'#PBS -M npablant@pppl.gov\n');
    end
    if strcmp(outputtype,'recon')
        fprintf(fid,'# Set Cluster\n');
        fprintf(fid,'#PBS -q kruskal\n');
        fprintf(fid,'# Nodes, Processors, Memory\n');
        if lmse_vac
            fprintf(fid,'#PBS -l nodes=2\n');
        else
            fprintf(fid,'#PBS -l nodes=%-3d\n',nnodes);
        end
    else
        fprintf(fid,'# Nodes, Processors, Memory\n');
        fprintf(fid,'#PBS -l nodes=1\n');
    end 
    fprintf(fid,'# Walltime\n');
    fprintf(fid,'#PBS -l walltime=144:00:00\n');
    fprintf(fid,'# Don''t rerun on fail\n');
    fprintf(fid,'#PBS -r n\n');
    fprintf(fid,'# Export Environment Vars\n');
    fprintf(fid,'#PBS -V\n');
    fprintf(fid,'# Standard Error files\n');
    fprintf(fid,'#PBS -j oe\n');
    fprintf(fid,'# ------------------------------------------------------------\n');
    fprintf(fid,'# Log interesting information\n');
    fprintf(fid,'#\n');
    fprintf(fid,'echo " "\n');
    fprintf(fid,'echo "-------------------"\n');
    fprintf(fid,'echo "This is a $PBS_ENVIRONMENT job"\n');
    fprintf(fid,'echo "This job was submitted to the queue: $PBS_QUEUE"\n');
    fprintf(fid,'echo "The jobs id is: $PBS_JOBID"\n');
    fprintf(fid,'echo "-------------------"\n');
    fprintf(fid,'echo "The master node of this job is: $PBS_O_HOST"\n');
    fprintf(fid,'set NPROCS=`wc -l < $PBS_NODEFILE`\n');
    fprintf(fid,'set NNODES=`uniq $PBS_NODEFILE | wc -l`\n');
    fprintf(fid,'set OMP_NUM_THREADS=$NPROCS\n');
    fprintf(fid,'echo "This job is using $NPROCS CPU(s) on the following $NNODES node(s):"\n');
    fprintf(fid,'echo "-----------------------"\n');
    fprintf(fid,'uniq $PBS_NODEFILE | sort\n');
    fprintf(fid,'echo "-----------------------"\n');
    fprintf(fid,'# ------------------------------------------------------------\n');
    fprintf(fid,'# Setup execution variables\n');
    fprintf(fid,'# --- PBS_O_WORKDIR is the current working directory from\n');
    fprintf(fid,'#     which this job was submitted using ''qsub''\n');
    fprintf(fid,'echo The working directory is $PBS_O_WORKDIR \n');
    fprintf(fid,'#\n');
    fprintf(fid,'# --- the path of your executable\n');
    if strcmp(outputtype,'recon')
        fprintf(fid,'set EXEPATH="/u/slazerso/bin_847/xstellopt"\n');
        %fprintf(fid,'set EXEPATH="/p/lhd/bin/xstellopt"\n');
    else
        fprintf(fid,'set EXEPATH="/p/lhd/bin/xvmec2000"\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'set ARGS="%s"\n',strtrim([shotnum '_' timestamp file_ext]));
    fprintf(fid,'# --- echo the command syntax\n');
    fprintf(fid,'echo "The command syntax for this job is:"\n');
    %fprintf(fid,'echo mpirun --mca btl_openib_want_fork_support 1 -np $NPROCS $EXEPATH $ARGS\n'); % Fork Support in OPENMPI
    fprintf(fid,'echo mpirun -np $NPROCS $EXEPATH $ARGS\n');
    fprintf(fid,'echo " "\n');
    fprintf(fid,'cd $PBS_O_WORKDIR\n');
    fprintf(fid,'echo -n ''Started job at : '' ; date\n');
    fprintf(fid,'time mpirun -np $NPROCS $EXEPATH $ARGS\n');
    fprintf(fid,' \n');
    fprintf(fid,'echo -n ''Ended job at  : '' ; date\n');
    fprintf(fid,'echo " " \n');
    fprintf(fid,'exit\n');
    fclose(fid);
    % Create output string
    string='  ';
    string=[string num2str(timedex1(i),'%3.2f')];
    string=[string '  ' num2str(extcur(1),percf)];
    string=[string '  ' num2str(extcur(2),percf)];
    string=[string '  ' num2str(extcur(3),percf)];
    string=[string '  ' num2str(extcur(4),percf)];
    string=[string '  ' num2str(extcur(5),percf)];
    string=[string '  ' num2str(extcur(6),percf)];
    string=[string '  ' num2str(signfactor*curtor/1000.,percf)];
    string=[string '  ' num2str(target_eplasma/1000.,percf)];
    disp(string);
end
% Write diagno file
if lwrite_diagno && strcmp(outputtype,'recon')
    fid=fopen('diagno.control','w+');
    fprintf(fid,'&diagno_in\n');
    fprintf(fid,'nu = 360\n');
    fprintf(fid,'nv = 72\n');
    if (ldiamag)
        fprintf(fid,'units = 1.\n');  % Diamagnetic loop is in [M]
    else
        fprintf(fid,'units = 1000.\n');
    end
    if (ldiamag)
        fprintf(fid,'flux_diag_file = ''/p/lhd/probes/LHD_diamag.diagno''\n');
        fprintf(fid,'flux_turns = 1\n');
    else
        fprintf(fid,'flux_diag_file = ''/p/lhd/probes/LHD_saddle.diagno''\n');
        fprintf(fid,'flux_turns = -1 1 -1 1 -1 1 -1 1 1 1 -1 1 -1 -1 -1 -1 -1 -1 -1 1 -1 1 -1 1\n');
    end
    fprintf(fid,'lcomp_dia = .false.\n');
    fprintf(fid,'lflux_comb = .false.\n');
    fprintf(fid,'ltrace_progress = .true.\n');
    fprintf(fid,'lwrpl_surf = .false.\n');
    fprintf(fid,'int_type = ''simpson''\n');
    fprintf(fid,'int_step = 2\n');
    fprintf(fid,'/\n');
    fprintf(fid,'&END\n');
    fclose(fid);
end
% Turn all MATLAB warnings back on
warning on

end

