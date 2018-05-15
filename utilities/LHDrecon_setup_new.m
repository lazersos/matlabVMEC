function LHDrecon_setup_new(data,timedex1,timedex2,varargin)
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
%           'linene'        Line Integrated Electron Density
%           'phiedge'       Vary enclosed toroidal flux.
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
tunits_fir=1.0;         % FIR Time units (s)
ne_units=1e16;          % Electron Number Density Units
ne_norm=1e18;           % Normalization used in STELLOPT
ec=1.60217653e-19;      % NIST 2008 pg 34
phi_press=-18.0;        % Phi Location of Pressure Profile (degrees)
Runits_press=0.001;     % units/meter for Radial Pressure Profile. [mm/m]
z_press=0.0;            % Z Location of Pressure Profile
scale_curtor=1000.;     % Toroidal Current is in [kA]
scale_extcur=1000.;     % Coil Currents are in [kA]
scale_eplasma=1000.;    % Scaling Factor for total stored energy [kJ]
scale_fdia=1.0e-4;      % Diamagnetic loop in (10^-4) [Wb]
fir_units=1e19;         % FIR singal in 10^19 [m^-2]
mu0=4*pi*10^-7;         % Mu0
luse_mags=0;             % Switch to use DIAGNO
luse_mse=0;             % Switch to use MSE data
lmse_vac=0;             % Switch to run a single iteration to calculate VMEC MSE vacuum response
luse_pres_spline=0;     % Switch to use pressure splines (v8.47)
luse_curr_spline=0;     % Switch to use current splines (v8.47)
luse_iota_spline=0;     % Switch to use iota splines (v8.47)
lcurtor=0;              % Switch to target total toroidal current
leplasma=0;             % Switch to target total stored energy
ldiamag=0;              % Switch to use diagno for diamagnetic loop
lphiedge=0;
lpres_opt=1;
liota_opt=0;            % Switch to optimize iota instead of I'
lisote=0;
lcoil_opt=0;            % Switch to optimize coil currents
lcurr_opt=0;
llinene=0;
npspline=15;            % Number of datapoints in Pressure spline
njspline=15;            % Number of datapoints in Current spline
nnodes=0;                % Node counter
lefield=0;
lbootstrap=0;           % Boostrap Current
ndiagno_flx = 24;
flux_vals=[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
flux_vals=[0.0 0.1 0.25 0.5 0.75 0.9 1.0];
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
                    luse_mags=1;
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
                    luse_mags=1;
                case {'phiedge'}
                    lphiedge=1;
                case {'current'}
                    lcurr_opt=1;
                case {'coil_opt'}
                    lcoil_opt=1;
                case {'bootstrap'}
                    lbootstrap=1;
                case {'linene'}
                    llinene=1;
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
vmec_input=vmec_namelist_init('INDATA');
vmec_input.lasym=0;
vmec_input.nfp=10;
vmec_input.ntor=6;
vmec_input.ncurr=1;
vmec_input.niter=20000;
vmec_input.nstep=200;
vmec_input.nvacskip=6;
vmec_input.mpol=8;
vmec_input.nzeta=36;
vmec_input.ntheta=2*vmec_input.mpol+6;
vmec_input.ns_array=[9  29  49  99];
vmec_input.ftol_array=[1.0E-6  1.0E-8  1.0E-10  1.0E-12];
vmec_input.delt=1.0;
vmec_input.tcon0=1.0;
vmec_input.phiedge=1.00;
vmec_input.spres_ped=1.0;
vmec_input.bloat=1.0;
vmec_input.pres_scale=2.5;
vmec_input.am=[1 -1 0 0 -1 1];
vmec_input.ac=[1 -1 0 0 -1 1];
vmec_input.am_aux_s=flux_vals;
vmec_input.am_aux_f=1.0E04.*polyval(fliplr(vmec_input.am),vmec_input.am_aux_s);
vmec_input.ac_aux_s=flux_vals;
vmec_input.ac_aux_f=polyval(fliplr(vmec_input.ac),vmec_input.ac_aux_s);
vmec_input.ai_aux_s=flux_vals;
vmec_input.pmass_type='AKIMA_SPLINE';
vmec_input.pcurr_type='AKIMA_SPLINE_IP';
%vmec_input.pcurr_type='TWO_POWER';
%vmec_input.ac=[1 1 1];
vmec_input.lfreeb=1;
vmec_input.mgrid_file='/p/lhd/coils/lhd_ys4/mgrid_lhd_ys4_847_nzeta36.nc';
vmec_input.rbc=zeros(2*vmec_input.ntor+1,vmec_input.mpol);
vmec_input.zbs=zeros(2*vmec_input.ntor+1,vmec_input.mpol);
vmec_input.rbc(vmec_input.ntor+1,1) =  3.7;
vmec_input.rbc(vmec_input.ntor+1,2) =  0.5;
vmec_input.rbc(vmec_input.ntor+2,2) = -0.13;
vmec_input.zbs(vmec_input.ntor+1,2) = -0.5;
vmec_input.zbs(vmec_input.ntor+2,2) = -0.14;
vmec_input.raxis    = [];
vmec_input.zaxis    = [];
vmec_input.raxis(1) = 3.8;
vmec_input.zaxis(1) = 0.0;
vmec_input.datatype='VMEC_input';
% Initialize STELLOPT
stel_input.niter_opt=5000;
stel_input.ftol=1.0E-4;
stel_input.xtol=1.0E-4;
stel_input.gtol=1.0E-8;
stel_input.epsfcn=1.0E-4;
stel_input.factor=100.0;
stel_input.nopt_alg=0;
stel_input.nopt_boundary=0;
stel_input.lreset_opt=1;
stel_input.ldiag_opt=0;
stel_input.lkeep_mins=1;
stel_input.lphiedge=lphiedge;
stel_input.liota_prof_opt=liota_opt;
stel_input.lcur_prof_opt=lcurr_opt;
stel_input.lcur_opt_edge0=1;
stel_input.ldiagno_opt=luse_mags;
stel_input.lpres_prof_opt=lpres_opt;
stel_input.lpres_opt_edge0=1;
stel_input.lpres_opt_edgegr0=0;
stel_input.lp_prof_incl_edge=1;
stel_input.factor_p_prof=1.5;
stel_input.lcurtor_opt=lcurtor;
stel_input.lpscale_opt=1;
stel_input.mboz = 8*vmec_input.mpol;
stel_input.nboz = 6*vmec_input.ntor;

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
% Plot FIR Data
if isfield(data,'fircall_nel') && makeplots
    plot(data.fircall_time.*tunits_fir,data.fircall_nel);
    legend(num2str(data.fircall_r));
    title(['FIR Line Integrated Density (shot:' shotnum ')']);
    xlabel('Time [s]');
    ylabel('Line Ne (x10^19) [m^{-2}]')
    saveas(gca,['fir_' shotnum '.fig']);
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

if makeplots && ~strcmp(outputtype,'recon')
    return
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
    vmec_input.extcur=[mean(ppval(hci,t_temp)) mean(ppval(hcm,t_temp)) mean(ppval(hco,t_temp))...
        -mean(ppval(pov,t_temp)) -mean(ppval(pis,t_temp)) mean(ppval(piv,t_temp))].*scale_extcur;
    stel_input.sigma_extcur=[std(ppval(hci,t_temp)) std(ppval(hcm,t_temp)) std(ppval(hco,t_temp))...
        -std(ppval(pov,t_temp)) -std(ppval(pis,t_temp)) std(ppval(piv,t_temp))].*scale_extcur;
    Ip_spline=pchip(t,data.Ip*scale_curtor);
    % Handle reconstructions
    if strcmp(outputtype,'recon')
        % Thomson Data
        np_prof       = size(data.R_TS,1);
        r_p_prof      = data.R_TS.*Runits_press;
        z_p_prof      = z_press*ones(1,np_prof);
        phi_p_prof    = phi_press.*ones(1,np_prof);
        % Create Thomson Splines
        ne_spline=pchip(data.t_TS*tunits_press,permute(data.n_e.*ne_units,[2 1]));
        te_spline=pchip(data.t_TS*tunits_press,permute(data.T_e,[2 1]));
        % Calculate Thomson Pressure
        ne_prof=mean(permute(ppval(ne_spline,t_temp), [2 1]));
        te_prof=mean(permute(ppval(te_spline,t_temp), [2 1]));
        % Calculate Standard Deviations
        % Note the we define sigma_p_prof by p=n*kB*(T_e+T_i+T_h) so
        % sigma_p=dn*kB*(T_e+T_i+T_h)+n*kB*(dT_e+dT_i+dT_h)
        std_ne=std(permute(ppval(ne_spline,t_temp), [2 1]),1);
        std_te=std(permute(ppval(te_spline,t_temp), [2 1]),1);
        sigma_ne_prof = 0.2*ne_prof;  % 20% error bar on ne
        %sigma_ne_prof = std_ne;
        sigma_te_prof = std_te;
        te_cutoff     = 60.;
        % Filter NE by cutoff
        ne_prof(te_prof < te_cutoff) = 0.0;
        sigma_ne_prof(te_prof < te_cutoff) = bigno;
        % Filter TE
        te_prof(te_prof < te_cutoff) = 0.0;
        sigma_te_prof(te_prof == 0.0)= 1.0;
        % Use sigmas to filter
        te_prof(sigma_te_prof > te_prof) = 0.0;
        sigma_te_prof(te_prof == 0.0)= 1.0;
        ne_prof(sigma_ne_prof > ne_prof) = 0.0;
        sigma_ne_prof(ne_prof == 0.0) = bigno;
        ne_prof = smooth(smooth(ne_prof));
        % Remove errors
        for j=2:length(te_prof)-1
            if (te_prof(j-1) == 0) && (te_prof(j+1) == 0)
                te_prof(j) = 0.0;
                sigma_te_prof(j) = 1.0;
            end
        end
        % Filter
        dex = 0;
        j=1;
        while (dex == 0) && (j < length(te_prof))
            num = sum(te_prof(j:j+9)> 0.0);
            if (num == 9)
                dex = j-1;
            end
            j = j+1;
        end
        te_prof(1:dex) = 0.0;
        sigma_te_prof(1:dex) = 1.0;
        ne_prof(1:dex) = 0.0;
        sigma_ne_prof(1:dex) = bigno;
        dex = 0;
        j=length(te_prof);
        while (dex == 0) && (j > 1)
            num = sum(te_prof(j-9:j)> 0.0);
            if (num == 9)
                dex = j+1;
            end
            j = j-1;
        end
        te_prof(dex:end) = 0.0;
        sigma_te_prof(dex:end) = 1.0;
        ne_prof(dex:end) = 0.0;
        sigma_ne_prof(dex:end) = bigno;
        % Remove bad NE points using TE
        sigma_ne_prof(te_prof == 0.0) = bigno;
        
        % Line Ne (FIR)
        if llinene
            data.fircall_nel(data.fircall_nel<0) = 0.0;
            data.fircall_nel(data.fircall_nel.*fir_units<5e17) = 0.0;
            ne_line_spline = pchip(data.fircall_time.*tunits_fir,data.fircall_nel'.*fir_units);
            ne_line        = mean(permute(ppval(ne_line_spline,t_temp),[2 1]));
            std_ne_line    = std(permute(ppval(ne_line_spline,t_temp),[2 1]));
            %sigma_ne_line  = 10.*std_ne_line;
            sigma_ne_line  = ne_line.*0.2;
            sigma_ne_line(ne_line == 0.0) = bigno;
            
        end
        % Extract Toroidal Current
        vmec_input.curtor=mean(ppval(Ip_spline,t_temp));
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
        if luse_mags
            % Extract Flux values (use delta flux if no vacdata supplied)
            sigma_diagno_flx=sigma_diagno_flx_default;
            for j=1:ndiagno_flx
                flux_spline=pchip(t,data.(['FluxLoop' num2str(j-1)]));
                target_diagno_flx(j)=mean(ppval(flux_spline,t_temp));
                if (sigma_diagno_flx(j) < bigno)
                    sigma_diagno_flx(j)=std(ppval(flux_spline,t_temp),1);
                end
            end
            if (ndiagno_flx > 0)
                percent=sigma_diagno_flx./target_diagno_flx;
                sigma_diagno_flx(percent < 0.20)=0.20.*target_diagno_flx(percent < 0.20);
            end
        end
        if ldiamag
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
            nmse_chords=length(data.R_MSE);
            coil_data=coil_biot_prep(coil_data);
            mse_r=data.R_MSE;
            mse_phi=pi*(data.phi_MSE+108)/180; %Note all angles in degrees (LHD assumes MSE is in PHI angles = CAD - 105)
            mse_z=0.0*data.R_MSE;
            vb     = 1.0;
            mse_a1 = cosd(data.THETA_MSE)';
            mse_a2 = sind(data.alpha_MSE)'.*cosd(data.beta_MSE)';
            mse_a3 = cosd(data.alpha_MSE)'.*cosd(data.beta_MSE)';
            mse_a4 = sind(data.THETA_MSE)'.*sind(data.beta_MSE)';
            mse_a5 = cosd(data.alpha_MSE+data.THETA_MSE)'./vb;
            mse_a6 =-cosd(data.beta_MSE)'./vb;
            mse_a7 = sind(data.beta_MSE)'.*sind(data.alpha_MSE+data.THETA_MSE)'./vb;
            for k=1:nmse_chords
                x=mse_r(k).*cos(mse_phi(k));
                y=mse_r(k).*sin(mse_phi(k));
                z=mse_z(k);
                [bx(k) by(k) bz(k)]=coil_biot(coil_data,x,y,z,vmec_input.extcur);
                br(k)   = bx(k)*cos(mse_phi(k)) + by(k)*sin(mse_phi(k));
                bphi(k) = by(k)*cos(mse_phi(k)) - bx(k)*sin(mse_phi(k));
                %disp(sprintf('%20.10E %20.10E %20.10E %20.10E %20.10E %20.10E %20.10E %20.10E',x,y,z,bx(k),by(k),bz(k),br(k),bphi(k)));
                mse_pol(k) = mean(pchip(data.t_MSE*tunits_mse,pi.*data.Pol_Angle(:,k)./180,t_temp));
                std_mse_pol(k) = std(pchip(data.t_MSE*tunits_mse,pi.*data.Pol_Angle(:,k)./180,t_temp));
                err_mse_pol(k) = mean(pchip(data.t_MSE*tunits_mse,pi.*data.dPol_Angle(:,k)./180,t_temp));
                sigma_mse_pol(k) = max([std_mse_pol(k) err_mse_pol(k)]);
            end
            mse_pol_vac = atan(bz.*mse_a1./(br.*mse_a3+bphi.*mse_a2+bz.*mse_a4));
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
                ac_aux_f=polyval(fliplr(vmec_input.ac),ac_aux_s);
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
        if stel_input.lpres_prof_opt, 
            if lisote
                file_ext=[file_ext 't'];
            else
                file_ext=[file_ext 'p'];
            end
        end
        if lcurr_opt, file_ext=[file_ext 'j']; end
        if liota_opt, file_ext=[file_ext 'i']; end
        if luse_mags, file_ext=[file_ext 'm']; end
        if ldiamag, file_ext=[file_ext 'd']; end
        if leplasma, file_ext=[file_ext 'e']; end
        if lcurtor, file_ext=[file_ext 'I']; end
        if luse_mse, file_ext=[file_ext 'M']; end
        file_ext=[file_ext '_'];
        if luse_pres_spline
            file_ext=[file_ext 'sm'];
        elseif lmse_vac
            file_ext=[file_ext 'vac'];
        else
            file_ext=[file_ext 'am'];
        end
        if luse_curr_spline
            file_ext=[file_ext 'sc'];
        elseif luse_iota_spline
            file_ext=[file_ext 'si'];
        elseif liota_opt
            file_ext=[file_ext 'ai'];
        elseif lcurr_opt
            file_ext=[file_ext 'ac'];
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
    write_vmec_input([filename timestamp file_ext],vmec_input);
    % Reconstruction Namelist
    if strcmp(outputtype,'recon')
        fid=fopen([filename timestamp file_ext],'a');
        fseek(fid,0,'eof');
        % Begin Writing Namelist
        fprintf(fid,'&OPTIMUM\n'); 
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          OPTIMIZER RUN CONTROL PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        write_namelist_int(fid,'NFUNC_MAX',stel_input.niter_opt);
        fprintf(fid,'  EQUIL_TYPE = ''VMEC2000''\n');
        fprintf(fid,'  OPT_TYPE   = ''LMDIF''\n');
        write_namelist_flt(fid,'FTOL',stel_input.ftol);
        write_namelist_flt(fid,'XTOL',stel_input.xtol);
        write_namelist_flt(fid,'GTOL',stel_input.gtol);
        write_namelist_flt(fid,'EPSFCN',stel_input.epsfcn);
        write_namelist_flt(fid,'FACTOR',stel_input.factor);
        write_namelist_boo(fid,'LKEEP_MINS',stel_input.lkeep_mins);
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          PHYSICS MODULES\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        write_namelist_boo(fid,'LPHIEDGE_OPT',stel_input.lphiedge);
        write_namelist_boo(fid,'LCURTOR_OPT',stel_input.lcurtor_opt);
        write_namelist_boo(fid,'LPSCALE_OPT',stel_input.lpscale_opt);
        fprintf(fid,'  LNE_F_OPT =');
        for j = 1: length(flux_vals)
            fprintf(fid,' T ');
        end
        fprintf(fid,'\n');
        fprintf(fid,'  LTE_F_OPT =');
        for j = 1: length(flux_vals)-1
            fprintf(fid,' T ');
        end
        fprintf(fid,' F ');
        fprintf(fid,'\n');
        if (lcurr_opt)
            fprintf(fid,'  LAC_F_OPT =');
            for j = 1: length(flux_vals)
                fprintf(fid,' T ');
            end
        end
        fprintf(fid,'\n');
        if (lefield)
            fprintf(fid,'  LPHI_F_OPT =');
            for j = 1: length(flux_vals)-1
                fprintf(fid,' T ');
            end
            fprintf(fid,' F ');
            fprintf(fid,'\n');
        end
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          PROFILE PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        % Output Ne
        stel_input.ne_aux_s=vmec_input.am_aux_s;
        stel_input.ne_aux_f=mean(ne_prof).*polyval([1],stel_input.ne_aux_s);
        if (llinene)
            stel_input.ne_aux_f=(0.10).*max(ne_line).*polyval([1],stel_input.ne_aux_s);
        end
        dex=find(stel_input.ne_aux_s == max(stel_input.ne_aux_s(4:length(stel_input.ne_aux_s))),1,'first');
        fprintf(fid,'  NE_TYPE = ''two_power''\n');
        write_namelist_vec(fid,'NE_OPT',[stel_input.ne_aux_f(1) 4.0 1.0]);
        write_namelist_vec(fid,'NE_AUX_S',stel_input.ne_aux_s(1:dex));
        write_namelist_vec(fid,'NE_AUX_F',stel_input.ne_aux_f(1:dex)./ne_norm);
        % Output Te
        stel_input.te_aux_s=vmec_input.am_aux_s;
        stel_input.te_aux_f=mean(te_prof).*polyval([1],stel_input.te_aux_s);
        dex=find(stel_input.te_aux_s == max(stel_input.te_aux_s(4:length(stel_input.te_aux_s))),1,'first');
        fprintf(fid,'  TE_TYPE = ''two_power''\n');
        write_namelist_vec(fid,'TE_OPT',[stel_input.te_aux_f(1) 1.0 1.0]);
        write_namelist_vec(fid,'TE_AUX_S',stel_input.te_aux_s(1:dex));
        write_namelist_vec(fid,'TE_AUX_F',stel_input.te_aux_f(1:dex));
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          EQUILIBRIUM AND GEOMETRY OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        if (lcurtor)
            write_namelist_flt(fid,'TARGET_CURTOR',vmec_input.curtor);
            write_namelist_flt(fid,'SIGMA_CURTOR',abs(std(ppval(Ip_spline,t_temp))));
        end
        if (leplasma)
            write_namelist_flt(fid,'TARGET_WP',target_eplasma);
            write_namelist_flt(fid,'SIGMA_WP',sigma_eplasma);
        end
        if (lcoil_opt)
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            fprintf(fid,'!          COIL OPTIMIZATION PARAMETERS\n');
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            for j=1:length(vmec_input.extcur)
                fprintf(fid,'  LEXTCUR(%2.2d) = T  TARGET_EXTCUR(%2.2d) = %20.10f  SIGMA_EXTCUR(%2.2d) = %20.10f \n',...
                    j,j,vmec_input.extcur(j),j,abs(stel_input.sigma_extcur(j)));
            end
        end
        if (luse_mags)
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            fprintf(fid,'!          DIAGNO OPTIMIZATION PARAMETERS\n');
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            itemp=1;
            if (ndiagno_flx > 0)
                for j=1:ndiagno_flx
                    fprintf(fid,'  TARGET_FLUXLOOP(%2.2d) = %20.10f  SIGMA_FLUXLOOP(%2.2d) = %20.10f \n',...
                        j,signfactor*target_diagno_flx(j),j,abs(sigma_diagno_flx(j)));
                    itemp=itemp+1;
                end
            end
            if (ldiamag)
                j=itemp;
                fprintf(fid,'  TARGET_FLUXLOOP(%2.2d) = %20.10f  SIGMA_FLUXLOOP(%2.2d) = %20.10f \n',...
                    j,target_diamag,j,abs(sigma_diamag));
            end
        end
        % Output Ne
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          THOMSON (NE) PROFILE OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        for k=1:length(ne_prof)
            if (ne_prof(k) > 0.0)
                fprintf(fid,['  R_NE(%3.3d) = %20.10E  PHI_NE(%3.3d) = %20.10E' ...
                    '  Z_NE(%3.3d) = %20.10E  TARGET_NE(%3.3d) = %20.10E'...
                    '  SIGMA_NE(%3.3d) = %20.10E\n'],...
                    k,r_p_prof(k),k,pi*phi_p_prof(k)/180,...
                    k,z_p_prof(k),k,ne_prof(k),...
                    k,abs(sigma_ne_prof(k)));
            end
        end
        % Output Te
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          THOMSON (TE) PROFILE OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        for k=1:length(te_prof)
            fprintf(fid,['  R_TE(%3.3d) = %20.10E  PHI_TE(%3.3d) = %20.10E' ...
                '  Z_TE(%3.3d) = %20.10E  TARGET_TE(%3.3d) = %20.10E'...
                '  SIGMA_TE(%3.3d) = %20.10E\n'],...
                k,r_p_prof(k),k,pi*phi_p_prof(k)/180,...
                k,z_p_prof(k),k,te_prof(k),...
                k,abs(sigma_te_prof(k)));
        end
        % Separatrix
        dex = find(te_prof > 0,1,'first');
        k = dex - 1;
        sig_sep = abs(r_p_prof(k) - r_p_prof(dex));
        fprintf(fid,['  R_SEPARATRIX(1,1) = %20.10E  PHI_SEPARATRIX(1,1) = %20.10E' ...
            '  Z_SEPARATRIX(1,1) = %20.10E  TARGET_SEPARATRIX(1,1) = %20.10E'...
            '  SIGMA_SEPARATRIX(1,1) = %20.10E\n'],...
            r_p_prof(k),pi*phi_p_prof(k)/180,...
            z_p_prof(k),te_prof(k),...
            sig_sep);
        dex = find(te_prof > 0,1,'last');
        k = dex + 1;
        sig_sep = abs(r_p_prof(k) - r_p_prof(dex));
        fprintf(fid,['  R_SEPARATRIX(1,2) = %20.10E  PHI_SEPARATRIX(1,2) = %20.10E' ...
            '  Z_SEPARATRIX(1,2) = %20.10E  TARGET_SEPARATRIX(1,2) = %20.10E'...
            '  SIGMA_SEPARATRIX(1,2) = %20.10E\n'],...
            r_p_prof(k),pi*phi_p_prof(k)/180,...
            z_p_prof(k),te_prof(k),...
            sig_sep);
        %LINE NE
        if llinene
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            fprintf(fid,'!          LINE NE OPTIMIZATION PARAMETERS\n');
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            for k = 1 : length(ne_line)
                r0_fir = data.fircall_r(k);  r1_fir = data.fircall_r(k);
                phi0_fir = 2.0*pi*0.5;      phi1_fir = phi0_fir;% 8.5U/L 270 degrees from Noon
                z0_fir     = -1.25;           z1_fir = 1.25;
                fprintf(fid,['  R0_NE_LINE(%3.3d) = %20.10E  PHI0_NE_LINE(%3.3d) = %20.10E  Z0_NE_LINE(%3.3d) = %20.10E'...
                    '  R1_NE_LINE(%3.3d) = %20.10E  PHI1_NE_LINE(%3.3d) = %20.10E  Z1_NE_LINE(%3.3d) = %20.10E'...
                    '  TARGET_NE_LINE(%3.3d) = %20.10E'...
                    '  SIGMA_NE_LINE(%3.3d) = %20.10E\n'],...
                    k,r0_fir,k,phi0_fir,k,z0_fir,...
                    k,r1_fir,k,phi1_fir,k,z1_fir,...
                    k,ne_line(k),k,sigma_ne_line(k));
                
            end
        end
        %MSE
        if luse_mse
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            fprintf(fid,'!          MSE PROFILE OPTIMIZATION PARAMETERS\n');
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            for k=1:length(mse_r)
                fprintf(fid,'  R_MSE(%3d) = %20.10e ',k,mse_r(k));
                fprintf(fid,'  PHI_MSE(%3d) = %20.10e ',k,mse_phi(k));
                fprintf(fid,'  Z_MSE(%3d) = %20.10e ',k,mse_z(k));
                fprintf(fid,'  A1_MSE(%3d) = %20.10e ',k,mse_a1(k));
                fprintf(fid,'  A2_MSE(%3d) = %20.10e ',k,mse_a2(k));
                fprintf(fid,'  A3_MSE(%3d) = %20.10e ',k,mse_a3(k));
                fprintf(fid,'  A4_MSE(%3d) = %20.10e ',k,mse_a4(k));
                fprintf(fid,'  A5_MSE(%3d) = %20.10e ',k,mse_a5(k));
                fprintf(fid,'  A6_MSE(%3d) = %20.10e ',k,mse_a6(k));
                fprintf(fid,'  A7_MSE(%3d) = %20.10e ',k,mse_a7(k));
                fprintf(fid,'  TARGET_MSE(%3d) = %20.10e ',k,mse_pol(k));
                fprintf(fid,'  SIGMA_MSE(%3d) = %20.10e ',k,sigma_mse_pol(k));
                fprintf(fid,'  VAC_MSE(%3d) = %20.10e \n',k,mse_pol_vac(k));
            end
        end
        if lbootstrap
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            fprintf(fid,'!          BOOTSTRAP CURRENT OPTIMIZATION PARAMETERS\n');
            fprintf(fid,'!-----------------------------------------------------------------------\n');
            fprintf(fid,'  MBOZ = %3d \n',stel_input.mboz);
            fprintf(fid,'  NBOZ = %3d \n',stel_input.nboz);
            for k=2:max(vmec_data.ns_array)
                fprintf(fid,'  TARGET_BOOTSTRAP(%3d) = %20.10e ',k,0.0);
                fprintf(fid,'  SIGMA_BOOTSTRAP(%3d) = %20.10e \n',k,vmec_input.curtor*.1);
            end
        end
        fprintf(fid,'/\n'); % End Optimizer namelist
    end
    % Handle the DIAGNO input namelist
    if (luse_mags)
        fprintf(fid,'&DIAGNO_IN\n');
        write_namelist_int(fid,'NU',72);
        write_namelist_int(fid,'NV',36);
        fprintf(fid,'  FLUX_DIAG_FILE = ''/p/lhd/probes/LHD_saddle.diagno''\n');
        fprintf(fid,'  FLUX_TURNS = -1 1 -1 1 -1 1 -1 1 1 1 -1 1 -1 -1 -1 -1 -1 -1 -1 1 -1 1 -1 1\n');
        fprintf(fid,'  INT_TYPE = ''simpson''\n');
        write_namelist_int(fid,'INT_STEP',2);
        write_namelist_boo(fid,'LRPHIZ',0);
        write_namelist_flt(fid,'VC_ADAPT_TOL',1.0E-8);
        write_namelist_flt(fid,'VC_ADAPT_REL',1.0E-2);
        write_namelist_boo(fid,'LVC_FIELD',1);
        fprintf(fid,'/\n');
    end
    fprintf(fid,'&END\n'); % End Namelists
    time=clock;
    time_comment=sprintf('%2d-%2d-%4d %02d:%02d:%02d',time(2),time(3),time(1),time(4),time(5),round(time(6)));
    fprintf(fid,['!-----  Created by LHDrecon_setup_new (' time_comment ')\n']);    %Add signature
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
        fprintf(fid,'# Nodes, Processors, Memory\n');
        if lmse_vac
            fprintf(fid,'#PBS -l nodes=2\n');
        else
            fprintf(fid,'#PBS -l nodes=%-3d\n',nnodes);
            if (nnodes*2 > 99)
                fprintf(fid,'#PBS -l mem=%-3dgb\n',nnodes*2);
            else
                fprintf(fid,'#PBS -l mem=%-2dgb\n',nnodes*2);
            end
        end
    else
        fprintf(fid,'# Nodes, Processors, Memory\n');
        fprintf(fid,'#PBS -l nodes=1\n');
    end 
    fprintf(fid,'# Walltime\n');
    fprintf(fid,'#PBS -l walltime=200:00:00\n');
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
        fprintf(fid,'set EXEPATH="/u/slazerso/bin/xstelloptv2"\n');
        %fprintf(fid,'set EXEPATH="/p/lhd/bin/xstellopt"\n');
    else
        fprintf(fid,'set EXEPATH="/p/lhd/bin_847/xvmec2000"\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'set ARGS="%s"\n',strtrim([shotnum '_' timestamp file_ext]));
    %fprintf(fid,'# --- Set the group to LHD\n');
    %fprintf(fid,'newgrp lhd');
    fprintf(fid,'# --- echo the command syntax\n');
    fprintf(fid,'echo "The command syntax for this job is:"\n');
    fprintf(fid,'echo mpirun --mca btl ^openib -np $NPROCS $EXEPATH input.$ARGS\n');
    fprintf(fid,'echo " "\n');
    fprintf(fid,'cd $PBS_O_WORKDIR\n');
    fprintf(fid,'echo -n ''Started job at : '' ; date\n');
    fprintf(fid,['mkdir ' strtrim([shotnum '_' timestamp file_ext]) '\n']);
    fprintf(fid,['cd ' strtrim([shotnum '_' timestamp file_ext]) '\n']);
    fprintf(fid,['mv ../input.' strtrim([shotnum '_' timestamp file_ext]) ' . \n']);
    fprintf(fid,'time mpirun --mca btl ^openib -np $NPROCS $EXEPATH input.$ARGS >& log.$ARGS\n');
    fprintf(fid,'rm *_opt* \n');
    fprintf(fid,'chgrp -R lhd * \n');
    fprintf(fid,' \n');
    fprintf(fid,'echo -n ''Ended job at  : '' ; date\n');
    fprintf(fid,'echo " " \n');
    fprintf(fid,'exit\n');
    fclose(fid);
    % Create output string
    string='  ';
    string=[string num2str(timedex1(i),'%3.2f')];
    string=[string '  ' num2str(vmec_input.extcur(1),percf)];
    string=[string '  ' num2str(vmec_input.extcur(2),percf)];
    string=[string '  ' num2str(vmec_input.extcur(3),percf)];
    string=[string '  ' num2str(vmec_input.extcur(4),percf)];
    string=[string '  ' num2str(vmec_input.extcur(5),percf)];
    string=[string '  ' num2str(vmec_input.extcur(6),percf)];
    string=[string '  ' num2str(signfactor*vmec_input.curtor/1000.,percf)];
%    string=[string '  ' num2str(target_eplasma/1000.,percf)];
    disp(string);
    end
% Turn all MATLAB warnings back on
warning on

end

