function [ output_args ] = DIIIDrecon_setup(shotnum,timedex1,timedex2,varargin)
%DIIID_recon_setup(shotnum,t1,t2) Produces a STELLOPT input file for DIIID.
%   DIIID_recon_setup(shotnum,t1,t2) Produces a STELLOPT input file where
%   data is averaged from t1 to t2.  The values t1 and t2 may be arrays in
%   which case multiple files are created spanning t1 to t2 where the value
%   of t1 is utilized in the shot name.
%       Optional Arguments:
%       coil_data
%       vessel


warning off
% Defaults
%shotnum=['DIIID' num2str(shotnum,'%d')];
filename=['input.DIIID' num2str(shotnum,'%d') '_'];
outputtype='VMEC';
percf='%+7.4f';
perct='%04.4i';
percd='%6.6d';
coil_file='/Volumes/slazerso/Sims/DIIID/coils/coils.d3d_efbic';
user='XXX';
bigno=1e30;
vesdata=[];
coil_data=[];
mse_vac=[];
makeplots=0;
pcurr_type = 'power_series';
piota_type = 'power_series';
pmass_type = 'power_series';
ndt=100000;              % Number of slices used in averageing between t1 and t2
tunits_press=0.001;     % Thomson Time units (ms)
tunits_mse=1.0;         % MSE Time units (s)
ne_units=1e16;          % Electron Number Density Units
ne_norm=1e18;           % STELLOPT ne normalization
ec=1.60217653e-19;      % NIST 2008 pg 34
phi_press=-18.0;        % Phi Location of Pressure Profile (degrees)
Runits_press=0.001;     % units/meter for Radial Pressure Profile. [mm/m]
z_press=0.0;            % Z Location of Pressure Profile
scale_curtor=1000.;     % Toroidal Current is in [kA]
scale_extcur=1000.;     % Coil Currents are in [kA]
scale_eplasma=1000.;    % Scaling Factor for total stored energy [kJ]
scale_co2=0.5;          % Scaling Factor for CO2 data [m/m^-2] double pass
scale_fdia=1.0e-4;      % Diamagnetic loop in (10^-4) [Wb]
mse_factor=10.;         % Factor to scale sigma for the MSE measurement
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
lcurr_opt=0;
lcoil_opt=0;            % Switch to optimize coil currents
npspline=10;            % Number of datapoints in Pressure spline
njspline=10;            % Number of datapoints in Current spline
lefield=0;              % Controls E-field calcualtion
luse_thomson=0;         % Controls use of Thomson
luse_co2=0;             % Controls use of C02 Interferrometer data
luse_sxr=0;             % Controls use of Soft X-Ray data
luse_cer=0;             % Controls use of Charge Exchange
luse_efit=0;            % Controls use of EFIT01 data
numprocs=40;            % Default number of processors
lhmode=0;               % Use H-mode like profiles
flux_vals=[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];
cer_ana_type='USER';
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
        %elseif strcmp(varargin{i},'shotnum')
        %    i=i+1;
        %    shotnum=varargin{i};
        %    filename=['input.' shotnum '_'];
        elseif strcmp(varargin{i},'user') || strcmp(varargin{i},'USER')
            i=i+1;
            user=varargin{i};
        else
            switch lower(varargin{i})
                case {'plots'}
                    makeplots=1;
                case {'magdiags'}
                    luse_mags=1;
                case {'msedata'}
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
                case {'diamag'}
                    ldiamag=1;
                case {'phiedge'}
                    lphiedge=1;
                case {'current'}
                    lcurr_opt=1;
                case {'efield'}
                    lefield=1;
                case {'coil_opt'}
                    lcoil_opt=1;
                case {'thomson'}
                    luse_thomson=1;
                case {'co2'}
                    luse_co2=1;
                case {'sxr'}
                    luse_sxr=1;
                case {'cer'}
                    luse_cer=1;
                case {'efit'}
                    luse_efit=1;
                case {'hmode'}
                    lhmode=1;
            end
        end
        i=i+1;
    end
end

% Initialize VMEC
vmec_input=vmec_namelist_init('INDATA');
vmec_input=rmfield(vmec_input,'bcrit');
vmec_input.lasym=1;
vmec_input.nfp=1;
vmec_input.ntor=3;
vmec_input.ncurr=1;
vmec_input.niter=20000;
vmec_input.nstep=200;
vmec_input.nvacskip=6;
vmec_input.mpol=12;
vmec_input.nzeta=48;
vmec_input.ntheta=2*vmec_input.mpol+6;
vmec_input.ns_array=[16  32  64  128];
vmec_input.niter_array=[2500  5000 10000 20000];
vmec_input.ftol_array=[1.0E-06  1.0E-08  1.0E-10  1.0E-12];
vmec_input.delt=1.0;
vmec_input.tcon0=1.0;
vmec_input.phiedge=-3.0;
vmec_input.spres_ped=1.0;
vmec_input.bloat=1.0;
vmec_input.pres_scale=1.0;
vmec_input.am=[1  1.5  1];
vmec_input.ac=[1 -1  1];
vmec_input.am_aux_s=flux_vals;
vmec_input.am_aux_f=1.0E04.*polyval(fliplr(vmec_input.am),vmec_input.am_aux_s);
vmec_input.ac_aux_s=flux_vals;
vmec_input.ac_aux_f=polyval(fliplr(vmec_input.ac),vmec_input.ac_aux_s);
vmec_input.ai_aux_s=flux_vals;
vmec_input.pmass_type='AKIMA_SPLINE';
if (lhmode)
    vmec_input.pcurr_type='AKIMA_SPLINE_IP';
else
    vmec_input.pcurr_type='TWO_POWER';
end
vmec_input.lfreeb=1;
vmec_input.mgrid_file='/u/slazerso/Sims/DIIID/coils/mgrid_d3d_efic_kp48.nc';
vmec_input.mgrid_file='/p/pies/Lazerson/Sims/MGRIDS/mgrid_d3d_efbic_kp48.nc'; % B Coils
vmec_input.rbc=zeros(2*vmec_input.ntor+1,vmec_input.mpol);
vmec_input.zbs=zeros(2*vmec_input.ntor+1,vmec_input.mpol);
vmec_input.rbc(vmec_input.ntor+1,1) =  1.7;
vmec_input.rbc(vmec_input.ntor+1,2) =  0.5;
vmec_input.zbs(vmec_input.ntor+1,2) =  1.0;
vmec_input.rbc(vmec_input.ntor+1,3) =  0.05;
vmec_input.zbs(vmec_input.ntor+1,3) = -0.05;
vmec_input.raxis    = [];
vmec_input.zaxis    = [];
vmec_input.raxis(1) = 1.7;
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


% Open a connection to the D3D server
import MdsPlus.*;   % So we can call MdsPlus java calls
server = 'atlas.gat.com';  % Server you wish to connect to
server_port = 8000;        % Server port
mds_server=MdsPlus(server,server_port);  % Open the connection

% Open Tree
tree     = 'd3d';              % Tree you wish to open
tree_mse = 'MSE';              % Tree you wish to open
tree_end = tree;
%shot_num=142603;           % Shot Number
mds_server.OpenTree(tree,shotnum);  % Open the Tree

% Now get the coil_currents
disp('  - Requesting Coil Currents');
disp('       - F Coils');
% F-Coil Currents
for i=1:9
    coil_str=['PTDATA("F' num2str(i,'%1d') 'A")'];
    coil_time_str=['DIM_OF(PTDATA("F' num2str(i,'%1d') 'A"))'];
    coil.(['F' num2str(i,'%1d') 'A'])=mds_server.Value(coil_str).Double;
    coil.(['F' num2str(i,'%1d') 'A_time'])=mds_server.Value(coil_time_str).Double;
end
for i=1:9
    coil_str=['PTDATA("F' num2str(i,'%1d') 'B")'];
    coil_time_str=['DIM_OF(PTDATA("F' num2str(i,'%1d') 'B"))'];
    coil.(['F' num2str(i,'%1d') 'B'])=mds_server.Value(coil_str).Double;
    coil.(['F' num2str(i,'%1d') 'B_time'])=mds_server.Value(coil_time_str).Double;
end
% E Coil
disp('       - E Coils');
coil_str='PTDATA("ECOILA")';
coil_time_str='DIM_OF(PTDATA("ECOILA"))';
coil.('ECA')=mds_server.Value(coil_str).Double;
coil.('ECA_time')=mds_server.Value(coil_time_str).Double;
coil_str='PTDATA("ECOILB")';
coil_time_str='DIM_OF(PTDATA("ECOILB"))';
coil.('ECB')=mds_server.Value(coil_str).Double;
coil.('ECB_time')=mds_server.Value(coil_time_str).Double;
% B Coil (I_line = 144*BCOIL)
disp('       - B Coil');
coil_str='PTDATA("BCOIL")';
coil_time_str='DIM_OF(PTDATA("BCOIL"))';
%coil.('B')=144.*mds_server.Value(coil_str).Double; % for EFIC coil
coil.('B')=-6.*mds_server.Value(coil_str).Double; % for EFBIC coils
coil.('B_time')=mds_server.Value(coil_time_str).Double;
% C Coils
disp('       - C Coils');
coil_str='PTDATA("C139")';
coil_time_str='DIM_OF(PTDATA("C139"))';
coil.('C139')=mds_server.Value(coil_str).Double;
coil.('C139_time')=mds_server.Value(coil_time_str).Double;
coil_str='PTDATA("C79")';
coil_time_str='DIM_OF(PTDATA("C79"))';
coil.('C79')=mds_server.Value(coil_str).Double;
coil.('C79_time')=mds_server.Value(coil_time_str).Double;
coil_str='PTDATA("C199")';
coil_time_str='DIM_OF(PTDATA("C199"))';
coil.('C199')=mds_server.Value(coil_str).Double;
coil.('C199_time')=mds_server.Value(coil_time_str).Double;
% I Coils
disp('       - I Coils');
icoil_names={'IU330' 'IU270' 'IU210' 'IU150' 'IU90' 'IU30'...
    'IL330' 'IL270' 'IL210' 'IL150' 'IL90' 'IL30'};
for i =1:length(icoil_names);
    coil_str=['PTDATA("' upper(icoil_names{i}) '")'];
    coil_time_str=['DIM_OF(PTDATA("' upper(icoil_names{i}) '"))'];
    coil_name=upper(icoil_names{i});
    coil.(coil_name)=mds_server.Value(coil_str).Double;
    coil.([coil_name '_time'])=mds_server.Value(coil_time_str).Double;
end
% Get Net toroidal Current
disp('  - Requesting Toroidal Current');
ip_str='PTDATA("IP")';
ip_time_str='DIM_OF(PTDATA("IP"))';
ip=mds_server.Value(ip_str).Double;
ip_time=mds_server.Value(ip_time_str).Double;
% Get Magneticts
if (luse_mags)
    disp('  - Requesting Magnetics');
    % Use this list for 3D probes
    bprobe_names={'MPI11M067', 'MPI2A067', 'MPI1B157', 'MPI2B067', 'MPI6NB157',...
        'MPI8A322', 'MPI89A322', 'MPI9A322', 'MPI79FA322', 'MPI79NA322',...
        'MPI7FA322', 'MPI7NA322', 'MPI67A322', 'MPI6FA322', 'MPI6NA322',...
        'MPI66M322', 'MPI6NB322', 'MPI6FB322', 'MPI67B322', 'MPI7NB322',...
        'MPI7FB322', 'MPI79B322', 'MPI9B322', 'MPI89B322', 'MPI8B322',...
        'MPI5B322', 'MPI4B322', 'MPI3B322', 'MPI2B322', 'MPI1B322',...
        'MPI11M322', 'MPI1A322', 'MPI2A322', 'MPI3A322', 'MPI4A322',...
        'MPI5A322', 'MPI1U157', 'MPI2U157', 'MPI3U157', 'MPI4U157',...
        'MPI5U157', 'MPI6U157', 'MPI7U157', 'MPI1L180', 'MPI2L180',...
        'MPI3L180', 'MPI1L020', 'MPI2L020', 'MPI1L050', 'MPI1L110',...
        'MPI1L230', 'MPI1L320', 'MPI66M067', 'MPI66M097', 'MPI66M127',...
        'MPI66M137', 'MPI66M157', 'MPI66M247', 'MPI66M277', 'MPI66M307',...
        'MPI66M312', 'MPI66M340', 'MPI67A1', 'MPI67A2', 'MPI67A3', 'MPI67A4',...
        'MPI67A5', 'MPI67A6', 'MPI67B1', 'MPI67B2', 'MPI67B3', 'MPI67B4', 'MPI67B5',...
        'MPI67B6'};
    %Sorted_names
    floop_names={'psf1a','psi12a','psf2a','psi23a','psf3a',...
        'psi34a','psf4a','psi45a','psf5a','psi58a',...
        'psf8a','psf9a','psi9a','psf7fa','psi7a',...
        'psf7na','psf6fa','psi6a','psf6na','psf6nb',...
        'psi6b','psf6fb','psf7nb','psi7b','psf7fb',...
        'psi3l','psi9b','psi2l','psf9b','psi89nb',...
        'psi1l','psi89fb','psf8b','psi58b','psf5b',...
        'psi45b','psf4b','psi34b','psf3b','psi23b',...
        'psf2b','psi12b','psf1b','psi11m',...
        'DSL1U180', 'DSL2U180', 'DSL3U180', 'DSL4U157',...
        'DSL5U157', 'DSL6U157', 'DSL12A067', 'DSL34A067', 'DSL59A067',...
        'DSL79A067', 'DSL67A067', 'DSL66M052', 'DSL67B067', 'DSL79B067',...
        'DSL59B067', 'DSL34B067', 'DSL12B067', 'DSL12A157', 'DSL34A157',...
        'DSL59A157', 'DSL79A157', 'DSL67A157', 'DSL66M152', 'DSL67B157',...
        'DSL79B157', 'DSL59B157', 'DSL34B157', 'DSL12B157', 'SL67FA345',...
        'SL67NA345', 'SL66A132', 'SL66B132', 'SL66A312', 'SL66B312',...
        'SL67NB015', 'SL67FB015', 'ESL019', 'ESL079', 'ESL139', 'ESL199',...
        'ESL259', 'ESL319', 'MISL1', 'MISL2', 'MISL3', 'MISL4', 'MISL5',...
        'MISL6', 'UISL1', 'UISL2', 'UISL3', 'UISL4', 'UISL5', 'UISL6', 'LISL1',...
        'LISL2', 'LISL3', 'LISL4', 'LISL5', 'LISL6'};
    % Flux loop scaling [Wb/rad] => [Wb]
    floop_scales=[6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     1.4849254993E-01...
        1.4716084060E-01     1.8889468502E-01     3.8361191333E-02     2.0396870728E-01     2.1800034317E-01...
        5.2359877560E-01     5.2359877560E-01     5.1836278784E-01     5.1661745859E-01     5.1836278784E-01...
        4.3284165449E-01     5.2185344635E-01     5.2185344635E-01     5.2185344635E-01     5.2359877560E-01...
        5.2359877560E-01     5.2359877560E-01     5.2359877560E-01     5.2359877560E-01     5.2185344635E-01...
        5.2534410485E-01     5.0789081233E-01     5.2359877560E-01     5.2359877560E-01     5.2359877560E-01...
        5.2359877560E-01     5.2359877560E-01     3.6651914292E-01     3.6651914292E-01     2.2863813201E-01...
        2.2863813201E-01     2.2340214426E-01     2.2340214426E-01     2.6179938780E-01     2.6179938780E-01...
        3.0959826973E+00     3.0959826973E+00     3.0959826973E+00     3.0959826973E+00     3.0959826973E+00...
        3.0959826973E+00     2.0468601827E+00     1.7006251269E+00     2.3863063159E+00     2.0230989534E+00...
        1.7006251269E+00     2.3625450865E+00     1.6780217639E+00     1.3342319391E+00     1.9399568685E+00...
        1.6370944038E+00     1.3342319391E+00     1.8990295084E+00     1.6289089317E+00     1.3342319391E+00...
        1.9399568685E+00     1.6370944038E+00     1.3342319391E+00     1.9481423405E+00];
    % Reordered 3D Floops
    floop_names={'psf1a','psi12a','psf2a','psi23a','psf3a',...
        'psi34a','psf4a','psi45a','psf5a','psi58a',...
        'psf8a','psf9a','psi9a','psf7fa','psi7a',...
        'psf7na','psf6fa','psi6a','psf6na','psf6nb',...
        'psi6b','psf6fb','psf7nb','psi7b','psf7fb',...
        'psi3l','psi9b','psi2l','psf9b','psi89nb',...
        'psi1l','psi89fb','psf8b','psi58b','psf5b',...
        'psi45b','psf4b','psi34b','psf3b','psi23b',...
        'psf2b','psi12b','psf1b','psi11m',...
        'DSL6U157', 'DSL5U157', 'DSL1U180', 'DSL2U180', 'DSL3U180' ...
        'DSL4U157', 'DSL12A067', 'DSL12A157', 'DSL34A067', 'DSL34A157' ...
        'DSL59A067', 'DSL59A157', 'DSL79A067', 'DSL79A157', 'DSL67A067' ...
        'DSL67A157', 'DSL66M052', 'DSL66M152', 'DSL67B067', 'DSL67B157' ...
        'DSL79B067', 'DSL79B157', 'DSL59B067', 'DSL59B157', 'DSL34B067' ...
        'DSL34B157', 'DSL12B067', 'DSL12B157', 'SL67FA345', 'SL67NA345' ...
        'SL66A132', 'SL66B132', 'SL66A312', 'SL66B312', 'SL67NB015' ...
        'SL67FB015', 'ESL019', 'ESL079', 'ESL139', 'ESL199' ...
        'ESL259', 'ESL319', 'UISL1', 'UISL2', 'UISL3' ...
        'UISL4', 'UISL5', 'UISL6', 'MISL1', 'MISL2' ...
        'MISL3', 'MISL4', 'MISL5', 'MISL6', 'LISL1' ...
        'LISL2', 'LISL3', 'LISL4', 'LISL5', 'LISL6'};
    floop_scales=[6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00...
        6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     6.2831853072E+00     2.1800034317E-01...
        2.0396870728E-01     1.4849254993E-01     1.4716084060E-01     1.8889468502E-01...
        3.8361191333E-02     5.2359877560E-01     5.2359877560E-01     5.2359877560E-01     5.2359877560E-01...
        5.1836278784E-01     5.2359877560E-01     5.1661745859E-01     5.2185344635E-01     5.1836278784E-01...
        5.2534410485E-01     4.3284165449E-01     5.0789081233E-01     5.2185344635E-01     5.2359877560E-01...
        5.2185344635E-01     5.2359877560E-01     5.2185344635E-01     5.2359877560E-01     5.2359877560E-01...
        5.2359877560E-01     5.2359877560E-01     5.2359877560E-01     3.6651914292E-01     3.6651914292E-01...
        2.2863813201E-01     2.2863813201E-01     2.2340214426E-01     2.2340214426E-01     2.6179938780E-01...
        2.6179938780E-01     3.0959826973E+00     3.0959826973E+00     3.0959826973E+00     3.0959826973E+00...
        3.0959826973E+00     3.0959826973E+00     1.6780217639E+00     1.3342319391E+00     1.9399568685E+00...
        1.6370944038E+00     1.3342319391E+00     1.8990295084E+00     2.0468601827E+00     1.7006251269E+00...
        2.3863063159E+00     2.0230989534E+00     1.7006251269E+00     2.3625450865E+00     1.6289089317E+00...
        1.3342319391E+00     1.9399568685E+00     1.6370944038E+00     1.3342319391E+00     1.9481423405E+00]; 
    disp('       - B Probes');
    bprobe=[];
    for j=1:length(bprobe_names)
        bprobe_str=['PTDATA("' upper(bprobe_names{j}) '")'];
        bprobe_time_str=['DIM_OF(PTDATA("' upper(bprobe_names{j}) '"))'];
        try
            bprobe.(upper(bprobe_names{j}))=mds_server.Value(bprobe_str).Double;
            bprobe.([upper(bprobe_names{j}) '_TIME'])=mds_server.Value(bprobe_time_str).Double;
        catch temp_err
            bprobe.(upper(bprobe_names{j}))=[0.0 0.0];
            bprobe.([upper(bprobe_names{j}) '_TIME'])=[-1E30 1.0E30];
        end
    end
    disp('       - Flux Loops');
    floop=[];
    floop.PSF1A=0;  % Do this because PSF1A is the offset loop and read first.
    for j=1:length(floop_names)
        floop_str=['PTDATA("' upper(floop_names{j}) '")'];
        floop_time_str=['DIM_OF(PTDATA("' upper(floop_names{j}) '"))'];
        try
            if strfind(floop_names{j},'ps')
                floop.(upper(floop_names{j}))=floop_scales(j).*mds_server.Value(floop_str).Double+floop.PSF1A; % Note we need to multiply by 2*pi to get proper signal because signal in [V*s/rad]=[Wb/rad]
                %floop.(upper(floop_names{j}))=floop_scales(j).*mds_server.Value(floop_str).Double; % Note we need to multiply by 2*pi to get proper signal because signal in [V*s/rad]=[Wb/rad]
                %floop.(upper(floop_names{j}))=mds_server.Value(floop_str).Double;
            else
                floop.(upper(floop_names{j}))=floop_scales(j).*mds_server.Value(floop_str).Double;
            end
            floop.([upper(floop_names{j}) '_TIME'])=mds_server.Value(floop_time_str).Double;
        catch temp_err
            floop.(upper(floop_names{j}))=[0.0 0.0];
            floop.([upper(floop_names{j}) '_TIME'])=[-1E30 1.0E30];
        end
    end
end
% Get C02 Data (correct toroidal angle)
if (luse_co2)
    disp('  - Requesting CO2 Inteferrometer Data');
    co2r0_den_str='\ELECTRONS::NELINE_R0';
    co2v1_den_str='\ELECTRONS::NELINE_V1';
    co2v2_den_str='\ELECTRONS::NELINE_V2';
    co2v3_den_str='\ELECTRONS::NELINE_V3';
    co2_t_str=['DIM_OF(' co2r0_den_str ')'];
    co2_t = mds_server.Value(co2_t_str).Double;
    co2r0_den = mds_server.Value(co2r0_den_str).Double.*scale_co2;
    co2v1_den = mds_server.Value(co2v1_den_str).Double.*scale_co2;
    co2v2_den = mds_server.Value(co2v2_den_str).Double.*scale_co2;
    co2v3_den = mds_server.Value(co2v3_den_str).Double.*scale_co2;
    co2r0_r1=2.365;  co2r0_r2=1.0;
    co2r0_phi1=pi*(360-225)/180.; co2r0_phi2=co2r0_phi1;
    co2r0_z1=0.0;    co2r0_z2=0.0;
    co2v1_r1=1.48;  co2v1_r2=1.48;
    co2v1_phi1=pi*(360-240)/180.; co2v1_phi2=co2v1_phi1;
    co2v1_z1=-1.25;    co2v1_z2=1.20;
    co2v2_r1=1.94;  co2v2_r2=1.94;
    co2v2_phi1=pi*(360-240)/180.; co2v2_phi2=co2v2_phi1;
    co2v2_z1=-1.10;    co2v2_z2=1.07;
    co2v3_r1=2.10;  co2v3_r2=2.10;
    co2v3_phi1=pi*(360-240)/180.; co2v3_phi2=co2v3_phi1;
    co2v3_z1=-1.00;    co2v3_z2=1.03;
end
% Get SXR Data
if (luse_sxr)
    disp('  - Requesting Soft X-Ray');
    sxr45r1s_r0(1:12) = [ 2.2599013657E+00      2.2646056886E+00      2.2667260830E+00      2.2693982716E+00      2.2722341143E+00      2.2757290892E+00...
        2.2790790486E+00      2.2830942621E+00      2.2899305029E+00      2.2985566921E+00      2.3089854102E+00      2.3225209339E+00 ];
    sxr45r1s_r1(1:12) = [ 1.0160000000E+00      1.1236166179E+00      1.1902688708E+00      1.2868898374E+00      1.3798436892E+00      1.5276495135E+00...
        1.6114021748E+00      1.7011925070E+00      1.8620884358E+00      2.0045283105E+00      2.1312485066E+00      2.2620084266E+00 ];
    sxr45r1s_phi0(1:12) = [ 5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00 ];
    sxr45r1s_phi1(1:12) = [ 5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00      5.4977871438E+00 ];
    sxr45r1s_z0(1:12) = [ 6.7270180144E-01      6.6124255330E-01      6.5607749013E-01      6.4956831282E-01      6.4266049079E-01      6.3414709051E-01 ...
        6.2598693293E-01      6.1620628472E-01      5.9955390327E-01      5.7854139115E-01      5.5313810340E-01      5.2016695594E-01 ];
    sxr45r1s_z1(1:12) = [ -1.1800947158E+00     -1.3329731862E+00     -1.3630000000E+00     -1.3630000000E+00     -1.3630000000E+00     -1.2500000000E+00...
        -1.2500000000E+00     -1.2500000000E+00     -1.1300523690E+00     -1.0477810620E+00     -9.7458922466E-01     -6.6535834922E-01 ];
    sxr165r1s_r0(1:12) = [ 2.2599013657E+00      2.2646056886E+00      2.2667260830E+00      2.2693982716E+00      2.2722341143E+00      2.2757290892E+00...
        2.2790790486E+00      2.2830942621E+00      2.2899305029E+00      2.2985566921E+00      2.3089854102E+00      2.3225209339E+00 ];
    sxr165r1s_r1(1:12) = [ 1.0160000000E+00      1.1236166179E+00      1.1902688708E+00      1.2868898374E+00      1.3798436892E+00      1.5276495135E+00...
        1.6114021748E+00      1.7011925070E+00      1.8620884358E+00      2.0045283105E+00      2.1312485066E+00      2.2620084266E+00 ];
    sxr165r1s_phi0(1:12) = [ 3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00 ];
    sxr165r1s_phi1(1:12) = [ 3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00      3.4033920414E+00 ];
    sxr165r1s_z0(1:12) = [ 6.7270180144E-01      6.6124255330E-01      6.5607749013E-01      6.4956831282E-01      6.4266049079E-01      6.3414709051E-01...
        6.2598693293E-01      6.1620628472E-01      5.9955390327E-01      5.7854139115E-01      5.5313810340E-01      5.2016695594E-01 ];
    sxr165r1s_z1(1:12) = [ -1.1800947158E+00     -1.3329731862E+00     -1.3630000000E+00     -1.3630000000E+00     -1.3630000000E+00     -1.2500000000E+00...
        -1.2500000000E+00     -1.2500000000E+00     -1.1300523690E+00     -1.0477810620E+00     -9.7458922466E-01     -6.6535834922E-01 ];
    sxr195r1s_r0(1:12) = [ 2.2599013657E+00      2.2646056886E+00      2.2667260830E+00      2.2693982716E+00      2.2722341143E+00      2.2757290892E+00...
        2.2790790486E+00      2.2830942621E+00      2.2899305029E+00      2.2985566921E+00      2.3089854102E+00      2.3225209339E+00 ];
    sxr195r1s_r1(1:12) = [ 1.0160000000E+00      1.1236166179E+00      1.1902688708E+00      1.2868898374E+00      1.3798436892E+00      1.5276495135E+00...
        1.6114021748E+00      1.7011925070E+00      1.8620884358E+00      2.0045283105E+00      2.1312485066E+00      2.2620084266E+00 ];
    sxr195r1s_phi0(1:12) = [ 2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00 ];
    sxr195r1s_phi1(1:12) = [ 2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00      2.8797932658E+00 ];
    sxr195r1s_z0(1:12) = [ 6.7270180144E-01      6.6124255330E-01      6.5607749013E-01      6.4956831282E-01      6.4266049079E-01      6.3414709051E-01...
        6.2598693293E-01      6.1620628472E-01      5.9955390327E-01      5.7854139115E-01      5.5313810340E-01      5.2016695594E-01 ];
    sxr195r1s_z1(1:12) = [ -1.1800947158E+00     -1.3329731862E+00     -1.3630000000E+00     -1.3630000000E+00     -1.3630000000E+00     -1.2500000000E+00...
        -1.2500000000E+00     -1.2500000000E+00     -1.1300523690E+00     -1.0477810620E+00     -9.7458922466E-01     -6.6535834922E-01 ];
    sxr90p1s_r0(1:32) = [ 2.3227210819E+00      2.3140393069E+00      2.3064047447E+00      2.2996386294E+00      2.2936007243E+00      2.2881794564E+00...
        2.2832849437E+00      2.2788439658E+00      2.2747962731E+00      2.2710918362E+00      2.2655474117E+00      2.2623125714E+00      2.2576489615E+00...
        2.2531979800E+00      2.2489453968E+00      2.2448778802E+00      2.2409726658E+00      2.2385830805E+00      2.2263422216E+00      2.2238414565E+00...
        2.2212775970E+00      2.2186482249E+00      2.2159507966E+00      2.2131826352E+00      2.2103409215E+00      2.2059338972E+00      2.1997713110E+00...
        2.1960383759E+00      2.1925374774E+00      2.1887291083E+00      2.1845709022E+00      2.1800123344E+00 ];
    sxr90p1s_r1(1:32) = [ 2.0337005161E+00      1.9213063376E+00      1.7954676094E+00      1.6537811058E+00      1.5312302153E+00      1.3494741251E+00...
        1.2085868316E+00      1.0966202366E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00...
        1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00...
        1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00...
        1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.1654980847E+00 ];
    sxr90p1s_phi0(1:32) = [ 4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00 ];
    sxr90p1s_phi1(1:32) = [ 4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00 ];
    sxr90p1s_z0(1:32) = [ 5.1967941596E-01      5.4082732940E-01      5.5942433983E-01      5.7590590273E-01      5.9061362035E-01      6.0381927280E-01...
        6.1574180372E-01      6.2655957051E-01      6.3641933475E-01      6.4544296320E-01      6.5894861248E-01      6.6682835175E-01      6.7818842719E-01...
        6.8903056164E-01      6.9938941816E-01      7.0929642756E-01      7.1877575142E-01      7.2457611238E-01      7.5428896631E-01      7.6035919968E-01...
        7.6658258505E-01      7.7296499259E-01      7.7951259633E-01      7.8623189411E-01      7.9312972908E-01      8.0382712126E-01      8.1878587747E-01...
        8.2784701920E-01      8.3634492675E-01      8.4558917308E-01      8.5568259635E-01      8.6674783777E-01 ];
    sxr90p1s_z1(1:32) = [ -1.0309315984E+00     -1.0958489257E+00     -1.1685316394E+00     -1.2500000000E+00     -1.2500000000E+00     -1.3630000000E+00...
        -1.3630000000E+00     -1.3053856432E+00     -1.2041616619E+00     -1.0075602341E+00     -7.3922208967E-01     -5.9535681673E-01     -4.0245364412E-01...
        -2.3264002118E-01     -8.2004109645E-02      5.2529267837E-02      1.7341115071E-01      2.4385108053E-01      2.0119519058E-01      2.7472392487E-01...
        3.4725729120E-01      4.1881536583E-01      4.8941768871E-01      5.5908328125E-01      6.2783066345E-01      7.2926941507E-01      8.6147178098E-01...
        9.3657419267E-01      1.0038764007E+00      1.0738752099E+00      1.1467360014E+00      1.1769297763E+00 ];
    sxr90m1s_r0(1:32) = [ 2.0649408855E+00      1.9646968086E+00      1.8589732099E+00      1.7489474262E+00      1.6349912968E+00      2.290000E+00...
        1.3276821317E+00      1.2739313654E+00      1.2427848095E+00      1.2017784639E+00      1.1477147742E+00      1.0825906890E+00      1.0128671111E+00...
        1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00...
        1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00...
        1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0160000000E+00      1.0813369619E+00 ];
    sxr90m1s_r1(1:32) = [ 2.3257488361E+00      2.3165187178E+00      2.3084376303E+00      2.3013034867E+00      2.2949590136E+00      1.4139506063E+00...
        2.2866580764E+00      2.2841668355E+00      2.2817967061E+00      2.2795390656E+00      2.2773860908E+00      2.2753306679E+00      2.2733663134E+00...
        2.2696876263E+00      2.2654807693E+00      2.2622434600E+00      2.2590976247E+00      2.2545420741E+00      2.2501717421E+00      2.2459755524E+00...
        2.2419432953E+00      2.2380655451E+00      2.2165679882E+00      2.2095945088E+00      2.2053764856E+00      2.2009750524E+00      2.1963779806E+00...
        2.1915719298E+00      2.1935738799E+00      2.1884648009E+00      2.1826271110E+00      2.1740447085E+00 ];
    sxr90m1s_phi0(1:32) = [ 4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00 ];
    sxr90m1s_phi1(1:32) = [ 4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00      4.7123889804E+00 ];
    sxr90m1s_z0(1:32) = [ 1.0406567973E+00      1.0536709406E+00      1.0673964605E+00      1.0770000000E+00      1.0825424784E+00      -0.60006E+00...
        1.3480000000E+00      1.3163638813E+00      1.2453080159E+00      1.1965838479E+00      1.1716931586E+00      1.1638322514E+00      1.1610979999E+00...
        9.6655067791E-01      7.6557250545E-01      6.2159396936E-01      4.8973038238E-01      3.1150402050E-01      1.5325132970E-01      1.1792695964E-02...
        -1.1541059981E-01     -2.3041011034E-01     -4.0572484238E-01     -5.9468370535E-01     -6.9963401089E-01     -8.0246671742E-01     -9.0324527641E-01...
        -1.0020306293E+00     -9.6165765761E-01     -1.0626099725E+00     -1.1700492102E+00     -1.2897676982E+00 ];
    sxr90m1s_z1(1:32) = [ -5.1217152151E-01     -5.3435419262E-01     -5.5377540706E-01     -5.7092083854E-01     -5.8616846104E-01      1.3100000000E+00...
        -6.0611803858E-01     -6.1210521831E-01     -6.1780133177E-01     -6.2322710169E-01     -6.2840132904E-01     -6.3334111080E-01     -6.3806202873E-01...
        -6.4690298872E-01     -6.5701329526E-01     -6.6479349522E-01     -6.7235385677E-01     -6.8330217590E-01     -6.9380536056E-01     -7.0389003051E-01...
        -7.1358072230E-01     -7.2290008925E-01     -7.7456499944E-01     -7.9132430815E-01     -8.0146145026E-01     -8.1203938021E-01     -8.2308748693E-01...
        -8.3463783136E-01     -8.2982656030E-01     -8.4210516988E-01     -8.5613484433E-01     -8.7676086509E-01 ];
    sxr_time_str='\SPECTROSCOPY::SX90RP1S:TIMEBASE';
    sxr90p1s_time=mds_server.Value(sxr_time_str).Double;
    sxr_time_str='\SPECTROSCOPY::SX90RM1S:TIMEBASE';
    sxr90m1s_time=mds_server.Value(sxr_time_str).Double;
    sxr_time_str='\SPECTROSCOPY::SX45R1S:TIMEBASE';
    sxr45r1s_time=mds_server.Value(sxr_time_str).Double;
    sxr_time_str='\SPECTROSCOPY::SX165R1S:TIMEBASE';
    sxr165r1s_time=mds_server.Value(sxr_time_str).Double;
    sxr_time_str='\SPECTROSCOPY::SX195R1S:TIMEBASE';
    sxr195r1s_time=mds_server.Value(sxr_time_str).Double;
    for i=1:32
        sxr_str = ['\SPECTROSCOPY::SX90RP1S:SX90RP1S' num2str(i,'%2.2d')];
        try
            sxr90p1s(i,:)=mds_server.Value(sxr_str).Double;
        catch
            sxr90p1s(i,:)=0.0;
        end
        sxr_str = ['\SPECTROSCOPY::SX90RM1S:SX90RM1S' num2str(i,'%2.2d')];
        try
            sxr90m1s(i,:)=mds_server.Value(sxr_str).Double;
        catch
            sxr90m1s(i,:)=0.0;
        end
        sxr_str = ['\SPECTROSCOPY::SX45R1S:SX45R1S' num2str(i,'%2.2d')];
        try
            sxr45r1s(i,:)=mds_server.Value(sxr_str).Double;
        catch
            sxr45r1s(i,:)=0.0;
        end
        sxr_str = ['\SPECTROSCOPY::SX165R1S:SX165R1S' num2str(i,'%2.2d')];
        try
            sxr165r1s(i,:)=mds_server.Value(sxr_str).Double;
        catch
            sxr165r1s(i,:)=0.0;
        end
        sxr_str = ['\SPECTROSCOPY::SX195R1S:SX195R1S' num2str(i,'%2.2d')];
        try
            sxr195r1s(i,:)=mds_server.Value(sxr_str).Double;
        catch
            sxr195r1s(i,:)=0.0;
        end
    end
end
% Get Thomson
if (luse_thomson)
    disp('  - Requesting Core Thomson');
    thomson_r_core_str='\D3D::TOP.ELECTRONS.TS.BLESSED.CORE:R';
    thomson_z_core_str='\D3D::TOP.ELECTRONS.TS.BLESSED.CORE:Z';
    thomson_phi_core_str='\D3D::TOP.ELECTRONS.TS.BLESSED.CORE:PHI';
    thomson_time_core_str='\D3D::TOP.ELECTRONS.TS.BLESSED.CORE:TIME';
    thomson_ne_core_str='\D3D::TOP.ELECTRONS.TS.BLESSED.CORE:DENSITY';
    thomson_signe_core_str='\D3D::TOP.ELECTRONS.TS.BLESSED.CORE:DENSITY_E';
    thomson_te_core_str='\D3D::TOP.ELECTRONS.TS.BLESSED.CORE:TEMP';
    thomson_sigte_core_str='\D3D::TOP.ELECTRONS.TS.BLESSED.CORE:TEMP_E';
    thom_r_core=mds_server.Value(thomson_r_core_str).Double;
    thom_z_core=mds_server.Value(thomson_z_core_str).Double;
    thom_phi_core=mds_server.Value(thomson_phi_core_str).Double;
    thom_time_core=mds_server.Value(thomson_time_core_str).Double;
    thom_ne_core=mds_server.Value(thomson_ne_core_str).Double;
    thom_signe_core=mds_server.Value(thomson_signe_core_str).Double;
    thom_te_core=mds_server.Value(thomson_te_core_str).Double;
    thom_sigte_core=mds_server.Value(thomson_sigte_core_str).Double;
    thom_ne_core = reshape(thom_ne_core,[length(thom_time_core) length(thom_ne_core)/length(thom_time_core)]);
    thom_te_core = reshape(thom_te_core,[length(thom_time_core) length(thom_te_core)/length(thom_time_core)]);
    thom_signe_core = reshape(thom_signe_core,[length(thom_time_core) length(thom_signe_core)/length(thom_time_core)]);
    thom_sigte_core = reshape(thom_sigte_core,[length(thom_time_core) length(thom_sigte_core)/length(thom_time_core)]);
    % Divertor Thomson
    disp('  - Requesting Divertor Thomson');
    thomson_r_div_str='\D3D::TOP.ELECTRONS.TS.BLESSED.DIVERTOR:R';
    thomson_z_div_str='\D3D::TOP.ELECTRONS.TS.BLESSED.DIVERTOR:Z';
    thomson_phi_div_str='\D3D::TOP.ELECTRONS.TS.BLESSED.DIVERTOR:PHI';
    thomson_time_div_str='\D3D::TOP.ELECTRONS.TS.BLESSED.DIVERTOR:TIME';
    thomson_ne_div_str='\D3D::TOP.ELECTRONS.TS.BLESSED.DIVERTOR:DENSITY';
    thomson_signe_div_str='\D3D::TOP.ELECTRONS.TS.BLESSED.DIVERTOR:DENSITY_E';
    thomson_te_div_str='\D3D::TOP.ELECTRONS.TS.BLESSED.DIVERTOR:TEMP';
    thomson_sigte_div_str='\D3D::TOP.ELECTRONS.TS.BLESSED.DIVERTOR:TEMP_E';
    thom_r_div=mds_server.Value(thomson_r_div_str).Double;
    thom_z_div=mds_server.Value(thomson_z_div_str).Double;
    thom_phi_div=mds_server.Value(thomson_phi_div_str).Double;
    thom_time_div=mds_server.Value(thomson_time_div_str).Double;
    thom_ne_div=mds_server.Value(thomson_ne_div_str).Double;
    thom_signe_div=mds_server.Value(thomson_signe_div_str).Double;
    thom_te_div=mds_server.Value(thomson_te_div_str).Double;
    thom_sigte_div=mds_server.Value(thomson_sigte_div_str).Double;
    thom_ne_div = reshape(thom_ne_div,[length(thom_time_div) length(thom_ne_div)/length(thom_time_div)]);
    thom_te_div = reshape(thom_te_div,[length(thom_time_div) length(thom_te_div)/length(thom_time_div)]);
    thom_signe_div = reshape(thom_signe_div,[length(thom_time_div) length(thom_signe_div)/length(thom_time_div)]);
    thom_sigte_div = reshape(thom_sigte_div,[length(thom_time_div) length(thom_sigte_div)/length(thom_time_div)]);
    % Tangential Thomson
    disp('  - Requesting Tangential Thomson');
    thomson_r_tan_str='\D3D::TOP.ELECTRONS.TS.BLESSED.TANGENTIAL:R';
    thomson_z_tan_str='\D3D::TOP.ELECTRONS.TS.BLESSED.TANGENTIAL:Z';
    thomson_phi_tan_str='\D3D::TOP.ELECTRONS.TS.BLESSED.TANGENTIAL:PHI';
    thomson_time_tan_str='\D3D::TOP.ELECTRONS.TS.BLESSED.TANGENTIAL:TIME';
    thomson_ne_tan_str='\D3D::TOP.ELECTRONS.TS.BLESSED.TANGENTIAL:DENSITY';
    thomson_signe_tan_str='\D3D::TOP.ELECTRONS.TS.BLESSED.TANGENTIAL:DENSITY_E';
    thomson_te_tan_str='\D3D::TOP.ELECTRONS.TS.BLESSED.TANGENTIAL:TEMP';
    thomson_sigte_tan_str='\D3D::TOP.ELECTRONS.TS.BLESSED.TANGENTIAL:TEMP_E';
    thom_r_tan=mds_server.Value(thomson_r_tan_str).Double;
    thom_z_tan=mds_server.Value(thomson_z_tan_str).Double;
    thom_phi_tan=mds_server.Value(thomson_phi_tan_str).Double;
    thom_time_tan=mds_server.Value(thomson_time_tan_str).Double;
    thom_ne_tan=mds_server.Value(thomson_ne_tan_str).Double;
    thom_signe_tan=mds_server.Value(thomson_signe_tan_str).Double;
    thom_te_tan=mds_server.Value(thomson_te_tan_str).Double;
    thom_sigte_tan=mds_server.Value(thomson_sigte_tan_str).Double;
    thom_ne_tan = reshape(thom_ne_tan,[length(thom_time_tan) length(thom_ne_tan)/length(thom_time_tan)]);
    thom_te_tan = reshape(thom_te_tan,[length(thom_time_tan) length(thom_te_tan)/length(thom_time_tan)]);
    thom_signe_tan = reshape(thom_signe_tan,[length(thom_time_tan) length(thom_signe_tan)/length(thom_time_tan)]);
    thom_sigte_tan = reshape(thom_sigte_tan,[length(thom_time_tan) length(thom_sigte_tan)/length(thom_time_tan)]);
end
if luse_cer
    disp('  - Requesting Charge Exchange');
    [cer.r cer.phi cer.z cer.ti cer.v cer.sig_ti cer.t]=get_d3d_cer(shotnum,mds_server,cer_ana_type);
end
if luse_efit
    disp('  - Requesting EFIT');
    efit_nbdry_str='\D3D::TOP.MHD.EFIT.EFIT01.RESULTS.GEQDSK:NBDRY';
    efit_bdry_str='\D3D::TOP.MHD.EFIT.EFIT01.RESULTS.GEQDSK:BDRY';
    efit_time_str='\D3D::TOP.MHD.EFIT.EFIT01.RESULTS.GEQDSK:GTIME';
    efit_rmag_str='\D3D::TOP.MHD.EFIT.EFIT01.RESULTS.GEQDSK:RMAXIS';
    efit_zmag_str='\D3D::TOP.MHD.EFIT.EFIT01.RESULTS.GEQDSK:ZMAXIS';
    efit_nbdry=mds_server.Value(efit_nbdry_str).Double;
    efit_bdry=mds_server.Value(efit_bdry_str).Double;
    efit_time=mds_server.Value(efit_time_str).Double;
    efit_rmag=mds_server.Value(efit_rmag_str).Double;
    efit_zmag=mds_server.Value(efit_zmag_str).Double;
    i1=1;
    nbdry=round(max(efit_nbdry));
    efit_r=zeros(nbdry,length(efit_nbdry));
    efit_z=zeros(nbdry,length(efit_nbdry));
    for i=1:length(efit_nbdry)
        i2=i1+2*nbdry-1;
        efit_r(1:nbdry,i)=efit_bdry(i1:2:i2);
        efit_z(1:nbdry,i)=efit_bdry(i1+1:2:i2);
        i1= i2+1;
    end
end

% Output file
nfiles=length(timedex1);
for i=1:nfiles
    % Create File name
    tfile_name=[filename num2str(timedex1(i)*1000.,'%4.4d') '_' user];
    % Create Timewindow
    timewindow=timedex1(i):(timedex2(i)-timedex1(i))/ndt:timedex2(i);
    % Get the coil currents
    vmec_input.extcur   = zeros(1,36);
    vmec_input.extcur(1)=sum(pchip(coil.F1A_time./1000,coil.F1A,timewindow))./(ndt+1);
    vmec_input.extcur(3)=sum(pchip(coil.F2A_time./1000,coil.F2A,timewindow))./(ndt+1);
    vmec_input.extcur(5)=sum(pchip(coil.F3A_time./1000,coil.F3A,timewindow))./(ndt+1);
    vmec_input.extcur(7)=sum(pchip(coil.F4A_time./1000,coil.F4A,timewindow))./(ndt+1);
    vmec_input.extcur(9)=sum(pchip(coil.F5A_time./1000,coil.F5A,timewindow))./(ndt+1);
    vmec_input.extcur(11)=sum(pchip(coil.F6A_time./1000,coil.F6A,timewindow))./(ndt+1);
    vmec_input.extcur(13)=sum(pchip(coil.F7A_time./1000,coil.F7A,timewindow))./(ndt+1);
    vmec_input.extcur(15)=sum(pchip(coil.F8A_time./1000,coil.F8A,timewindow))./(ndt+1);
    vmec_input.extcur(17)=sum(pchip(coil.F9A_time./1000,coil.F9A,timewindow))./(ndt+1);
    vmec_input.extcur(2)=sum(pchip(coil.F1B_time./1000,coil.F1B,timewindow))./(ndt+1);
    vmec_input.extcur(4)=sum(pchip(coil.F2B_time./1000,coil.F2B,timewindow))./(ndt+1);
    vmec_input.extcur(6)=sum(pchip(coil.F3B_time./1000,coil.F3B,timewindow))./(ndt+1);
    vmec_input.extcur(8)=sum(pchip(coil.F4B_time./1000,coil.F4B,timewindow))./(ndt+1);
    vmec_input.extcur(10)=sum(pchip(coil.F5B_time./1000,coil.F5B,timewindow))./(ndt+1);
    vmec_input.extcur(12)=sum(pchip(coil.F6B_time./1000,coil.F6B,timewindow))./(ndt+1);
    vmec_input.extcur(14)=sum(pchip(coil.F7B_time./1000,coil.F7B,timewindow))./(ndt+1);
    vmec_input.extcur(16)=sum(pchip(coil.F8B_time./1000,coil.F8B,timewindow))./(ndt+1);
    vmec_input.extcur(18)=sum(pchip(coil.F9B_time./1000,coil.F9B,timewindow))./(ndt+1);
    vmec_input.extcur(19)=sum(pchip(coil.B_time./1000,coil.B,timewindow))./(ndt+1);
    vmec_input.extcur(20)=sum(pchip(coil.ECA_time./1000,coil.ECA,timewindow))./(ndt+1);
    vmec_input.extcur(21)=sum(pchip(coil.ECB_time./1000,coil.ECB,timewindow))./(ndt+1);
    vmec_input.extcur(22)=sum(pchip(coil.IU330_time./1000,coil.IU330,timewindow))./(ndt+1);
    vmec_input.extcur(23)=sum(pchip(coil.IU270_time./1000,coil.IU270,timewindow))./(ndt+1);
    vmec_input.extcur(24)=sum(pchip(coil.IU210_time./1000,coil.IU210,timewindow))./(ndt+1);
    vmec_input.extcur(25)=sum(pchip(coil.IU150_time./1000,coil.IU150,timewindow))./(ndt+1);
    vmec_input.extcur(26)=sum(pchip(coil.IU90_time./1000,coil.IU90,timewindow))./(ndt+1);
    vmec_input.extcur(27)=sum(pchip(coil.IU30_time./1000,coil.IU30,timewindow))./(ndt+1);
    vmec_input.extcur(28)=sum(pchip(coil.IL330_time./1000,coil.IL330,timewindow))./(ndt+1);
    vmec_input.extcur(29)=sum(pchip(coil.IL270_time./1000,coil.IL270,timewindow))./(ndt+1);
    vmec_input.extcur(30)=sum(pchip(coil.IL210_time./1000,coil.IL210,timewindow))./(ndt+1);
    vmec_input.extcur(31)=sum(pchip(coil.IL150_time./1000,coil.IL150,timewindow))./(ndt+1);
    vmec_input.extcur(32)=sum(pchip(coil.IL90_time./1000,coil.IL90,timewindow))./(ndt+1);
    vmec_input.extcur(33)=sum(pchip(coil.IL30_time./1000,coil.IL30,timewindow))./(ndt+1);
    vmec_input.extcur(34)=sum(pchip(coil.C79_time./1000,coil.C79,timewindow))./(ndt+1);
    vmec_input.extcur(35)=sum(pchip(coil.C139_time./1000,coil.C139,timewindow))./(ndt+1);
    vmec_input.extcur(36)=sum(pchip(coil.C199_time./1000,coil.C199,timewindow))./(ndt+1);
    % Get the VMEC boundary guess
    if luse_efit
        for j=1:size(efit_r,1)
            rbdry(j)=sum(pchip(efit_time./1000,efit_r(j,:),timewindow))./(ndt+1);
            zbdry(j)=sum(pchip(efit_time./1000,efit_z(j,:),timewindow))./(ndt+1);
            rax=sum(pchip(efit_time./1000,efit_rmag,timewindow))./(ndt+1);
            zax=sum(pchip(efit_time./1000,efit_zmag,timewindow))./(ndt+1);
        end
        zbdry = zbdry(rbdry > 0.0); % Get rid of zero dimesions (z first then r since we test by r)
        rbdry = rbdry(rbdry > 0.0); % Get rid of zero dimesions
        xbdry=rbdry-rax;
        ybdry=zbdry-zax;
        rhobdry=sqrt(xbdry.*xbdry+ybdry.*ybdry);
        thetabdry=atan2(ybdry,xbdry);
        thetabdry(thetabdry<0.0) = thetabdry(thetabdry<0.0) + 2*pi;
        [thetabdry, idex] = sort(thetabdry);
        rbdry2=rbdry(idex);
        zbdry2=zbdry(idex);
        ntheta=360;
        theta2=0:2*pi/ntheta:2*pi;
        r2=pchip(thetabdry,rbdry2,theta2);
        z2=pchip(thetabdry,zbdry2,theta2);
        L=sum(r2 > 0.0);
        NFFT=2^nextpow2(L);
        rfft=fft(r2,NFFT)/L;
        zfft=fft(z2,NFFT)/L;
        for j=1:vmec_input.mpol+1
            vmec_input.rbc(vmec_input.ntor+1,j)=real(rfft(j));
            vmec_input.zbs(vmec_input.ntor+1,j)=imag(zfft(j));
            vmec_input.rbs(vmec_input.ntor+1,j)=imag(rfft(j));
            vmec_input.zbc(vmec_input.ntor+1,j)=real(zfft(j));
        end
        vmec_input.raxis_cc(1)=rax;
        vmec_input.zaxis_cc(1)=zax;
    end
    % Get Itor
    vmec_input.curtor=sum(pchip(ip_time./1000,ip,timewindow))./(ndt+1);
    % Write input namelist
    write_vmec_input(tfile_name,vmec_input);
    % Open file and go to end
    fid=fopen(tfile_name,'a');
    fseek(fid,0,'eof');
    % Write the recon nameslist
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
    %write_namelist_flt(fid,'EPSFCN',stel_input.epsfcn);
    %write_namelist_int(fid,'NITER_OPT',stel_input.niter_opt);
    %write_namelist_int(fid,'NOPT_ALG',stel_input.nopt_alg);
    %write_namelist_int(fid,'NOPT_BOUNDARY',stel_input.nopt_boundary);
    %write_namelist_boo(fid,'LRESET_OPT',stel_input.lreset_opt);
    %write_namelist_boo(fid,'LDIAG_OPT',stel_input.ldiag_opt);
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    fprintf(fid,'!          OPTIMIZED QUANTITIES\n');
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    write_namelist_boo(fid,'LPHIEDGE_OPT',stel_input.lphiedge);
    write_namelist_flt(fid,'DPHIEDGE_OPT',1.0);
    write_namelist_boo(fid,'LCURTOR_OPT',stel_input.lcurtor_opt);
    write_namelist_flt(fid,'DCURTOR_OPT',0.5);
    write_namelist_boo(fid,'LPSCALE_OPT',stel_input.lpscale_opt);
    write_namelist_flt(fid,'DPSCALE_OPT',1.0);
    fprintf(fid,'  LNE_F_OPT =');
    for j = 1: length(flux_vals)-1
        fprintf(fid,'    T ');
    end
    fprintf(fid,'    F ');
    fprintf(fid,'\n');
    fprintf(fid,'  DNE_F_OPT =');
    for j = 1: length(flux_vals)
        fprintf(fid,' 0.25 ');
    end
    fprintf(fid,'\n');
    fprintf(fid,'  LTE_F_OPT =');
    for j = 1: length(flux_vals)-1
        fprintf(fid,'    T ');
    end
    fprintf(fid,'    F ');
    fprintf(fid,'\n');
    fprintf(fid,'  DTE_F_OPT =');
    for j = 1: length(flux_vals)-1
        fprintf(fid,' 0.25 ');
    end
    fprintf(fid,'\n');
    fprintf(fid,'  LTI_F_OPT =');
    for j = 1: length(flux_vals)-1
        fprintf(fid,'    T ');
    end
    fprintf(fid,'    F ');
    fprintf(fid,'\n');
    fprintf(fid,'  DTI_F_OPT =');
    for j = 1: length(flux_vals)
        fprintf(fid,' 0.25 ');
    end
    fprintf(fid,'\n');
    if (lefield)
        fprintf(fid,'  LPHI_F_OPT =');
        for j = 1: length(flux_vals)
            fprintf(fid,'    T ');
        end
        fprintf(fid,'\n');
    end
    if (lhmode)
        fprintf(fid,'  LBEAMJ_F_OPT = F T F\n');
        fprintf(fid,'  LBOOTJ_F_OPT = T T F\n');
    else
        fprintf(fid,'  LAC_OPT = F T F\n');
    end
    %write_namelist_boo(fid,'LPHIEDGE',stel_input.lphiedge);
    %write_namelist_boo(fid,'LIOTA_PROF_OPT',stel_input.liota_prof_opt);
    %write_namelist_boo(fid,'LCUR_PROF_OPT',stel_input.lcur_prof_opt);
    %write_namelist_boo(fid,'LCUR_OPT_EDGE0',stel_input.lcur_opt_edge0);
    %write_namelist_boo(fid,'LDIAGNO_OPT',stel_input.ldiagno_opt);
    %write_namelist_boo(fid,'LPRES_PROF_OPT',stel_input.lpres_prof_opt);
    %write_namelist_boo(fid,'LPRES_OPT_EDGE0',stel_input.lpres_opt_edge0);
    %write_namelist_boo(fid,'LPRES_OPT_EDGEGR0',stel_input.lpres_opt_edgegr0);
    %write_namelist_boo(fid,'LP_PROF_INCL_EDGE',stel_input.lp_prof_incl_edge);
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    fprintf(fid,'!          PROFILE DATA\n');
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    if (luse_thomson)
        % Filter
        cutoff=80;
        te_core=sum(pchip(thom_time_core./1000,thom_te_core',timewindow),2)./(ndt+1);
        ne_core=sum(pchip(thom_time_core./1000,thom_ne_core',timewindow),2)./(ndt+1);
        sigma_te_core=std(pchip(thom_time_core./1000,thom_te_core',timewindow),0,2);
        sigma_ne_core=std(pchip(thom_time_core./1000,thom_ne_core',timewindow),0,2);
        te_core(te_core < cutoff) = 0.0;
        ne_core(te_core < cutoff) = 0.0;
        sigma_ne_core(ne_core == 0.0) = 1.0E30;
        sigma_te_core(te_core == 0.0) = 1.0E-4;
        cutoff=80;
        te_div=sum(pchip(thom_time_div./1000,thom_te_div',timewindow),2)./(ndt+1);
        ne_div=sum(pchip(thom_time_div./1000,thom_ne_div',timewindow),2)./(ndt+1);
        sigma_te_div=std(pchip(thom_time_div./1000,thom_te_div',timewindow),0,2);
        sigma_ne_div=std(pchip(thom_time_div./1000,thom_ne_div',timewindow),0,2);
        te_div(te_div < cutoff) = 0.0;
        ne_div(te_div < cutoff) = 0.0;
        sigma_ne_div(ne_div == 0.0) = 1.0E30;
        sigma_te_div(te_div == 0.0) = 1.0E-4;
        cutoff=80;
        te_tan=sum(pchip(thom_time_tan./1000,thom_te_tan',timewindow),2)./(ndt+1);
        ne_tan=sum(pchip(thom_time_tan./1000,thom_ne_tan',timewindow),2)./(ndt+1);
        sigma_te_tan=std(pchip(thom_time_tan./1000,thom_te_tan',timewindow),0,2);
        sigma_ne_tan=std(pchip(thom_time_tan./1000,thom_ne_tan',timewindow),0,2);
        te_tan(te_tan < cutoff) = 0.0;
        ne_tan(te_tan < cutoff) = 0.0;
        sigma_ne_tan(ne_tan == 0.0) = 1.0E30;
        sigma_te_tan(te_tan == 0.0) = 1.0E30;
        % Output Profiles
        stel_input.ne_aux_s=vmec_input.am_aux_s;
        stel_input.ne_aux_f=0.9*max(ne_core).*polyval([-1 0 0 0 0 1],stel_input.ne_aux_s);
        dex=find(stel_input.ne_aux_s == max(stel_input.ne_aux_s(4:length(stel_input.ne_aux_s))),1,'first');
        stel_input.ne_aux_f(dex) = 0.0;
        fprintf(fid,'  NE_TYPE   = ''akima_spline''\n');
        write_namelist_vec(fid,'NE_AUX_S',stel_input.ne_aux_s(1:dex));
        write_namelist_vec(fid,'NE_AUX_F',stel_input.ne_aux_f(1:dex)./ne_norm);
        stel_input.te_aux_s=vmec_input.am_aux_s;
        stel_input.te_aux_f=0.95*max(te_core).*polyval([-1 0 1],stel_input.te_aux_s);
        dex=find(stel_input.te_aux_s == max(stel_input.te_aux_s(4:length(stel_input.te_aux_s))),1,'first');
        %stel_input.te_aux_f(1:dex-1) = 0.75.*max(te_core);
        stel_input.te_aux_f(dex) = 0.0;
        fprintf(fid,'  TE_TYPE   = ''akima_spline''\n');
        write_namelist_vec(fid,'TE_AUX_S',stel_input.te_aux_s(1:dex));
        write_namelist_vec(fid,'TE_AUX_F',stel_input.te_aux_f(1:dex));
    end
    if luse_cer
        stel_input.ti_aux_s=vmec_input.am_aux_s;
        for j=1:length(cer.r)
            if isempty(cer.ti{j})
                ti_temp(j) = 0.0;
            else
                ti_temp(j) = sum(pchip(cer.t{j}./1000,cer.ti{j},timewindow))./(ndt+1);
            end
            if (ti_temp(j) < 0.0) || (ti_temp(j) > 1.0E4)
                ti_temp(j) = 0.0;
            end
        end
        stel_input.ti_aux_f=max(ti_temp).*polyval([-1 1],stel_input.ti_aux_s);
        dex=find(stel_input.ti_aux_s == max(stel_input.ti_aux_s(4:length(stel_input.ti_aux_s))),1,'first');
        fprintf(fid,'  TI_TYPE   = ''akima_spline''\n');
        write_namelist_vec(fid,'TI_AUX_S',stel_input.ti_aux_s(1:dex));
        write_namelist_vec(fid,'TI_AUX_F',stel_input.ti_aux_f(1:dex));
    end
    if (lhmode)
        fprintf(fid,'  BEAMJ_TYPE   = ''two_power''\n');
        write_namelist_vec(fid,'BEAMJ_AUX_S',[0.0 0.5 1.0]);
        write_namelist_vec(fid,'BEAMJ_AUX_F',[1.0 1.0 1.0]);
        fprintf(fid,'  BOOTJ_TYPE   = ''bump''\n');
        write_namelist_vec(fid,'BOOTJ_AUX_S',[0.0 0.5 1.0]);
        write_namelist_vec(fid,'BOOTJ_AUX_F',[0.95 0.5 0.0]);
    end
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    fprintf(fid,'!          EQUILIBRIUM AND GEOMETRY OPTIMIZATION PARAMETERS\n');
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    if (lcurtor)
        write_namelist_flt(fid,'TARGET_CURTOR',vmec_input.curtor);
        write_namelist_flt(fid,'SIGMA_CURTOR',abs(std(pchip(ip_time./1000,ip,timewindow))));
    end
    %write_namelist_flt(fid,'FACTOR_P_PROF',stel_input.factor_p_prof);
    if (lcoil_opt)
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          COIL OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        stel_input.sigma_extcur(1)=std(pchip(coil.F1A_time./1000,coil.F1A,timewindow));
        stel_input.sigma_extcur(3)=std(pchip(coil.F2A_time./1000,coil.F2A,timewindow));
        stel_input.sigma_extcur(5)=std(pchip(coil.F3A_time./1000,coil.F3A,timewindow));
        stel_input.sigma_extcur(7)=std(pchip(coil.F4A_time./1000,coil.F4A,timewindow));
        stel_input.sigma_extcur(9)=std(pchip(coil.F5A_time./1000,coil.F5A,timewindow));
        stel_input.sigma_extcur(11)=std(pchip(coil.F6A_time./1000,coil.F6A,timewindow));
        stel_input.sigma_extcur(13)=std(pchip(coil.F7A_time./1000,coil.F7A,timewindow));
        stel_input.sigma_extcur(15)=std(pchip(coil.F8A_time./1000,coil.F8A,timewindow));
        stel_input.sigma_extcur(17)=std(pchip(coil.F9A_time./1000,coil.F9A,timewindow));
        stel_input.sigma_extcur(2)=std(pchip(coil.F1B_time./1000,coil.F1B,timewindow));
        stel_input.sigma_extcur(4)=std(pchip(coil.F2B_time./1000,coil.F2B,timewindow));
        stel_input.sigma_extcur(6)=std(pchip(coil.F3B_time./1000,coil.F3B,timewindow));
        stel_input.sigma_extcur(8)=std(pchip(coil.F4B_time./1000,coil.F4B,timewindow));
        stel_input.sigma_extcur(10)=std(pchip(coil.F5B_time./1000,coil.F5B,timewindow));
        stel_input.sigma_extcur(12)=std(pchip(coil.F6B_time./1000,coil.F6B,timewindow));
        stel_input.sigma_extcur(14)=std(pchip(coil.F7B_time./1000,coil.F7B,timewindow));
        stel_input.sigma_extcur(16)=std(pchip(coil.F8B_time./1000,coil.F8B,timewindow));
        stel_input.sigma_extcur(18)=std(pchip(coil.F9B_time./1000,coil.F9B,timewindow));
        stel_input.sigma_extcur(19)=std(pchip(coil.B_time./1000,coil.B,timewindow));
        stel_input.sigma_extcur(20)=std(pchip(coil.ECA_time./1000,coil.ECA,timewindow));
        stel_input.sigma_extcur(21)=std(pchip(coil.ECB_time./1000,coil.ECB,timewindow));
        stel_input.sigma_extcur(22)=std(pchip(coil.IU330_time./1000,coil.IU330,timewindow));
        stel_input.sigma_extcur(23)=std(pchip(coil.IU270_time./1000,coil.IU270,timewindow));
        stel_input.sigma_extcur(24)=std(pchip(coil.IU210_time./1000,coil.IU210,timewindow));
        stel_input.sigma_extcur(25)=std(pchip(coil.IU150_time./1000,coil.IU150,timewindow));
        stel_input.sigma_extcur(26)=std(pchip(coil.IU90_time./1000,coil.IU90,timewindow));
        stel_input.sigma_extcur(27)=std(pchip(coil.IU30_time./1000,coil.IU30,timewindow));
        stel_input.sigma_extcur(28)=std(pchip(coil.IL330_time./1000,coil.IL330,timewindow));
        stel_input.sigma_extcur(29)=std(pchip(coil.IL270_time./1000,coil.IL270,timewindow));
        stel_input.sigma_extcur(30)=std(pchip(coil.IL210_time./1000,coil.IL210,timewindow));
        stel_input.sigma_extcur(31)=std(pchip(coil.IL150_time./1000,coil.IL150,timewindow));
        stel_input.sigma_extcur(32)=std(pchip(coil.IL90_time./1000,coil.IL90,timewindow));
        stel_input.sigma_extcur(33)=std(pchip(coil.IL30_time./1000,coil.IL30,timewindow));
        stel_input.sigma_extcur(34)=std(pchip(coil.C79_time./1000,coil.C79,timewindow));
        stel_input.sigma_extcur(35)=std(pchip(coil.C139_time./1000,coil.C139,timewindow));
        stel_input.sigma_extcur(36)=std(pchip(coil.C199_time./1000,coil.C199,timewindow));
        for j=1:length(vmec_input.extcur)-3
            fprintf(fid,'  LEXTCUR_OPT(%2.2d) = T  DEXTCUR_OPT(%2.2d) = 0.01  TARGET_EXTCUR(%2.2d) = %20.10f  SIGMA_EXTCUR(%2.2d) = %20.10e \n',...
                j,j,j,vmec_input.extcur(j),j,abs(stel_input.sigma_extcur(j)));
        end
        for j=34:36
            fprintf(fid,'  LEXTCUR_OPT(%2.2d) = F  DEXTCUR_OPT(%2.2d) = 0.01  TARGET_EXTCUR(%2.2d) = %20.10f  SIGMA_EXTCUR(%2.2d) = %20.10e \n',...
                j,j,j,vmec_input.extcur(j),j,abs(stel_input.sigma_extcur(j)));
        end
    end
    if (luse_mags)
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          DIAGNO OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'  MAGDIAG_COIL = ''/u/slazerso/Sims/DIIID/coils/coils.d3d_efbic''\n');
        for j=1:length(bprobe_names)
            bprobe_val=bprobe.(upper(bprobe_names{j}));
            bprobe_time=bprobe.([upper(bprobe_names{j}) '_TIME']);
            bprobe_sigma=std(pchip(bprobe_time./1000,bprobe_val,timewindow));
            bprobe_val=sum(pchip(bprobe_time./1000,bprobe_val,timewindow))./(ndt+1);
            if bprobe_sigma < 0.0
                bprobe_sigma = 1.0E30;
            end
            if abs(bprobe_sigma) < 0.05*abs(bprobe_val)
                bprobe_sigma = 0.05*abs(bprobe_val);
            end
            if (abs(bprobe_val) < 1.0E-3)
                bprobe_sigma = 1.0E30;
            end
            if (bprobe_val == 0.0)
                bprobe_sigma = 1.0E30;
            end
            fprintf(fid,'  TARGET_BPROBE(%2.2d) = %20.10E  SIGMA_BPROBE(%2.2d) = %20.10E\n',...
                    j,bprobe_val,j,bprobe_sigma);
        end
        for j=1:length(floop_names)
            floop_val=floop.(upper(floop_names{j}));
            floop_time=floop.([upper(floop_names{j}) '_TIME']);
            floop_sigma=std(pchip(floop_time./1000,floop_val,timewindow));
            floop_val=sum(pchip(floop_time./1000,floop_val,timewindow))./(ndt+1);
            if (floop_sigma < 0.0)
                floop_sigma = 1.0E30;
            end
            if (floop_sigma < 0.2*abs(floop_val))
                floop_sigma = 0.2*abs(floop_val);
            end
            if (abs(floop_val) < 1.0E-3)
                floop_sigma = 1.0E30;
            end
            if (floop_val == 0.0)
                floop_sigma = 1.0E30;
            end
            fprintf(fid,'  TARGET_FLUXLOOP(%2.2d) = %20.10E  SIGMA_FLUXLOOP(%2.2d) = %20.10E\n',...
                j,floop_val,j,floop_sigma);
        end
    end
    % Core Thomson
    if (luse_thomson)
        % Output Ne
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          THOMSON (NE) PROFILE OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        % Fix the value of phi (waiting on Ed's response to verify)
        thom_phi_core = (360-thom_phi_core);
        thom_phi_div = (360-thom_phi_div);
        thom_phi_tan = (360-thom_phi_tan);
        %fprintf(fid,'  NNE_PROF = %3.3d\n',length(thom_r_core)+length(thom_r_div)+length(thom_r_tan));
        for j=1:length(thom_r_core)
            fprintf(fid,['  R_NE(%3.3d) = %20.10E  PHI_NE(%3.3d) = %20.10E' ...
                '  Z_NE(%3.3d) = %20.10E  TARGET_NE(%3.3d) = %20.10E'...
                '  SIGMA_NE(%3.3d) = %20.10E\n'],...
                j,thom_r_core(j),j,pi*thom_phi_core(j)/180,...
                j,thom_z_core(j),j,ne_core(j),...
                j,sigma_ne_core(j));
        end
        joff=length(thom_r_core);
        for j=1:length(thom_r_div)
            fprintf(fid,['!  R_NE(%3.3d) = %20.10E  PHI_NE(%3.3d) = %20.10E' ...
                '  Z_NE(%3.3d) = %20.10E  TARGET_NE(%3.3d) = %20.10E'...
                '  SIGMA_NE(%3.3d) = %20.10E\n'],...
                j+joff,thom_r_div(j),j+joff,pi*thom_phi_div(j)/180,...
                j+joff,thom_z_div(j),j+joff,ne_div(j),...
                j+joff,sigma_ne_div(j));
        end
        joff=length(thom_r_core)+length(thom_r_div);
        for j=1:length(thom_r_tan)
            fprintf(fid,['  R_NE(%3.3d) = %20.10E  PHI_NE(%3.3d) = %20.10E' ...
                '  Z_NE(%3.3d) = %20.10E  TARGET_NE(%3.3d) = %20.10E'...
                '  SIGMA_NE(%3.3d) = %20.10E\n'],...
                j+joff,thom_r_tan(j),j+joff,pi*thom_phi_tan(j)/180,...
                j+joff,thom_z_tan(j),j+joff,ne_tan(j),...
                j+joff,sigma_ne_tan(j));
        end
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          THOMSON (TE) PROFILE OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        %fprintf(fid,'  NTE_PROF = %3.3d\n',length(thom_r_core)+length(thom_r_div)+length(thom_r_tan));
        % Output Te
        for j=1:length(thom_r_core)
            if abs(sigma_te_core(j)) < 0.05*abs(te_core(j))
                sigma_te_core(j) = 0.05*abs(te_core(j));
            end
            fprintf(fid,['  R_TE(%3.3d) = %20.10E  PHI_TE(%3.3d) = %20.10E' ...
                '  Z_TE(%3.3d) = %20.10E  TARGET_TE(%3.3d) = %20.10E'...
                '  SIGMA_TE(%3.3d) = %20.10E\n'],...
                j,thom_r_core(j),j,pi*thom_phi_core(j)/180,...
                j,thom_z_core(j),j,te_core(j),...
                j,sigma_te_core(j));
        end
        joff=length(thom_r_core);
        for j=1:length(thom_r_div)
            if abs(sigma_te_div(j)) < 0.05*abs(te_div(j))
                sigma_te_div(j) = 0.05*abs(te_div(j));
            end
            fprintf(fid,['!  R_TE(%3.3d) = %20.10E  PHI_TE(%3.3d) = %20.10E' ...
                '  Z_TE(%3.3d) = %20.10E  TARGET_TE(%3.3d) = %20.10E'...
                '  SIGMA_TE(%3.3d) = %20.10E\n'],...
                j+joff,thom_r_div(j),j+joff,pi*thom_phi_div(j)/180,...
                j+joff,thom_z_div(j),j+joff,te_div(j),...
                j+joff,sigma_te_div(j));
        end
        joff=length(thom_r_core)+length(thom_r_div);
        for j=1:length(thom_r_tan)
            if abs(sigma_te_tan(j)) < 0.05*abs(te_tan(j))
                sigma_te_tan(j) = 0.05*abs(te_tan(j));
            end
            fprintf(fid,['  R_TE(%3.3d) = %20.10E  PHI_TE(%3.3d) = %20.10E' ...
                '  Z_TE(%3.3d) = %20.10E  TARGET_TE(%3.3d) = %20.10E'...
                '  SIGMA_TE(%3.3d) = %20.10E\n'],...
                j+joff,thom_r_tan(j),j+joff,pi*thom_phi_tan(j)/180,...
                j+joff,thom_z_tan(j),j+joff,te_tan(j),...
                j+joff,sigma_te_tan(j));
        end
        j = find(te_core == 0.0,1,'last');
        r = 0.5*(thom_r_core(j+1)+thom_r_core(j-1));
        phi = 0.5*(thom_phi_core(j+1)+thom_phi_core(j-1));
        z = 0.5*(thom_z_core(j+1)+thom_z_core(j-1));
        dz = (thom_z_core(j+1)-thom_z_core(j-1));
        fprintf(fid,['  R_SEPARATRIX(1,1) = %20.10E  PHI_SEPARATRIX(1,1) = %20.10E' ...
            '  Z_SEPARATRIX(1,1) = %20.10E  TARGET_SEPARATRIX(1,1) = %20.10E'...
            '  SIGMA_SEPARATRIX(1,1) = %20.10E\n'],...
            r,pi*phi/180,...
            z,0.0,...
            abs(dz));
        j = find(te_div == 0.0,1,'last');
        r = 0.5*(thom_r_div(j+1)+thom_r_div(j-1));
        phi = 0.5*(thom_phi_div(j+1)+thom_phi_div(j-1));
        z = 0.5*(thom_z_div(j+1)+thom_z_div(j-1));
        dz = 2*(thom_z_div(j+1)-thom_z_div(j-1));
        fprintf(fid,['  R_SEPARATRIX(1,2) = %20.10E  PHI_SEPARATRIX(1,2) = %20.10E' ...
            '  Z_SEPARATRIX(1,2) = %20.10E  TARGET_SEPARATRIX(1,2) = %20.10E'...
            '  SIGMA_SEPARATRIX(1,2) = %20.10E\n'],...
            r,pi*phi/180,...
            z,0.0,...
            abs(dz));
        
    end
    if luse_co2
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          CO2 Laser Interferrometer\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        line_ne_target=sum(pchip(co2_t./1000,co2r0_den,timewindow))./(ndt+1);
        line_ne_sigma=std(pchip(co2_t./1000,co2r0_den,timewindow));
        if line_ne_target > 0
            fprintf(fid,['  R0_NE_LINE(%3.3d) = %20.10E  PHI0_NE_LINE(%3.3d) = %20.10E  Z0_NE_LINE(%3.3d) = %20.10E'...
                '  R1_NE_LINE(%3.3d) = %20.10E  PHI1_NE_LINE(%3.3d) = %20.10E  Z1_NE_LINE(%3.3d) = %20.10E'...
                '  TARGET_NE_LINE(%3.3d) = %20.10E'...
                '  SIGMA_NE_LINE(%3.3d) = %20.10E\n'],...
                1,co2r0_r1,1,co2r0_phi1,1,co2r0_z1,...
                1,co2r0_r2,1,co2r0_phi2,1,co2r0_z2,...
                1,line_ne_target,1,line_ne_sigma);
        end
        line_ne_target=sum(pchip(co2_t./1000,co2v1_den,timewindow))./(ndt+1);
        line_ne_sigma=std(pchip(co2_t./1000,co2v1_den,timewindow));
        if line_ne_target > 0
            fprintf(fid,['  R0_NE_LINE(%3.3d) = %20.10E  PHI0_NE_LINE(%3.3d) = %20.10E  Z0_NE_LINE(%3.3d) = %20.10E'...
                '  R1_NE_LINE(%3.3d) = %20.10E  PHI1_NE_LINE(%3.3d) = %20.10E  Z1_NE_LINE(%3.3d) = %20.10E'...
                '  TARGET_NE_LINE(%3.3d) = %20.10E'...
                '  SIGMA_NE_LINE(%3.3d) = %20.10E\n'],...
                2,co2v1_r1,2,co2v1_phi1,2,co2v1_z1,...
                2,co2v1_r2,2,co2v1_phi2,2,co2v1_z2,...
                2,line_ne_target,2,line_ne_sigma);
        end
        line_ne_target=sum(pchip(co2_t./1000,co2v2_den,timewindow))./(ndt+1);
        line_ne_sigma=std(pchip(co2_t./1000,co2v2_den,timewindow));
        if line_ne_target > 0
            fprintf(fid,['  R0_NE_LINE(%3.3d) = %20.10E  PHI0_NE_LINE(%3.3d) = %20.10E  Z0_NE_LINE(%3.3d) = %20.10E'...
                '  R1_NE_LINE(%3.3d) = %20.10E  PHI1_NE_LINE(%3.3d) = %20.10E  Z1_NE_LINE(%3.3d) = %20.10E'...
                '  TARGET_NE_LINE(%3.3d) = %20.10E'...
                '  SIGMA_NE_LINE(%3.3d) = %20.10E\n'],...
                3,co2v2_r1,3,co2v2_phi1,3,co2v2_z1,...
                3,co2v2_r2,3,co2v2_phi2,3,co2v2_z2,...
                3,line_ne_target,3,line_ne_sigma);
        end
        line_ne_target=sum(pchip(co2_t./1000,co2v3_den,timewindow))./(ndt+1);
        line_ne_sigma=std(pchip(co2_t./1000,co2v3_den,timewindow));
        if line_ne_target > 0
            fprintf(fid,['  R0_NE_LINE(%3.3d) = %20.10E  PHI0_NE_LINE(%3.3d) = %20.10E  Z0_NE_LINE(%3.3d) = %20.10E'...
                '  R1_NE_LINE(%3.3d) = %20.10E  PHI1_NE_LINE(%3.3d) = %20.10E  Z1_NE_LINE(%3.3d) = %20.10E'...
                '  TARGET_NE_LINE(%3.3d) = %20.10E'...
                '  SIGMA_NE_LINE(%3.3d) = %20.10E\n'],...
                4,co2v3_r1,4,co2v3_phi1,4,co2v3_z1,...
                4,co2v3_r2,4,co2v3_phi2,4,co2v3_z2,...
                4,line_ne_target,4,line_ne_sigma);
        end
    end
    if luse_sxr
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          Soft X-Ray Emissions\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        stel_input.zeff_aux_s=vmec_input.am_aux_s;
        dex=find(stel_input.zeff_aux_s == max(stel_input.zeff_aux_s(4:length(stel_input.zeff_aux_s))),1,'first');
        stel_input.zeff_aux_f(1:dex) = 1.0;
        write_namelist_vec(fid,'ZEFF_AUX_S',stel_input.zeff_aux_s(1:dex));
        write_namelist_vec(fid,'ZEFF_AUX_F',stel_input.zeff_aux_f(1:dex));
        % TA Chords
        n=1;
        for j=1:12
            line_sxr_target=sum(pchip(sxr45r1s_time./1000,sxr45r1s(j,:),timewindow))./(ndt+1);
            line_sxr_sigma=std(pchip(sxr45r1s_time./1000,sxr45r1s(j,:),timewindow));
            if line_sxr_target > 0
                fprintf(fid,['  R0_SXR(%3.3d) = %20.10E  PHI0_SXR(%3.3d) = %20.10E  Z0_SXR(%3.3d) = %20.10E'...
                    '  R1_SXR(%3.3d) = %20.10E  PHI1_SXR(%3.3d) = %20.10E  Z1_SXR(%3.3d) = %20.10E'...
                    '  TARGET_SXR(%3.3d) = %20.10E'...
                    '  SIGMA_SXR(%3.3d) = %20.10E\n'],...
                    n,sxr45r1s_r0(j),n,sxr45r1s_phi0(j),n,sxr45r1s_z0(j),...
                    n,sxr45r1s_r1(j),n,sxr45r1s_phi1(j),n,sxr45r1s_z1(j),...
                    n,line_sxr_target,n,line_sxr_sigma);
                n=n+1;
            end
            line_sxr_target=sum(pchip(sxr165r1s_time./1000,sxr165r1s(j,:),timewindow))./(ndt+1);
            line_sxr_sigma=std(pchip(sxr165r1s_time./1000,sxr165r1s(j,:),timewindow));
            if line_sxr_target > 0
                fprintf(fid,['  R0_SXR(%3.3d) = %20.10E  PHI0_SXR(%3.3d) = %20.10E  Z0_SXR(%3.3d) = %20.10E'...
                    '  R1_SXR(%3.3d) = %20.10E  PHI1_SXR(%3.3d) = %20.10E  Z1_SXR(%3.3d) = %20.10E'...
                    '  TARGET_SXR(%3.3d) = %20.10E'...
                    '  SIGMA_SXR(%3.3d) = %20.10E\n'],...
                    n,sxr165r1s_r0(j),n,sxr165r1s_phi0(j),n,sxr165r1s_z0(j),...
                    n,sxr165r1s_r1(j),n,sxr165r1s_phi1(j),n,sxr165r1s_z1(j),...
                    n,line_sxr_target,n,line_sxr_sigma);
                n=n+1;
            end
            line_sxr_target=sum(pchip(sxr195r1s_time./1000,sxr195r1s(j,:),timewindow))./(ndt+1);
            line_sxr_sigma=std(pchip(sxr195r1s_time./1000,sxr195r1s(j,:),timewindow));
            if line_sxr_target > 0
                fprintf(fid,['  R0_SXR(%3.3d) = %20.10E  PHI0_SXR(%3.3d) = %20.10E  Z0_SXR(%3.3d) = %20.10E'...
                    '  R1_SXR(%3.3d) = %20.10E  PHI1_SXR(%3.3d) = %20.10E  Z1_SXR(%3.3d) = %20.10E'...
                    '  TARGET_SXR(%3.3d) = %20.10E'...
                    '  SIGMA_SXR(%3.3d) = %20.10E\n'],...
                    n,sxr195r1s_r0(j),n,sxr195r1s_phi0(j),n,sxr195r1s_z0(j),...
                    n,sxr195r1s_r1(j),n,sxr195r1s_phi1(j),n,sxr195r1s_z1(j),...
                    n,line_sxr_target,n,line_sxr_sigma);
                n=n+1;
            end
        end
        % PA Chords
        for j=1:32
            line_sxr_target=sum(pchip(sxr90p1s_time./1000,sxr90p1s(j,:),timewindow))./(ndt+1);
            line_sxr_sigma=std(pchip(sxr90p1s_time./1000,sxr90p1s(j,:),timewindow));
            if line_sxr_target > 0
                fprintf(fid,['!  R0_SXR(%3.3d) = %20.10E  PHI0_SXR(%3.3d) = %20.10E  Z0_SXR(%3.3d) = %20.10E'...
                    '  R1_SXR(%3.3d) = %20.10E  PHI1_SXR(%3.3d) = %20.10E  Z1_SXR(%3.3d) = %20.10E'...
                    '  TARGET_SXR(%3.3d) = %20.10E'...
                    '  SIGMA_SXR(%3.3d) = %20.10E\n'],...
                    n,sxr90p1s_r0(j),n,sxr90p1s_phi0(j),n,sxr90p1s_z0(j),...
                    n,sxr90p1s_r1(j),n,sxr90p1s_phi1(j),n,sxr90p1s_z1(j),...
                    n,line_sxr_target,n,line_sxr_sigma);
                n=n+1;
            end
            line_sxr_target=sum(pchip(sxr90m1s_time./1000,sxr90m1s(j,:),timewindow))./(ndt+1);
            line_sxr_sigma=std(pchip(sxr90m1s_time./1000,sxr90m1s(j,:),timewindow));
            if line_sxr_target > 0
                fprintf(fid,['!  R0_SXR(%3.3d) = %20.10E  PHI0_SXR(%3.3d) = %20.10E  Z0_SXR(%3.3d) = %20.10E'...
                    '  R1_SXR(%3.3d) = %20.10E  PHI1_SXR(%3.3d) = %20.10E  Z1_SXR(%3.3d) = %20.10E'...
                    '  TARGET_SXR(%3.3d) = %20.10E'...
                    '  SIGMA_SXR(%3.3d) = %20.10E\n'],...
                    n,sxr90m1s_r0(j),n,sxr90m1s_phi0(j),n,sxr90m1s_z0(j),...
                    n,sxr90m1s_r1(j),n,sxr90m1s_phi1(j),n,sxr90m1s_z1(j),...
                    n,line_sxr_target,n,line_sxr_sigma);
                n=n+1;
            end
        end
    end
    if luse_cer
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          CHARGE EXCHANGE (TI) OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        for j=1:length(cer.r)
            if(isempty(cer.ti{j}))
                ti=0.0;
                ti_sig=1.0E30;
            else
                ti=sum(pchip(cer.t{j}./1000,cer.ti{j},timewindow))./(ndt+1);
                ti_sig=std(pchip(cer.t{j}./1000,cer.ti{j},timewindow));
                ti_sig2=sum(pchip(cer.t{j}./1000,cer.sig_ti{j},timewindow))./(ndt+1);
                if (ti < 0.0) || (ti > 1.0E4) || (cer.r(j) >= 2.365)
                    ti = 0.0;
                    ti_sig=1.0E30;
                else
                    ti_sig=max([ti_sig ti_sig2]);
                    fprintf(fid,['  R_TI(%3.3d) = %20.10E  PHI_TI(%3.3d) = %20.10E' ...
                        '  Z_TI(%3.3d) = %20.10E  TARGET_TI(%3.3d) = %20.10E'...
                        '  SIGMA_TI(%3.3d) = %20.10E\n'],...
                        j,cer.r(j),j,cer.phi(j),...
                        j,cer.z(j),j,ti,...
                        j,ti_sig);
                end
            end
        end
    end
    if (luse_mse)
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          MOTIONAL STARK EFFECT (MSE) OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        disp('  - Requesting MSE');
        mds_server.CloseTree(tree,shotnum);
        mds_server.OpenTree(tree_mse,shotnum);  % Open the Tree
        tree_end=tree_mse; % so the last command works
        if shotnum == 126006
            geom = importdata('/Volumes/slazerso/Sims/DIIID/probes/mrz_43.dat');
            num_mse=45;
            mse_str_A='\\MSEP%-2d';
            mse_str_B=[];
            mse_str_C='DIM_OF(\\MSEP%-2d)';
        else
            geom = importdata('/Volumes/slazerso/Sims/DIIID/probes/mrz_44.dat');
            num_mse=69;
            mse_str_A='\\MSE::TOP.ANALYSIS_01:MSEP%2.2d:TWO_K';
            mse_str_B='\\MSE::TOP.ANALYSIS_01:MSEP%2.2d:TWO_K:STDDEV';
            mse_str_C='\\MSE::TOP.ANALYSIS_01:MSEP%2.2d:TWO_K:TIME';
        end
        mse_a1 = geom.data(:,3);
        mse_a2 = geom.data(:,4);
        mse_a3 = geom.data(:,5);
        mse_a4 = geom.data(:,6);
        mse_a5 = geom.data(:,7);
        mse_a6 = geom.data(:,8);
        mse_a7 = geom.data(:,9);
        geom = importdata('/Volumes/slazerso/Sims/DIIID/probes/mse_geo.txt'); % FROM DIII-D Website
        mse_r = geom.data(:,2);
        mse_x = geom.data(:,3);
        mse_y = geom.data(:,4);
        mse_z = geom.data(:,5);
        mse_phi = atan2(mse_y,mse_x);
        mse_phid = 180.*mse_phi/pi;
        mse_phid = mse_phid - 90; 
        mse_phid = 360+mse_phid;
        mse_phi = pi*mse_phid./180.;
        %coil_data=read_coils(coil_file);
        coil_data=coil_biot_prep(coil_data);
        if (lefield)
            stel_input.phi_aux_s=vmec_input.am_aux_s;
            stel_input.phi_aux_f=1.0*stel_input.phi_aux_f;  % Set to constant =1.0 so normalization works
            dex=find(stel_input.phi_aux_s == max(stel_input.phi_aux_s(4:length(stel_input.phi_aux_s))),1,'first');
            write_namelist_vec(fid,'PHI_AUX_S',stel_input.phi_aux_s(1:dex));
            write_namelist_vec(fid,'PHI_AUX_F',stel_input.phi_aux_f(1:dex));
        end
        fprintf(fid,'  LMSE_EXTCUR(1:36) = 36*F\n');
        for j=1:num_mse
            mse_gamma_str=num2str(j,mse_str_A);
            sigma_mse_gamma_str=num2str(j,mse_str_B);
            mse_time_str=num2str(j,mse_str_C);
            % These should work but aren't implemented
            mse_a1_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:A1'];
            mse_a2_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:A2'];
            mse_a3_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:A3'];
            mse_a4_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:A4'];
            mse_a5_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:A5'];
            mse_a6_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:A6'];
            mse_a7_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:A7'];
            mse_r_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:R'];
            mse_z_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:Z'];
            mse_phi_str=['\MSE::TOP.MS' num2str(j,'%2.2d') '.GEOMETRY:PHI_TOR'];
            try
                mse_gamma=mds_server.Value(mse_gamma_str).Double;
                if (~isempty(mse_str_B))
                    sigma_mse_gamma=mds_server.Value(sigma_mse_gamma_str).Double;
                else
                    sigma_mse_gamma=mse_gamma.*0.0;
                end
                mse_time=mds_server.Value(mse_time_str).Double;
                x=mse_r(j)*cos(mse_phi(j));
                y=mse_r(j)*sin(mse_phi(j));
                z=mse_z(j);
                [bx by bz]=coil_biot(coil_data,x,y,z,vmec_input.extcur);
                br=bx*cos(mse_phi(j))+by*sin(mse_phi(j));
                bphi=by*cos(mse_phi(j))-bx*sin(mse_phi(j));
                mse_pol_vac=atan((mse_a1(j)*bz)/(mse_a2(j)*bphi+mse_a3(j)*br+mse_a4(j)*bz));
                gamma=sum(pchip(mse_time./1000,mse_gamma,timewindow))./(ndt+1);
                sig=std(pchip(mse_time./1000,mse_gamma,timewindow))./(ndt+1);
                sig2=sum(pchip(mse_time./1000,sigma_mse_gamma,timewindow))./(ndt+1);
                sigma_mse = 0.05*pi.*gamma./180;
                mse_pol_vac = 0.0;
                fprintf(fid,'  R_MSE(%3d) = %20.10e ',j,mse_r(j));
                fprintf(fid,'  PHI_MSE(%3d) = %20.10e ',j,mse_phi(j));
                fprintf(fid,'  Z_MSE(%3d) = %20.10e ',j,mse_z(j));
                fprintf(fid,'  A1_MSE(%3d) = %20.10e ',j,mse_a1(j));
                fprintf(fid,'  A2_MSE(%3d) = %20.10e ',j,mse_a2(j));
                fprintf(fid,'  A3_MSE(%3d) = %20.10e ',j,mse_a3(j));
                fprintf(fid,'  A4_MSE(%3d) = %20.10e ',j,mse_a4(j));
                fprintf(fid,'  A5_MSE(%3d) = %20.10e ',j,mse_a5(j));
                fprintf(fid,'  A6_MSE(%3d) = %20.10e ',j,mse_a6(j));
                fprintf(fid,'  A7_MSE(%3d) = %20.10e ',j,mse_a7(j));
                fprintf(fid,'  TARGET_MSE(%3d) = %20.10e ',j,pi.*gamma./180);
                fprintf(fid,'  SIGMA_MSE(%3d) = %20.10e ',j,sigma_mse);
                fprintf(fid,'  VAC_MSE(%3d) = %20.10e \n',j,mse_pol_vac);
            catch
                gamma=0.0;
                sigma_mse=1.0E30;
                mse_time=1.0;
                mse_pol_vac=0.0;
            end
        end
    end
    fprintf(fid,'/\n');
    % Handle the DIAGNO input namelist
    if (luse_mags)
        fprintf(fid,'&DIAGNO_IN\n');
        write_namelist_int(fid,'NU',72);
        write_namelist_int(fid,'NV',72);
        fprintf(fid,'  FLUX_DIAG_FILE = ''/u/slazerso/Sims/DIIID/probes/fluxloop_d3d_all_sorted.diagno''\n');
        fprintf(fid,'  BPROBES_FILE = ''/u/slazerso/Sims/DIIID/probes/bprobes_d3d_nonaxi.diagno''\n');
        %fprintf(fid,'  MIRNOV_FILE = ''/u/slazerso/Sims/DIIID/probes''\n');
        %fprintf(fid,'  SEG_ROG_FILE = ''/u/slazerso/Sims/DIIID/probes''\n');
        %fprintf(fid,'  BFIELD_POINTS_FILE = ''/u/slazerso/Sims/DIIID/probes''\n');
        fprintf(fid,'  INT_TYPE = ''simpson''\n');
        write_namelist_int(fid,'INT_STEP',2);
        write_namelist_boo(fid,'LRPHIZ',0);
        write_namelist_flt(fid,'VC_ADAPT_TOL',0.0);
        write_namelist_flt(fid,'VC_ADAPT_REL',5.0E-3);
        write_namelist_boo(fid,'LVC_FIELD',1);
        fprintf(fid,'/\n');
    end
    fclose(fid);
    % Make batch file
    batchname=strtrim([num2str(timedex1(i)*1000.,'%4.4d') '_' user '.batch']);
    fid=fopen(batchname,'w+');
    fprintf(fid,'#!/bin/tcsh\n');  % NOTE TCSH
    fprintf(fid,'#Batchname\n');
    fprintf(fid,'#PBS -N %s\n',strtrim(['DIIID' num2str(shotnum,'%d') '_' num2str(timedex1(i)*1000.,'%4.4d') '_' user ]));
    fprintf(fid,'# Email\n');
    fprintf(fid,'#PBS -m aeb\n');
    switch user
        case 'SAL'
            fprintf(fid,'#PBS -M lazerson@pppl.gov\n');
        case 'NAP'
            fprintf(fid,'#PBS -M npablant@pppl.gov\n');
    end
    fprintf(fid,'# Set Cluster\n');
    %fprintf(fid,'#PBS -q kruskal\n');
    fprintf(fid,'# Nodes, Processors, Memory\n');
    fprintf(fid,'#PBS -l nodes=%-3d\n',numprocs);
    fprintf(fid,'#PBS -l mem=%-5dgb\n',2*numprocs);
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
    fprintf(fid,'set EXEPATH="/u/slazerso/bin/xstelloptv2"\n');
    fprintf(fid,'\n');
    dex=strfind(tfile_name,'.');
    fprintf(fid,'set ARGS="%s"\n',strtrim(tfile_name(dex+1:length(tfile_name))));
    fprintf(fid,'# --- echo the command syntax\n');
    fprintf(fid,'echo "The command syntax for this job is:"\n');
    fprintf(fid,'echo mpirun --mca btl ^openib -np $NPROCS $EXEPATH input.$ARGS\n');
    fprintf(fid,'echo " "\n');
    fprintf(fid,'cd $PBS_O_WORKDIR\n');
    fprintf(fid,'echo -n ''Started job at : '' ; date\n');
    fprintf(fid,'time mpirun --mca btl ^openib -np $NPROCS $EXEPATH input.$ARGS >& log.$ARGS\n');
    fprintf(fid,'rm *_opt* \n');
    fprintf(fid,' \n');
    fprintf(fid,'echo -n ''Ended job at  : '' ; date\n');
    fprintf(fid,'echo " " \n');
    fprintf(fid,'exit\n');
    fclose(fid);
end

% Close the Connection
mds_server.CloseTree(tree_end,shotnum);
mds_server.disconnect(server);

end

