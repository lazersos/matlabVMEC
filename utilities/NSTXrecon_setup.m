function [ output_args ] = NSTXrecon_setup(shotnum,timedex1,timedex2,varargin)
%NSTX_recon_setup(shotnum,t1,t2) Produces a STELLOPT input file for NSTX.
%   Detailed explanation goes here


warning off
% Defaults
filename=['input.NSTX' num2str(shotnum,'%d') '_'];
outputtype='VMEC';
percf='%+7.4f';
perct='%04.4i';
percd='%6.6d';
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
ne_units=1e6;          % Electron Number Density Units
te_units=1e3;          % Electron Temperature Units
ti_units=1e3;          % Ion Temperature Units
vt_units=1e3;          % Rotation Profile Units [km/s]
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
luse_eddy=1;            % Controls use of Eddy Currents
lfull_shot=0;           % Controls execution of full shot reconstructions
numprocs=40;            % Default number of processors
flux_vals=[0.0 0.1 0.25 0.50 0.75 1.0];
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
        elseif strcmp(varargin{i},'user') || strcmp(varargin{i},'USER')
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
                case {'CO2'}
                    luse_co2=1;
                case {'SXR'}
                    luse_sxr=1;
                case {'cer'}
                    luse_cer=1;
                case {'efit'}
                    luse_efit=1;
                case {'eddy'}
                    luse_eddy=1;
                case {'fullshot'}
                    lfull_shot=1;
            end
        end
        i=i+1;
    end
end

% Initialize VMEC
vmec_input=vmec_namelist_init('INDATA');
vmec_input.lasym=1;
vmec_input.nfp=1;
vmec_input.ntor=3;
vmec_input.ncurr=1;
vmec_input.niter=20000;
vmec_input.nstep=200;
vmec_input.nvacskip=6;
vmec_input.mpol=24;
vmec_input.nzeta=72;
vmec_input.ntheta=2*vmec_input.mpol+6;
vmec_input.ns_array=[16  32  64  99];
vmec_input.ftol_array=[1.0E-6  1.0E-8  1.0E-10  1.0E-12];
vmec_input.niter_array=[1000  2000  5000  20000];
vmec_input.delt=1.0;
vmec_input.tcon0=1.0;
vmec_input.phiedge=-0.9;
vmec_input.spres_ped=1.0;
vmec_input.bloat=1.0;
vmec_input.pres_scale=1.5;
vmec_input.am=[1 -1 0 0 -1 1];
%vmec_input.ac=[1 -1 0 0 -1 1];
vmec_input.ac=[1 2 1];
vmec_input.am_aux_s=flux_vals;
vmec_input.am_aux_f=1.0E04.*polyval(fliplr(vmec_input.am),vmec_input.am_aux_s);
vmec_input.ac_aux_s=flux_vals;
vmec_input.ac_aux_f=polyval(fliplr(vmec_input.ac),vmec_input.ac_aux_s);
vmec_input.ai_aux_s=flux_vals;
vmec_input.pmass_type='AKIMA_SPLINE';
vmec_input.pcurr_type='AKIMA_SPLINE_IP';
vmec_input.pcurr_type='two_power';
vmec_input.lfreeb=1;
vmec_input.mgrid_file='/u/slazerso/Sims/NSTX/coils/mgrid_nstx_nphi72.nc';
vmec_input.rbc=zeros(2*vmec_input.ntor+1,vmec_input.mpol);
vmec_input.zbs=zeros(2*vmec_input.ntor+1,vmec_input.mpol);
vmec_input.rbc(vmec_input.ntor+1,1) =  1.7;
vmec_input.rbc(vmec_input.ntor+1,2) =  0.5;
vmec_input.zbs(vmec_input.ntor+1,2) = -1.0;
vmec_input.rbc(vmec_input.ntor+1,3) =  0.05;
vmec_input.zbs(vmec_input.ntor+1,3) = -0.03;
vmec_input.raxis    = [];
vmec_input.zaxis    = [];
vmec_input.raxis(1) = 1.0;
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
if (lfull_shot)
    stel_input.niter_opt=1000;
    stel_input.ftol=1.0E-3;
    stel_input.xtol=1.0E-30;
    stel_input.gtol=1.0E-30;
    stel_input.epsfcn=1.0E-6;
    stel_input.factor=100.0;
    stel_input.lkeep_mins=0;
end
  RBC(1,0+1) = 8.3002527747E-01;  ZBS(1,0+1) = 0.0000000000E+00;
  RBS(1,0+1) = 0.0000000000E+00;  ZBC(1,0+1) = -1.0680834222E-02;
  RBC(1,1+1) = 5.9538381079E-01;  ZBS(1,1+1) = 1.0124985691E+00;
  RBS(1,1+1) = 4.8248985515E-04;  ZBC(1,1+1) = 5.1443526261E-04;
  RBC(1,2+1) = 4.7192755774E-02;  ZBS(1,2+1) = -2.7308877390E-02;
  RBS(1,2+1) = 2.7611255257E-04;  ZBC(1,2+1) = 1.3037472859E-04;
  RBC(1,3+1) = -3.0725772308E-03;  ZBS(1,3+1) = 1.2582732114E-03;
  RBS(1,3+1) = 2.3546089056E-05;  ZBC(1,3+1) = 3.4553275155E-05;
  RBC(1,4+1) = 7.1168564199E-03;  ZBS(1,4+1) = 4.5921408381E-03;
  RBS(1,4+1) = 1.5969232971E-04;  ZBC(1,4+1) = 6.8137948633E-05;
  RBC(1,5+1) = 1.6269729490E-03;  ZBS(1,5+1) = -1.8592654398E-03;
  RBS(1,5+1) = 1.6763224223E-04;  ZBC(1,5+1) = 8.9149842628E-05;
  RBC(1,6+1) = 1.4114687892E-03;  ZBS(1,6+1) = -1.1529497810E-03;
  RBS(1,6+1) = 1.0334076939E-04;  ZBC(1,6+1) = -3.5151219962E-07;
  RBC(1,7+1) = 1.5108419146E-03;  ZBS(1,7+1) = 7.8283644277E-04;
  RBS(1,7+1) = 8.9983231343E-05;  ZBC(1,7+1) = -3.7882946568E-05;
  RBC(1,8+1) = 8.7076089199E-04;  ZBS(1,8+1) = 5.6298956108E-04;
  RBS(1,8+1) = 4.0950427603E-05;  ZBC(1,8+1) = -1.0367630704E-05;
  RBC(1,9+1) = 3.1182425950E-04;  ZBS(1,9+1) = 8.6216540916E-06;
  RBS(1,9+1) = 7.9386823894E-06;  ZBC(1,9+1) = -6.3031582532E-06;
  RBC(1,10+1) = 1.4723430977E-05;  ZBS(1,10+1) = 2.0649390917E-04;
  RBS(1,10+1) = 5.9592699700E-06;  ZBC(1,10+1) = -7.5376420940E-06;
  RBC(1,11+1) = -8.2079887792E-05;  ZBS(1,11+1) = 3.4523623859E-04;
  RBS(1,11+1) = 6.5505755898E-07;  ZBC(1,11+1) = 9.0383169962E-07;
  vmec_input.rbc(vmec_input.ntor+1,1:12) = RBC(1,1:12);
  vmec_input.rbs(vmec_input.ntor+1,1:12) = RBS(1,1:12);
  vmec_input.zbc(vmec_input.ntor+1,1:12) = ZBC(1,1:12);
  vmec_input.zbs(vmec_input.ntor+1,1:12) = ZBS(1,1:12);
%vmec_input.rbc(vmec_input.ntor+1,1) =  8.3E-01;
%vmec_input.rbc(vmec_input.ntor+1,2) =  0.5;
%vmec_input.zbs(vmec_input.ntor+1,2) = -1.0;
%vmec_input.rbc(vmec_input.ntor+1,3) =  0.05;
%vmec_input.zbs(vmec_input.ntor+1,3) = -0.03;


% Open a connection to the D3D server
import MdsPlus.*;   % So we can call MdsPlus java calls
server = 'skylark.pppl.gov';  % Server you wish to connect to
server_port = 8501;        % Server port
mds_server=MdsPlus(server,server_port);  % Open the connection

% Open Tree
%tree_mse = 'MSE';              % Tree you wish to open
%tree_end = tree;
%shot_num=142603;           % Shot Number
%mds_server.OpenTree(tree,shotnum);  % Open the Tree

% Now get the coil_currents
disp('  - Requesting Coil Currents');
% Switch Tree
tree = 'NSTX';
mds_server.OpenTree(tree,shotnum);  % Open the Tree
coil_names={'TF','OH','PF1AU','PF2U','PF3U','PF4','PF5','PF3L','PF2L','PF1AL',...
    'PF1B','RWM1','RWM2','RWM3','RWM4','RWM5','RWM6'};
for i=1:length(coil_names)
    coil_str = ['\ENGINEERING::I' coil_names{i}];
    coil_time_str = ['DIM_OF(' coil_str ')'];
    try
        coil.(coil_names{i}) = mds_server.Value(coil_str).Double;
        coil.([coil_names{i} '_time']) = mds_server.Value(coil_time_str).Double;
        disp(['       - ' coil_names{i} ' OK']);
    catch
        coil.(coil_names{i}) = [0.0 0.0];
        coil.([coil_names{i} '_time']) = [-bigno bigno];
        disp(['       - ' coil_names{i} ' FAIL']);
    end
end
% Get Eddy
if (luse_eddy)
    disp('  - Requesting Eddy Currents');
    num_coils(1:17) = 0.0;
    num_coils(18:52) = [135  135   40  150   80 ...
                        125   75  280  340  220 ...
                        150   90   90   90   90 ...
                        150  220  340  250  100 ...
                         80  150   40  135  135 ...
                        290  290   40   40   55 ...
                         55   55   55   40   40];
    vmec_input.mgrid_file='/p/pies/Lazerson/Sims/MGRIDS/mgrid_nstx_060413_errf_nphi48.nc';
    try
        eddy = mds_server.Value('\EFIT01::CCBRSP').Double;
        gtime = mds_server.Value('\EFIT01::GTIME').Double;
    catch
        eddy = [0.0 0.0];
        gtime = [-bigno bigno];
    end
    niter = length(gtime);
    neddy = length(eddy)./niter;
    eddy = reshape(eddy,[neddy niter]);
    for i=18:neddy
        eddy(i,:)=smooth(eddy(i,:),8,'rlowess')./num_coils(i);
    end
end
% Get Net toroidal Current
disp('  - Requesting Toroidal Current');
% Switch Tree
mds_server.CloseTree(tree,shotnum);
tree = 'NSTX';
mds_server.OpenTree(tree,shotnum);  % Open the Tree
ip_str='\ENGINEERING::IP1';
ip_time_str='DIM_OF(\ENGINEERING::IP1)';
ip=mds_server.Value(ip_str).Double;
ip_time=mds_server.Value(ip_time_str).Double;
% Get Magneticts
if (luse_mags)
    % Switch Tree
    mds_server.CloseTree(tree,shotnum);
    tree = 'NSTX';
    mds_server.OpenTree(tree,shotnum);  % Open the Tree
    disp('  - Requesting Magnetics');
    floop_names={'FLEVVL2','FLEVVL3','FLEVVL4','FLEVVL5','FLEVVL6','FLEVVL7',...
        'FLEVVU2','FLEVVU3','FLEVVU4','FLEVVU5','FLEVVU6','FLEVVU7',...
        'FLIVVL1','FLIVVL2','FLIVVU1','FLIVVU2',...
        'FLOBDL1','FLOBDL2','FLOBDL3','FLOBDU1','FLOBDU2','FLOBDU3',...
        'FLOHL2','FLOHL4','FLOHM','FLOHU2','FLOHU4',...
        'FLPF1AL2','FLPF1AL4','FLPF1AU2','FLPF1AU4','FLPF1BL',...
        'FLPFAB1U1','FLPFAB1U2',...
        'FLPPPL1','FLPPPL2','FLPPPL3','FLPPPL4','FLPPPU1','FLPPPU2','FLPPPU3','FLPPPU4',...
        'FLSPPL1','FLSPPL2','FLSPPL3','FLSPPL4','FLSPPU1','FLSPPU2','FLSPPU3','FLSPPU4'};
    bprobe_names={'1DMCSCL1','1DMCSCL2','1DMCSCL3','1DMCSCL4','1DMCSCL5','1DMCSCL6',...
        '1DMCSCU1','1DMCSCU2','1DMCSCU3','1DMCSCU4','1DMCSCU5','1DMCSCU6',...
        '2DMCSCL2N','2DMCSCL2T','2DMCSCU1N','2DMCSCU1T',...
        '2DMIBDHL5N','2DMIBDHL5T','2DMIBDHL6N','2DMIBDHL6T',...
        '2DMIBDHU5N','2DMIBDHU5T','2DMIBDHU6N','2DMIBDHU6T',...
        '2DMIBDVL1N','2DMIBDVL3N','2DMIBDVU1N','2DMIBDVU3N',...
        '2DMOBDL1N','2DMOBDL1T','2DMOBDL2N','2DMOBDL2T','2DMOBDL3N','2DMOBDL3T','2DMOBDL4N','2DMOBDL4T','2DMOBDL5N','2DMOBDL5T',...
        '2DMOBDU1N','2DMOBDU1T','2DMOBDU2N','2DMOBDU2T','2DMOBDU3N','2DMOBDU3T','2DMOBDU4N','2DMOBDU4T','2DMOBDU5N','2DMOBDU5T',...
        'HCDOBDIRLAB','HCDOBDIRLCD','HCDOBDIRLEF','HCDOBDIRLGH','HCDOBDIRLIJ','HCDOBDIRLKL',...
        'HCDOBDORLAB','HCDOBDORLCD','HCDOBDORLEF','HCDOBDORLGH','HCDOBDORLIJ','HCDOBDORLKL',...
        'L1DMPPPGL1','L1DMPPPGL2','L1DMPPPGL3','L1DMPPPGL4','L1DMPPPGL5','L1DMPPPGL6','L1DMPPPGL7','L1DMPPPGL8',...
        'L1DMPPPGU1','L1DMPPPGU2','L1DMPPPGU3','L1DMPPPGU4','L1DMPPPGU5','L1DMPPPGU6','L1DMPPPGU7','L1DMPPPGU8',...
        'L1DMSPPGL1','L1DMSPPGL2','L1DMSPPGL3','L1DMSPPGL4','L1DMSPPGL5','L1DMSPPGL6','L1DMSPPGL7',...
        'L1DMSPPGU1','L1DMSPPGU2','L1DMSPPGU3','L1DMSPPGU4','L1DMSPPGU5','L1DMSPPGU6','L1DMSPPGU7',...
        };
    disp('       - B Probes');
    bprobe=[];
    for j=1:length(bprobe_names)
        bprobe_str=['\OPERATIONS::B_' upper(bprobe_names{j})];
        bprobe_time_str=['DIM_OF(' bprobe_str ')'];
        try
            bprobe.(['B' upper(bprobe_names{j})])=smooth(mds_server.Value(bprobe_str).Double,100,'moving');
            bprobe.(['B' upper(bprobe_names{j}) '_TIME'])=mds_server.Value(bprobe_time_str).Double;
        catch temp_err
            bprobe.(['B' upper(bprobe_names{j})])=[0.0 0.0];
            bprobe.(['B' upper(bprobe_names{j}) '_TIME'])=[-1E30 1.0E30];
        end
    end
    disp('       - Flux Loops');
    floop=[];
    floop.PSF1A=0;  % Do this because PSF1A is the offset loop and read first.
    for j=1:length(floop_names)
        floop_str=['\OPERATIONS::F_' upper(floop_names{j})];
        floop_time_str=['DIM_OF(' floop_str ')'];
        try
            floop.(upper(floop_names{j}))=smooth(mds_server.Value(floop_str).Double,100,'moving');
            floop.([upper(floop_names{j}) '_TIME'])=mds_server.Value(floop_time_str).Double;
        catch 
            disp(['       - ' floop_names{j} ' FAIL']);
            floop.(upper(floop_names{j}))=[0.0 0.0];
            floop.([upper(floop_names{j}) '_TIME'])=[-1E30 1.0E30];
        end
    end
end
% Get C02 Data
if (luse_co2)
    % Switch Tree
    mds_server.CloseTree(tree,shotnum);
    tree = 'MICROWAVE';
    mds_server.OpenTree(tree,shotnum);  % Open the Tree
    disp('  - Requesting CO2 Inteferrometer Data');
    for i=1:7
        co2_den_str = ['\NE_C' num2str(i,'%1d') '_INTEG'];
        try
            co2_den{i} = mds_server.Value(co2r0_den_str).Double.*scale_co2;
        catch
        end
    end
end
% Get SXR Data
if (luse_sxr)
end
% Get Thomson
if (luse_thomson)
    % Switch Tree
    mds_server.CloseTree(tree,shotnum);
    tree = 'NSTX';
    mds_server.OpenTree(tree,shotnum);  % Open the Tree
    disp('  - Requesting Core Thomson');
    thomson_time_core_str='\TS_TIMES';
    thomson_ne_core_str='\activespec::nef';
    thomson_signe_core_str='\activespec::dne';
    thomson_te_core_str='\activespec::tef';
    thomson_sigte_core_str='\activespec::dte';
    thomson_valid_core_str='\activespec::TS_BEST:VALID';
    %thom_r_core=mds_server.Value(thomson_r_core_str).Double;
    thom_r_core=[2.7598384023e-01    4.6741470695e-01    6.8640893698e-01...
        9.1536599398e-01    1.0880001783e+00    1.2516412735e+00...
        1.3306819201e+00    1.4202042818e+00    1.4781734943e+00...
        1.5299917459e+00    3.8919675350e-01    6.1787593365e-01...
        8.0930978060e-01    1.0074237585e+00    1.1916993856e+00...
        1.3056254387e+00    1.3774180412e+00    1.4595996141e+00...
        1.5133252144e+00    1.5617003441e+00    5.4459947348e-01...
        7.5014024973e-01    1.1591422558e+00    1.2224814892e+00...
        1.2793152332e+00    1.3545830250e+00    1.3992675543e+00...
        1.4358093739e+00    1.4447805882e+00    1.4960665703e+00];
    thom_z_core=0.0*thom_r_core;
    thom_phi_core=[1.0787397766e+02    1.5511195374e+02    1.6747312927e+02...
        1.7358239746e+02    1.7641368103e+02    1.7835482788e+02...
        1.7911724854e+02    1.7987615967e+02    1.8031764221e+02...
        1.8068347168e+02    1.4625003052e+02    1.6467547607e+02...
        1.7121081543e+02    1.7521780396e+02    1.7770741272e+02...
        1.7888578796e+02    1.7952601624e+02    1.8018008423e+02...
        1.8056858826e+02    1.8089515686e+02    1.6079321289e+02...
        1.6957609558e+02    1.7732691956e+02    1.7804801941e+02...
        1.7863273621e+02    1.7932987976e+02    1.7970755005e+02...
        1.8003327942e+02    1.8003327942e+02    1.8044688416e+02];
    [thom_r_core, idex]=sort(thom_r_core);
    thom_phi_core = thom_phi_core(idex);
    thom_phi_core = thom_phi_core+90;
    %thomson_r_core = mds_server.Value('\RADII').Double;
    thom_time_core=mds_server.Value(thomson_time_core_str).Double;
    thom_ne_core=mds_server.Value(thomson_ne_core_str).Double;
    thom_signe_core=mds_server.Value(thomson_signe_core_str).Double;
    thom_te_core=mds_server.Value(thomson_te_core_str).Double;
    thom_sigte_core=mds_server.Value(thomson_sigte_core_str).Double;
    thom_valid = mds_server.Value(thomson_valid_core_str).Double;
    % Filter
    thom_ne_core(thom_valid<0) = thom_ne_core(thom_valid<0).*0.0;
    thom_te_core(thom_valid<0) = thom_te_core(thom_valid<0).*0.0;
    thom_signe_core(thom_valid<0) = thom_signe_core(thom_valid<0).*0.0;
    thom_sigte_core(thom_valid<0) = thom_sigte_core(thom_valid<0).*0.0;
    %
    thom_ne_core = reshape(thom_ne_core,[length(thom_time_core) length(thom_ne_core)/length(thom_time_core)]);
    thom_te_core = reshape(thom_te_core,[length(thom_time_core) length(thom_te_core)/length(thom_time_core)]);
    thom_signe_core = reshape(thom_signe_core,[length(thom_time_core) length(thom_signe_core)/length(thom_time_core)]);
    thom_sigte_core = reshape(thom_sigte_core,[length(thom_time_core) length(thom_sigte_core)/length(thom_time_core)]);
end
if luse_cer
    % Switch Tree
    mds_server.CloseTree(tree,shotnum);
    tree = 'NSTX';
    mds_server.OpenTree(tree,shotnum);  % Open the Tree
    disp('  - Requesting Charge Exchange');
    chers_time_str='\ACTIVESPEC::CHERS_BEST:TIME';
    chers_ti_str='\ACTIVESPEC::CHERS_BEST:TI';
    chers_dti_str='\ACTIVESPEC::CHERS_BEST:DTI';
    chers_vt_str='\ACTIVESPEC::CHERS_BEST:VT';
    chers_dvt_str='\ACTIVESPEC::CHERS_BEST:VT';
    chers_valid_str='\ACTIVESPEC::CHERS_BEST:VALID';
    chers_time=mds_server.Value(chers_time_str).Double;
    chers_ti=mds_server.Value(chers_ti_str).Double;
    chers_ti_err=mds_server.Value(chers_dti_str).Double;
    chers_vt=mds_server.Value(chers_vt_str).Double;
    chers_vt_err=mds_server.Value(chers_dvt_str).Double;
    chers_valid=mds_server.Value(chers_valid_str).Double;
    % FILTER
    chers_ti(chers_valid == 0) = 0.0;
    chers_ti_err(chers_valid == 0) = bigno;
    % Get BEAM and Geometry
    mds_server.CloseTree(tree,shotnum);
    tree = 'NBI';
    mds_server.OpenTree(tree,shotnum);  % Open the Tree
    nbi_powera_str = '\N1A_POWER';
    nbi_powerb_str = '\N1B_POWER';
    nbi_powerc_str = '\N1C_POWER';
    nbi_powera=max(mds_server.Value(nbi_powera_str).Double);
    nbi_powerb=max(mds_server.Value(nbi_powerb_str).Double);
    nbi_powerc=max(mds_server.Value(nbi_powerc_str).Double);
    [~, dex]=max([nbi_powera nbi_powerb nbi_powerc]);
    if (dex ==1)
        chers_r=[9.14061E+01  9.40513E+01  9.70520E+01  9.92994E+01...
            1.02670E+02  1.05087E+02  1.08039E+02  1.11660E+02...
            1.14145E+02  1.16434E+02  1.18882E+02  1.21226E+02...
            1.22876E+02  1.24537E+02  1.26445E+02  1.27910E+02...
            1.29503E+02  1.31015E+02  1.32219E+02  1.33390E+02...
            1.34342E+02  1.35526E+02  1.36663E+02  1.37508E+02...
            1.38419E+02  1.38983E+02  1.39594E+02  1.40201E+02...
            1.40767E+02  1.41346E+02  1.41872E+02  1.42466E+02...
            1.43068E+02  1.43671E+02  1.44171E+02  1.44644E+02...
            1.45162E+02  1.45787E+02  1.46354E+02  1.46926E+02...
            1.47443E+02  1.47954E+02  1.48907E+02  1.49796E+02...
            1.50543E+02  1.51379E+02  1.52553E+02  1.53534E+02...
            1.54661E+02  1.55917E+02  1.57285E+02];
        chers_phi=[1.02732E+02   1.00887E+02   9.89859E+01   9.76754E+01   9.58659E+01...
            9.46694E+01   9.33072E+01   9.17678E+01   9.07851E+01   8.99281E+01...
            8.90576E+01   8.82653E+01   8.77298E+01   8.72086E+01   8.66306E+01...
            8.62012E+01   8.57472E+01   8.53286E+01   8.50035E+01   8.46938E+01...
            8.44470E+01   8.41455E+01   8.38619E+01   8.36546E+01   8.34346E+01...
            8.33000E+01   8.31556E+01   8.30135E+01   8.28826E+01   8.27497E+01...
            8.26301E+01   8.24962E+01   8.23618E+01   8.22286E+01   8.21190E+01...
            8.20163E+01   8.19045E+01   8.17710E+01   8.16510E+01   8.15309E+01...
            8.14233E+01   8.13178E+01   8.11234E+01   8.09445E+01   8.07960E+01...
            8.06318E+01   8.04046E+01   8.02178E+01   8.00066E+01   7.97751E+01...
            7.95277E+01];
        chers_z=chers_r.*0.0;
    elseif (dex == 2)
        chers_r=[9.14316E+01   9.41709E+01   9.72439E+01   9.95265E+01   1.02927E+02...
            1.05351E+02   1.08301E+02   1.11905E+02   1.14373E+02   1.16642E+02...
            1.19066E+02   1.21387E+02   1.23021E+02   1.24665E+02   1.26553E+02...
            1.28003E+02   1.29581E+02   1.31078E+02   1.32271E+02   1.33433E+02...
            1.34376E+02   1.35551E+02   1.36680E+02   1.37519E+02   1.38423E+02...
            1.38984E+02   1.39591E+02   1.40195E+02   1.40758E+02   1.41334E+02...
            1.41857E+02   1.42449E+02   1.43049E+02   1.43649E+02   1.44148E+02...
            1.44619E+02   1.45136E+02   1.45759E+02   1.46325E+02   1.46896E+02...
            1.47412E+02   1.47923E+02   1.48875E+02   1.49764E+02   1.50511E+02...
            1.51348E+02   1.52525E+02   1.53509E+02   1.54639E+02   1.55901E+02...
            1.57277E+02];
        chers_phi=[9.76994E+01   9.63006E+01   9.48544E+01   9.38546E+01   9.24688E+01...
            9.15487E+01   9.04968E+01   8.93022E+01   8.85360E+01   8.78652E+01...
            8.71814E+01   8.65566E+01   8.61330E+01   8.57197E+01   8.52601E+01...
            8.49177E+01   8.45550E+01   8.42197E+01   8.39589E+01   8.37099E+01...
            8.35112E+01   8.32680E+01   8.30389E+01   8.28713E+01   8.26931E+01...
            8.25840E+01   8.24668E+01   8.23515E+01   8.22450E+01   8.21370E+01...
            8.20396E+01   8.19305E+01   8.18210E+01   8.17123E+01   8.16228E+01...
            8.15389E+01   8.14475E+01   8.13383E+01   8.12400E+01   8.11416E+01...
            8.10534E+01   8.09668E+01   8.08071E+01   8.06600E+01   8.05378E+01...
            8.04025E+01   8.02150E+01   8.00607E+01   7.98860E+01   7.96942E+01...
            7.94889E+01];
        chers_z=chers_r.*0.0;
    elseif (dex == 3)
        chers_r=[9.20307E+01   9.47718E+01   9.78312E+01   1.00095E+02   1.03458E+02...
            1.05850E+02   1.08756E+02   1.12302E+02   1.14729E+02   1.16961E+02...
            1.19345E+02   1.21629E+02   1.23237E+02   1.24856E+02   1.26716E+02...
            1.28145E+02   1.29702E+02   1.31181E+02   1.32359E+02   1.33507E+02...
            1.34440E+02   1.35603E+02   1.36720E+02   1.37552E+02   1.38449E+02...
            1.39005E+02   1.39607E+02   1.40207E+02   1.40765E+02   1.41337E+02...
            1.41857E+02   1.42446E+02   1.43042E+02   1.43639E+02   1.44135E+02...
            1.44605E+02   1.45120E+02   1.45740E+02   1.46304E+02   1.46874E+02...
            1.47389E+02   1.47898E+02   1.48849E+02   1.49737E+02   1.50485E+02...
            1.51322E+02   1.52500E+02   1.53486E+02   1.54619E+02   1.55886E+02...
            1.57269E+02];
        chers_phi=[9.33240E+01   9.22994E+01   9.12352E+01   9.04963E+01   8.94674E+01...
            8.87809E+01   8.79926E+01   8.70924E+01   8.65121E+01   8.60022E+01...
            8.54803E+01   8.50018E+01   8.46764E+01   8.43580E+01   8.40031E+01...
            8.37380E+01   8.34565E+01   8.31958E+01   8.29925E+01   8.27982E+01...
            8.26429E+01   8.24525E+01   8.22729E+01   8.21413E+01   8.20012E+01...
            8.19154E+01   8.18231E+01   8.17322E+01   8.16482E+01   8.15630E+01...
            8.14861E+01   8.13999E+01   8.13132E+01   8.12272E+01   8.11563E+01...
            8.10897E+01   8.10172E+01   8.09305E+01   8.08525E+01   8.07742E+01...
            8.07041E+01   8.06352E+01   8.05079E+01   8.03905E+01   8.02930E+01...
            8.01849E+01   8.00349E+01   7.99112E+01   7.97710E+01   7.96170E+01...
            7.94518E+01];
        chers_z=chers_r.*0.0;
    end
    % Now reformulate arrays
    n1 = length(chers_r);
    n2 = length(chers_time);
    chers_ti = ti_units*reshape(chers_ti,n1,n2);
    chers_ti_err = ti_units*reshape(chers_ti_err,n1,n2);
    chers_vt = vt_units*reshape(chers_vt,n1,n2);
    chers_vt_err = vt_units*reshape(chers_vt_err,n1,n2);
    vmec_input.gamma=5./3.;
    vmec_input.bcrit=0.3;
end
if luse_efit
    mds_server.CloseTree(tree,shotnum);
    tree = 'NSTX';
    mds_server.OpenTree(tree,shotnum);  % Open the Tree
    disp('  - Requesting EFIT');
    efit_nbdry = mds_server.Value('\EFIT01::NBDRY').Double;
    efit_nmass = mds_server.Value('\EFIT01::NMASS').Double;
    efit_gtime = mds_server.Value('\EFIT01::GTIME').Double;
    efit_zbdry = mds_server.Value('\EFIT01::ZBDRY').Double;
    efit_rbdry = mds_server.Value('\EFIT01::RBDRY').Double;
    efit_ffp = mds_server.Value('\EFIT01::FFPRIM').Double;
    efit_pp = mds_server.Value('\EFIT01::PPRIME').Double;
    efit_pres = mds_server.Value('\EFIT01::PRES').Double;
    efit_fpol = mds_server.Value('\EFIT01::FPOL').Double;
    efit_qpsi = mds_server.Value('\EFIT01::QPSI').Double;
    efit_psin = mds_server.Value('\EFIT01::PSIN').Double;
    efit_psi = mds_server.Value('\EFIT01::PSI').Double;
    efit_rax = mds_server.Value('\EFIT01::R0').Double;
    efit_zax = mds_server.Value('\EFIT01::Z0').Double;
    efit_btorvac = mds_server.Value('\EFIT01::BT0VAC').Double;
    efit_betap = mds_server.Value('\EFIT01::BETAP').Double;
    efit_betat = mds_server.Value('\EFIT01::BETAT').Double;
    i1=1;
    j1=1;
    nbdry = round(max(efit_nbdry));
    ntime = length(efit_gtime) ;
    nrho  = length(efit_ffp)/ntime;
    efit_ffp = reshape(efit_ffp, [nrho ntime]);
    efit_pres = reshape(efit_pres, [nrho ntime]);
    efit_fpol = reshape(efit_fpol, [nrho ntime]);
    efit_qpsi = reshape(efit_qpsi, [nrho ntime]);
    efit_psin = reshape(efit_psin, [nrho ntime]);
    efit_psi = reshape(efit_psi, [nrho ntime]);
    efit_pp = reshape(efit_pp, [nrho ntime]);
    efit_btorvac = reshape(efit_btorvac, [1 ntime]);
    efit_betap = reshape(efit_betap, [1 ntime]);
    efit_betat = reshape(efit_betat, [1 ntime]);
    for i=1:length(efit_nbdry)
        i2=i1+nbdry-1;
        efit_r(1:nbdry,i) = efit_rbdry(i1:i2);
        efit_z(1:nbdry,i) = efit_zbdry(i1:i2);
        i1 = i2 + 1;
    end
end
if luse_mse
    % Switch Tree
    mds_server.CloseTree(tree,shotnum);
    tree = 'MSE';
    mds_server.OpenTree(tree,shotnum);  % Open the Tree
    % MDS+ Strings
    mse_r1_str='.MSE_CIF.ANALYSIS.GEOMETRY:RADIUS';
    mse_r2_str='.MSE_CIF.ANALYSIS:RADIUS';
    mse_alpha_str='.MSE_CIF.ANALYSIS.GEOMETRY:ALPHA';
    mse_omega_str='.MSE_CIF.ANALYSIS.GEOMETRY:OMEGA';
    %mse_g_str='.MSE_CIF.ANALYSIS.GEOMETRY:G';
    mse_pol_str='\MSE_CIF_PA';
    %mse_pol_str='\MSE_CIF_PA_CORR_ER';
    mse_time_str='.MSE_CIF.ANALYSIS:TIME';
    % Get data
    mse_time=mds_server.Value(mse_time_str).Double;
    mse_r1=mds_server.Value(mse_r1_str).Double;
    mse_r2=mds_server.Value(mse_r2_str).Double;
    mse_alpha=mds_server.Value(mse_alpha_str).Double;
    mse_omega=mds_server.Value(mse_omega_str).Double;
    %mse_g=mds_server.Value(mse_g_str).Double;
    mse_pol=mds_server.Value(mse_pol_str).Double;
    % Adjust the data
    for l = 1:length(mse_r2)
        dex = find(mse_r1 == mse_r2(l),1,'first');
        mse_r_temp(l) = mse_r1(dex);
        mse_alpha_temp(l) = mse_alpha(dex);
        mse_omega_temp(l) = mse_omega(dex);
    end
    mse_r = mse_r_temp;
    mse_alpha = mse_alpha_temp;
    mse_omega = mse_omega_temp;
    mse_pol=reshape(mse_pol,[length(mse_r) length(mse_time)]);
    mse_z = mse_r.*0.0;
    % Calculate coefs
    mse_alpha=-mse_alpha;  % This gets things correct (not sure why)
    mse_theta=mse_alpha.*0.0;
    mse_a1=-cosd(mse_alpha+mse_omega);
    mse_a2=sind(mse_alpha).*cosd(mse_theta);
    mse_a3=cosd(mse_alpha).*cosd(mse_theta);
    mse_a4=sind(mse_alpha+mse_omega).*sind(mse_theta);
    mse_a5=0.0.*mse_a4;
    mse_a6=0.0.*mse_a4;
    mse_a7=0.0.*mse_a4;
    mds_server.CloseTree(tree,shotnum);
    % Now get the toroidal angles
    tree = 'NBI';
    mds_server.OpenTree(tree,shotnum);  % Open the Tree
    nbi_powera_str = '\N1A_POWER';
    nbi_powerb_str = '\N1B_POWER';
    nbi_powerc_str = '\N1C_POWER';
    nbi_powera=max(mds_server.Value(nbi_powera_str).Double);
    nbi_powerb=max(mds_server.Value(nbi_powerb_str).Double);
    nbi_powerc=max(mds_server.Value(nbi_powerc_str).Double);
    [~, dex]=max([nbi_powera nbi_powerb nbi_powerc]);
    if (dex ==1)
        chers_r=[9.14061E+01  9.40513E+01  9.70520E+01  9.92994E+01...
            1.02670E+02  1.05087E+02  1.08039E+02  1.11660E+02...
            1.14145E+02  1.16434E+02  1.18882E+02  1.21226E+02...
            1.22876E+02  1.24537E+02  1.26445E+02  1.27910E+02...
            1.29503E+02  1.31015E+02  1.32219E+02  1.33390E+02...
            1.34342E+02  1.35526E+02  1.36663E+02  1.37508E+02...
            1.38419E+02  1.38983E+02  1.39594E+02  1.40201E+02...
            1.40767E+02  1.41346E+02  1.41872E+02  1.42466E+02...
            1.43068E+02  1.43671E+02  1.44171E+02  1.44644E+02...
            1.45162E+02  1.45787E+02  1.46354E+02  1.46926E+02...
            1.47443E+02  1.47954E+02  1.48907E+02  1.49796E+02...
            1.50543E+02  1.51379E+02  1.52553E+02  1.53534E+02...
            1.54661E+02  1.55917E+02  1.57285E+02];
        chers_phi=[1.02732E+02   1.00887E+02   9.89859E+01   9.76754E+01   9.58659E+01...
            9.46694E+01   9.33072E+01   9.17678E+01   9.07851E+01   8.99281E+01...
            8.90576E+01   8.82653E+01   8.77298E+01   8.72086E+01   8.66306E+01...
            8.62012E+01   8.57472E+01   8.53286E+01   8.50035E+01   8.46938E+01...
            8.44470E+01   8.41455E+01   8.38619E+01   8.36546E+01   8.34346E+01...
            8.33000E+01   8.31556E+01   8.30135E+01   8.28826E+01   8.27497E+01...
            8.26301E+01   8.24962E+01   8.23618E+01   8.22286E+01   8.21190E+01...
            8.20163E+01   8.19045E+01   8.17710E+01   8.16510E+01   8.15309E+01...
            8.14233E+01   8.13178E+01   8.11234E+01   8.09445E+01   8.07960E+01...
            8.06318E+01   8.04046E+01   8.02178E+01   8.00066E+01   7.97751E+01...
            7.95277E+01];
    elseif (dex == 2)
        chers_r=[9.14316E+01   9.41709E+01   9.72439E+01   9.95265E+01   1.02927E+02...
            1.05351E+02   1.08301E+02   1.11905E+02   1.14373E+02   1.16642E+02...
            1.19066E+02   1.21387E+02   1.23021E+02   1.24665E+02   1.26553E+02...
            1.28003E+02   1.29581E+02   1.31078E+02   1.32271E+02   1.33433E+02...
            1.34376E+02   1.35551E+02   1.36680E+02   1.37519E+02   1.38423E+02...
            1.38984E+02   1.39591E+02   1.40195E+02   1.40758E+02   1.41334E+02...
            1.41857E+02   1.42449E+02   1.43049E+02   1.43649E+02   1.44148E+02...
            1.44619E+02   1.45136E+02   1.45759E+02   1.46325E+02   1.46896E+02...
            1.47412E+02   1.47923E+02   1.48875E+02   1.49764E+02   1.50511E+02...
            1.51348E+02   1.52525E+02   1.53509E+02   1.54639E+02   1.55901E+02...
            1.57277E+02];
        chers_phi=[9.76994E+01   9.63006E+01   9.48544E+01   9.38546E+01   9.24688E+01...
            9.15487E+01   9.04968E+01   8.93022E+01   8.85360E+01   8.78652E+01...
            8.71814E+01   8.65566E+01   8.61330E+01   8.57197E+01   8.52601E+01...
            8.49177E+01   8.45550E+01   8.42197E+01   8.39589E+01   8.37099E+01...
            8.35112E+01   8.32680E+01   8.30389E+01   8.28713E+01   8.26931E+01...
            8.25840E+01   8.24668E+01   8.23515E+01   8.22450E+01   8.21370E+01...
            8.20396E+01   8.19305E+01   8.18210E+01   8.17123E+01   8.16228E+01...
            8.15389E+01   8.14475E+01   8.13383E+01   8.12400E+01   8.11416E+01...
            8.10534E+01   8.09668E+01   8.08071E+01   8.06600E+01   8.05378E+01...
            8.04025E+01   8.02150E+01   8.00607E+01   7.98860E+01   7.96942E+01...
            7.94889E+01];
    elseif (dex == 3)
        chers_r=[9.20307E+01   9.47718E+01   9.78312E+01   1.00095E+02   1.03458E+02...
            1.05850E+02   1.08756E+02   1.12302E+02   1.14729E+02   1.16961E+02...
            1.19345E+02   1.21629E+02   1.23237E+02   1.24856E+02   1.26716E+02...
            1.28145E+02   1.29702E+02   1.31181E+02   1.32359E+02   1.33507E+02...
            1.34440E+02   1.35603E+02   1.36720E+02   1.37552E+02   1.38449E+02...
            1.39005E+02   1.39607E+02   1.40207E+02   1.40765E+02   1.41337E+02...
            1.41857E+02   1.42446E+02   1.43042E+02   1.43639E+02   1.44135E+02...
            1.44605E+02   1.45120E+02   1.45740E+02   1.46304E+02   1.46874E+02...
            1.47389E+02   1.47898E+02   1.48849E+02   1.49737E+02   1.50485E+02...
            1.51322E+02   1.52500E+02   1.53486E+02   1.54619E+02   1.55886E+02...
            1.57269E+02];
        chers_phi=[9.33240E+01   9.22994E+01   9.12352E+01   9.04963E+01   8.94674E+01...
            8.87809E+01   8.79926E+01   8.70924E+01   8.65121E+01   8.60022E+01...
            8.54803E+01   8.50018E+01   8.46764E+01   8.43580E+01   8.40031E+01...
            8.37380E+01   8.34565E+01   8.31958E+01   8.29925E+01   8.27982E+01...
            8.26429E+01   8.24525E+01   8.22729E+01   8.21413E+01   8.20012E+01...
            8.19154E+01   8.18231E+01   8.17322E+01   8.16482E+01   8.15630E+01...
            8.14861E+01   8.13999E+01   8.13132E+01   8.12272E+01   8.11563E+01...
            8.10897E+01   8.10172E+01   8.09305E+01   8.08525E+01   8.07742E+01...
            8.07041E+01   8.06352E+01   8.05079E+01   8.03905E+01   8.02930E+01...
            8.01849E+01   8.00349E+01   7.99112E+01   7.97710E+01   7.96170E+01...
            7.94518E+01];
    end
    phi_spl = pchip(chers_r./100,chers_phi);
    mse_phi = ppval(phi_spl,mse_r);
    % plot
    %hold on;
    %polar(2*pi*mse_phi/360,mse_r,'+r');
    %polar(2*pi*chers_phi/360,chers_r./100,'ob');
    %hold off
    
end
    
fid_sh=[];
if lfull_shot
   fid_sh=fopen('fullshot.sh','w');
   fprintf(fid_sh,'#!/bin/tcsh\n');
end
% Output file
nfiles=length(timedex1);
for i=1:nfiles
    % Create File name
    tfile_name=[filename num2str(round(timedex1(i)*1000.),'%4.4d') '_' user];
    % Create Timewindow
    timewindow=timedex1(i):(timedex2(i)-timedex1(i))/ndt:timedex2(i);
    % Get the coil currents
    vmec_input.extcur=zeros(1,16);
    for j=1:length(coil_names)
        vmec_input.extcur(j) = sum(pchip(coil.([coil_names{j} '_time']),coil.(coil_names{j}),timewindow))./(ndt+1);
    end
    vmec_input.extcur(1) = vmec_input.extcur(1)*36./12.; % Note negative sign for correct TF field convention
    if luse_eddy
        for j=18:52
            vmec_input.extcur(j) = sum(pchip(gtime,eddy(j,:),timewindow))./(ndt+1);
        end
    end
    % Get the VMEC boundary guess
    if luse_efit
        for j=1:size(efit_r,1)
            rbdry(j)=sum(pchip(efit_gtime,efit_r(j,:),timewindow))./(ndt+1);
            zbdry(j)=sum(pchip(efit_gtime,efit_z(j,:),timewindow))./(ndt+1);
            rax=sum(pchip(efit_gtime,efit_rax,timewindow))./(ndt+1);
            zax=sum(pchip(efit_gtime,efit_zax,timewindow))./(ndt+1);
            pres =sum(pchip(efit_gtime,efit_pres,timewindow))./(ndt+1);
            pp   =sum(pchip(efit_gtime,efit_pp,timewindow))./(ndt+1);
            fp   =sum(pchip(efit_gtime,efit_ffp,timewindow))./(ndt+1);
            q    =sum(pchip(efit_gtime,efit_qpsi,timewindow))./(ndt+1);
            psi  =sum(pchip(efit_gtime,efit_psi,timewindow))./(ndt+1);
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
        jcurv = pp + fp;
        iota  = 1./q;
        q_spl = pchip(psi,q);
        fun   = @(x_fun) ppval(q_spl,x_fun);
        for k=1:nrho
            tflux(k) = integral(fun,psi(1),psi(k));
        end
        dex2 = find(pres<0,1,'first');
        pres(dex2) = 0.0;
        vmec_data.phiedge = tflux(dex2)*2*pi;
        for k=1:99
            vmec_input.am_aux_s(k) = pchip(tflux./vmec_data.phiedge,tflux,(k-1)/98);
            vmec_input.am_aux_f(k) = pchip(tflux./vmec_data.phiedge,pres,(k-1)/98);
        end
        for k=1:99
            vmec_input.am_aux_s(k) = pchip(tflux./tflux(nrho),tflux,(k-1)/98);
            vmec_input.ac_aux_f(k) = pchip(tflux./tflux(nrho),jcurv,(k-1)/98);
            vmec_input.ai_aux_s(k) = pchip(tflux./tflux(nrho),tflux,(k-1)/98);
            vmec_input.ai_aux_f(k) = pchip(tflux./tflux(nrho),iota,(k-1)/98);
        end
    end
    % Setup FLOW
    vmec_input.bcrit = 0.3;
    vmec_input.ph_type='akima_spline';
    vmec_input.ah_aux_s = flux_vals;
    vmec_input.ah_aux_f = 1.0*(1-flux_vals);
    vmec_input.pt_type='akima_spline';
    vmec_input.at_aux_s = flux_vals;
    vmec_input.at_aux_f = 1.0*(1-flux_vals);
    % Get Itor
    vmec_input.curtor=sum(pchip(ip_time,ip,timewindow))./(ndt+1);
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
    %write_namelist_int(fid,'MODE',2);
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
    %write_namelist_flt(fid,'DPHIEDGE_OPT',0.5);
    write_namelist_boo(fid,'LCURTOR_OPT',stel_input.lcurtor_opt);
    %write_namelist_flt(fid,'DCURTOR_OPT',0.01);
    write_namelist_boo(fid,'LPSCALE_OPT',stel_input.lpscale_opt);
    %write_namelist_flt(fid,'DPSCALE_OPT',0.5);
    fprintf(fid,'  LNE_F_OPT =');
    for j = 1: length(flux_vals)
        fprintf(fid,' T ');
    end
    fprintf(fid,'  DNE_F_OPT =');
    for j = 1: length(flux_vals)
        fprintf(fid,' 1.0 ');
    end
    fprintf(fid,'\n');
    fprintf(fid,'  LTE_F_OPT =');
    for j = 1: length(flux_vals)-1
        fprintf(fid,' T ');
    end
    fprintf(fid,' F ');
    fprintf(fid,'  DTE_F_OPT =');
    for j = 1: length(flux_vals)
        fprintf(fid,' 1.0 ');
    end
    fprintf(fid,'\n');
    fprintf(fid,'  LTI_F_OPT =');
    for j = 1: length(flux_vals)
        fprintf(fid,' T ');
    end
    fprintf(fid,'  DTI_F_OPT =');
    for j = 1: length(flux_vals)
        fprintf(fid,' 1.0 ');
    end
    fprintf(fid,'\n');
    if (lefield)
        fprintf(fid,'  LPHI_F_OPT =');
        for j = 1: length(flux_vals)
            fprintf(fid,' T ');
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    fprintf(fid,'!          PROFILE PARAMETERS\n');
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    if (luse_thomson)
        te_core=te_units.*sum(pchip(thom_time_core,thom_te_core',timewindow),2)./(ndt+1);
        ne_core=ne_units.*sum(pchip(thom_time_core,thom_ne_core',timewindow),2)./(ndt+1);
        stel_input.ne_aux_s=vmec_input.am_aux_s;
        stel_input.ne_aux_f=max(ne_core).*polyval([1 -1 0 0 -1 1],stel_input.ne_aux_s);
        dex=find(stel_input.ne_aux_s == max(stel_input.ne_aux_s(4:length(stel_input.ne_aux_s))),1,'first');
        stel_input.ne_aux_f(1:dex) = 0.75.*max(ne_core);
        fprintf(fid,'  NE_TYPE = ''two_power''\n');
        write_namelist_vec(fid,'NE_OPT',[0.75*max(ne_core)./1.0E18 4.0 1.0]);
        write_namelist_vec(fid,'NE_AUX_S',stel_input.ne_aux_s(1:dex));
        write_namelist_vec(fid,'NE_AUX_F',stel_input.ne_aux_f(1:dex)./ne_norm);
        stel_input.te_aux_s=vmec_input.am_aux_s;
        stel_input.te_aux_f=max(te_core).*polyval([1 -1 0 0 -1 1],stel_input.te_aux_s);
        dex=find(stel_input.te_aux_s == max(stel_input.te_aux_s(4:length(stel_input.te_aux_s))),1,'first');
        stel_input.te_aux_f(1:dex-1) = 0.75.*max(te_core);
        stel_input.te_aux_f(dex) = 0.0;
        fprintf(fid,'  TE_TYPE = ''two_power''\n');
        write_namelist_vec(fid,'TE_OPT',[0.75*max(te_core) 1.1 1.0]);
        write_namelist_vec(fid,'TE_AUX_S',stel_input.te_aux_s(1:dex));
        write_namelist_vec(fid,'TE_AUX_F',stel_input.te_aux_f(1:dex));
    end
    if luse_cer
        ti = sum(pchip(chers_time,chers_ti,timewindow)'./(ndt+1));
        stel_input.ti_aux_s=vmec_input.am_aux_s;
        stel_input.ti_aux_f=max(ti).*polyval([1 -1 0 0 -1 1],stel_input.ti_aux_s);
        dex=find(stel_input.ti_aux_s == max(stel_input.ti_aux_s(4:length(stel_input.ti_aux_s))),1,'first');
        stel_input.ti_aux_f(1:dex-1) = 0.75.*max(ti);
        stel_input.ti_aux_f(dex) = 0.0;
        fprintf(fid,'  TI_TYPE = ''two_power''\n');
        write_namelist_vec(fid,'TI_OPT',[0.75*max(ti) 1.1 1.0]);
        write_namelist_vec(fid,'TI_AUX_S',stel_input.ti_aux_s(1:dex));
        write_namelist_vec(fid,'TI_AUX_F',stel_input.ti_aux_f(1:dex));
    end
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    fprintf(fid,'!          EQUILIBRIUM AND GEOMETRY OPTIMIZATION PARAMETERS\n');
    fprintf(fid,'!-----------------------------------------------------------------------\n');
    if (lcurtor)
        write_namelist_flt(fid,'TARGET_CURTOR',vmec_input.curtor);
        write_namelist_flt(fid,'SIGMA_CURTOR',abs(std(pchip(ip_time./1000,ip,timewindow))));
    end
    if (lcoil_opt)
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          COIL OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        for j=1:17
            stel_input.sigma_extcur(j) = std(pchip(coil.([coil_names{j} '_time']),coil.(coil_names{j}),timewindow))./(ndt+1);
            if (stel_input.sigma_extcur(j) == 0.0)
                stel_input.sigma_extcur(j) = bigno;
            elseif (abs(stel_input.sigma_extcur(j))  < 0.01*abs(vmec_input.extcur(j)))
                stel_input.sigma_extcur(j) =  0.01*abs(vmec_input.extcur(j));
            end
            if (stel_input.sigma_extcur(j) < bigno)
                fprintf(fid,'  LEXTCUR_OPT(%2.2d) = T  DEXTCUR_OPT(%2.2d) = 1.0E-3  TARGET_EXTCUR(%2.2d) = %20.10f  SIGMA_EXTCUR(%2.2d) = %20.10e \n',...
                    j,j,j,vmec_input.extcur(j),j,stel_input.sigma_extcur(j));
            end
        end
        if luse_eddy
            for j=18:52
                stel_input.sigma_extcur(j) = std(pchip(gtime,eddy(j,:),timewindow))./(ndt+1);
                if (stel_input.sigma_extcur(j) == 0.0)
                    stel_input.sigma_extcur(j) = bigno;
                elseif (abs(stel_input.sigma_extcur(j))  < 0.01*abs(vmec_input.extcur(j)))
                    stel_input.sigma_extcur(j) =  0.01*abs(vmec_input.extcur(j));
                end
                fprintf(fid,'  LEXTCUR_OPT(%2.2d) = T  DEXTCUR_OPT(%2.2d) = 1.0  TARGET_EXTCUR(%2.2d) = %20.10f  SIGMA_EXTCUR(%2.2d) = %20.10e \n',...
                    j,j,j,vmec_input.extcur(j),j,stel_input.sigma_extcur(j));
            end
        end
    end
    if (luse_mags)
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          DIAGNO OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'  MAGDIAG_COIL = ''/u/slazerso/Sims/NSTX/coils/coils.nstx_SAL_020613''\n');
        if (luse_eddy)
            fprintf(fid,'  MAGDIAG_COIL = ''/u/slazerso/Sims/NSTX/coils/eddy/coils.nstx_060413_errf''\n');
        end
        for j=1:length(bprobe_names)
            bprobe_val=bprobe.(['B' upper(bprobe_names{j})]);
            bprobe_time=bprobe.(['B' upper(bprobe_names{j}) '_TIME']);
            bprobe_sigma=std(pchip(bprobe_time,bprobe_val,timewindow));
            bprobe_val=sum(pchip(bprobe_time,bprobe_val,timewindow))./(ndt+1);
            if bprobe_sigma < 0.0
                bprobe_sigma = 1.0E30;
                bprobe_val   = 1.0; % So optimizer doesn't get confused.
            else
                bprobe_sigma = 5.0E-2;
            end
            if bprobe_val == 0.0
                bprobe_val = 1.0;
                bprobe_sigma = 1.0E30;
            end
            if (j>=49 || j<=60)
                bprobe_sigma = bigno;
            end
            if abs(bprobe_val) < 1.0E-4
                bprobe_sigma = bigno;
            end
            fprintf(fid,'  TARGET_BPROBE(%2.2d) = %20.10E  SIGMA_BPROBE(%2.2d) = %20.10E\n',...
                    j,bprobe_val,j,bprobe_sigma);
        end
        for j=1:length(floop_names)
            floop_val=floop.(upper(floop_names{j}));
            floop_time=floop.([upper(floop_names{j}) '_TIME']);
            floop_sigma=std(pchip(floop_time,floop_val,timewindow));
            floop_val=sum(pchip(floop_time,floop_val,timewindow))./(ndt+1);
            if floop_sigma < 0.0
                floop_sigma = 1.0E30;
                floop_val   = 1.0;
            else
                floop_sigma = 2.0E-01;  % This seems to work well
            end
            if abs(floop_val) < 1.0E-4
                floop_sigma = bigno;
            end
            fprintf(fid,'  TARGET_FLUXLOOP(%2.2d) = %20.10E  SIGMA_FLUXLOOP(%2.2d) = %20.10E\n',...
                j,floop_val,j,floop_sigma);
        end
    end
    % Core Thomson
    if (luse_thomson)
        % Filter
        cutoff=80;
        te_core=te_units.*sum(pchip(thom_time_core,thom_te_core',timewindow),2)./(ndt+1);
        ne_core=ne_units.*sum(pchip(thom_time_core,thom_ne_core',timewindow),2)./(ndt+1);
        sigma_te_core=te_units.*std(pchip(thom_time_core,thom_te_core',timewindow),0,2);
        sigma_ne_core=ne_units.*std(pchip(thom_time_core,thom_ne_core',timewindow),0,2);
        te_core(te_core < cutoff) = 0.0;
        ne_core(te_core < cutoff) = 0.0;
        sigma_ne_core(ne_core == 0.0) = 1.0E30;
        sigma_te_core(te_core == 0.0) = 1.0E-4;
        for j=2:length(te_core)-1
            if (ne_core(j-1) == 0.0) && (ne_core(j+1) == 0.0)
                ne_core(j) = 0.0;
                sigma_ne_core(j) = 1.0E+30;
            end
            if (te_core(j-1) == 0.0) && (te_core(j+1) == 0.0)
                te_core(j) = 0.0;
                sigma_te_core(j) = 1.0E-4;
            end
        end
        % Output Ne
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          THOMSON (NE) PROFILE OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        stel_input.ne_aux_s=vmec_input.am_aux_s;
        stel_input.ne_aux_f=max(ne_core).*polyval([1 -1 0 0 -1 1],stel_input.ne_aux_s);
        dex=find(stel_input.ne_aux_s == max(stel_input.ne_aux_s(4:length(stel_input.ne_aux_s))),1,'first');
        stel_input.ne_aux_f(1:dex) = 0.75.*max(ne_core);
        write_namelist_vec(fid,'NE_AUX_S',stel_input.ne_aux_s(1:dex));
        write_namelist_vec(fid,'NE_AUX_F',stel_input.ne_aux_f(1:dex)./ne_norm);
        %fprintf(fid,'  NNE_PROF = %3.3d\n',length(thom_r_core)+length(thom_r_div)+length(thom_r_tan));
        for j=1:length(thom_r_core)
            if (sigma_ne_core(j) < 0.05*ne_core(j))
                sigma_ne_core(j) = 0.05*ne_core(j);
            end
            fprintf(fid,['  R_NE(%3.3d) = %20.10E  PHI_NE(%3.3d) = %20.10E' ...
                '  Z_NE(%3.3d) = %20.10E  TARGET_NE(%3.3d) = %20.10E'...
                '  SIGMA_NE(%3.3d) = %20.10E\n'],...
                j,thom_r_core(j),j,pi*thom_phi_core(j)/180,...
                j,thom_z_core(j),j,ne_core(j),...
                j,sigma_ne_core(j));
        end
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          THOMSON (TE) PROFILE OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        stel_input.te_aux_s=vmec_input.am_aux_s;
        stel_input.te_aux_f=max(te_core).*polyval([1 -1 0 0 -1 1],stel_input.te_aux_s);
        dex=find(stel_input.te_aux_s == max(stel_input.te_aux_s(4:length(stel_input.te_aux_s))),1,'first');
        stel_input.te_aux_f(1:dex-1) = 0.75.*max(te_core);
        stel_input.te_aux_f(dex) = 0.0;
        write_namelist_vec(fid,'TE_AUX_S',stel_input.te_aux_s(1:dex));
        write_namelist_vec(fid,'TE_AUX_F',stel_input.te_aux_f(1:dex));
        %fprintf(fid,'  NTE_PROF = %3.3d\n',length(thom_r_core)+length(thom_r_div)+length(thom_r_tan));
        % Output Te
        for j=1:length(thom_r_core)
            if (sigma_te_core(j) < 0.05*te_core(j))
                sigma_te_core(j) = 0.05*te_core(j);
            end
            fprintf(fid,['  R_TE(%3.3d) = %20.10E  PHI_TE(%3.3d) = %20.10E' ...
                '  Z_TE(%3.3d) = %20.10E  TARGET_TE(%3.3d) = %20.10E'...
                '  SIGMA_TE(%3.3d) = %20.10E\n'],...
                j,thom_r_core(j),j,pi*thom_phi_core(j)/180,...
                j,thom_z_core(j),j,te_core(j),...
                j,sigma_te_core(j));
        end
        j = 1;
        fprintf(fid,['  R_SEPARATRIX(1,1) = %20.10E  PHI_SEPARATRIX(1,1) = %20.10E' ...
            '  Z_SEPARATRIX(1,1) = %20.10E  TARGET_SEPARATRIX(1,1) = %20.10E'...
            '  SIGMA_SEPARATRIX(1,1) = %20.10E\n'],...
            thom_r_core(j),pi*thom_phi_core(j)/180,...
            thom_z_core(j),0.0,...
            abs(thom_r_core(j)-thom_r_core(j+1)));
        j = find(te_core == min(te_core),1,'first');
        fprintf(fid,['  R_SEPARATRIX(1,2) = %20.10E  PHI_SEPARATRIX(1,2) = %20.10E' ...
            '  Z_SEPARATRIX(1,2) = %20.10E  TARGET_SEPARATRIX(1,2) = %20.10E'...
            '  SIGMA_SEPARATRIX(1,2) = %20.10E\n'],...
            thom_r_core(j),pi*thom_phi_core(j)/180,...
            thom_z_core(j),0.0,...
            abs(thom_r_core(j)-thom_r_core(j-1)));
        
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
        dex=find(stel_input.ne_aux_s == max(stel_input.zeff_aux_s(4:length(stel_input.ne_aux_s))),1,'first');
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
                    n,sxr45r1s_r0(i),n,sxr45r1s_phi0(i),n,sxr45r1s_z0(i),...
                    n,sxr45r1s_r1(i),n,sxr45r1s_phi1(i),n,sxr45r1s_z1(i),...
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
                    n,sxr165r1s_r0(i),n,sxr165r1s_phi0(i),n,sxr165r1s_z0(i),...
                    n,sxr165r1s_r1(i),n,sxr165r1s_phi1(i),n,sxr165r1s_z1(i),...
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
                    n,sxr195r1s_r0(i),n,sxr195r1s_phi0(i),n,sxr195r1s_z0(i),...
                    n,sxr195r1s_r1(i),n,sxr195r1s_phi1(i),n,sxr195r1s_z1(i),...
                    n,line_sxr_target,n,line_sxr_sigma);
                n=n+1;
            end
        end
        % PA Chords
        for j=1:32
            line_sxr_target=sum(pchip(sxr90p1s_time./1000,sxr90p1s(j,:),timewindow))./(ndt+1);
            line_sxr_sigma=std(pchip(sxr90p1s_time./1000,sxr90p1s(j,:),timewindow));
            if line_sxr_target > 0
                fprintf(fid,['  R0_SXR(%3.3d) = %20.10E  PHI0_SXR(%3.3d) = %20.10E  Z0_SXR(%3.3d) = %20.10E'...
                    '  R1_SXR(%3.3d) = %20.10E  PHI1_SXR(%3.3d) = %20.10E  Z1_SXR(%3.3d) = %20.10E'...
                    '  TARGET_SXR(%3.3d) = %20.10E'...
                    '  SIGMA_SXR(%3.3d) = %20.10E\n'],...
                    n,sxr90p1s_r0(i),n,sxr90p1s_phi0(i),n,sxr90p1s_z0(i),...
                    n,sxr90p1s_r1(i),n,sxr90p1s_phi1(i),n,sxr90p1s_z1(i),...
                    n,line_sxr_target,n,line_sxr_sigma);
                n=n+1;
            end
            line_sxr_target=sum(pchip(sxr90m1s_time./1000,sxr90m1s(j,:),timewindow))./(ndt+1);
            line_sxr_sigma=std(pchip(sxr90m1s_time./1000,sxr90m1s(j,:),timewindow));
            if line_sxr_target > 0
                fprintf(fid,['  R0_SXR(%3.3d) = %20.10E  PHI0_SXR(%3.3d) = %20.10E  Z0_SXR(%3.3d) = %20.10E'...
                    '  R1_SXR(%3.3d) = %20.10E  PHI1_SXR(%3.3d) = %20.10E  Z1_SXR(%3.3d) = %20.10E'...
                    '  TARGET_SXR(%3.3d) = %20.10E'...
                    '  SIGMA_SXR(%3.3d) = %20.10E\n'],...
                    n,sxr90m1s_r0(i),n,sxr90m1s_phi0(i),n,sxr90m1s_z0(i),...
                    n,sxr90m1s_r1(i),n,sxr90m1s_phi1(i),n,sxr90m1s_z1(i),...
                    n,line_sxr_target,n,line_sxr_sigma);
                n=n+1;
            end
        end
    end
    if luse_cer
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          CHARGE EXCHANGE (TI) OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        ti = sum(pchip(chers_time,chers_ti,timewindow)'./(ndt+1));
        ti_err = sum(pchip(chers_time,chers_ti_err,timewindow)'./(ndt+1));
        ti_std = std(pchip(chers_time,chers_ti,timewindow)');
        stel_input.ti_aux_s=vmec_input.am_aux_s;
        stel_input.ti_aux_f=max(ti).*polyval([1 -1 0 0 -1 1],stel_input.ti_aux_s);
        dex=find(stel_input.ti_aux_s == max(stel_input.ti_aux_s(4:length(stel_input.ti_aux_s))),1,'first');
        write_namelist_vec(fid,'TI_AUX_S',stel_input.ti_aux_s(1:dex));
        write_namelist_vec(fid,'TI_AUX_F',stel_input.ti_aux_f(1:dex));
        for j=1:length(ti)
            ti_sig = max([ti_std(j) ti_err(j)]);
            if (ti_sig > ti(j))
                ti_sig = ti(j);
            end
            if (ti(j) == 0.0)
                ti_sig = 1.0E-4;
            elseif (abs(ti_sig/ti(j)) < 0.05)
                ti_sig = ti(j)*0.05;
            end
            fprintf(fid,['  R_TI(%3.3d) = %20.10E  PHI_TI(%3.3d) = %20.10E' ...
                '  Z_TI(%3.3d) = %20.10E  TARGET_TI(%3.3d) = %20.10E'...
                '  SIGMA_TI(%3.3d) = %20.10E\n'],...
                j,chers_r(j)./100,j,pi*chers_phi(j)./180,...
                j,chers_z(j),j,ti(j),...
                j,ti_sig);
        end
        j = find(ti == min(ti),1,'first');
        fprintf(fid,['  R_SEPARATRIX(1,3) = %20.10E  PHI_SEPARATRIX(1,3) = %20.10E' ...
            '  Z_SEPARATRIX(1,3) = %20.10E  TARGET_SEPARATRIX(1,3) = %20.10E'...
            '  SIGMA_SEPARATRIX(1,3) = %20.10E\n'],...
                chers_r(j)./100,pi*chers_phi(j)./180,...
                chers_z(j),0.0,...
                abs(chers_r(j)-chers_r(j-1))./1000);
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          CHARGE EXCHANGE (VT) OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        vt = sum(pchip(chers_time,chers_vt,timewindow)'./(ndt+1));
        vt_err = sum(pchip(chers_time,chers_vt_err,timewindow)'./(ndt+1));
        vt_std = std(pchip(chers_time,chers_vt,timewindow)');
        %stel_input.ti_aux_s=vmec_input.am_aux_s;
        %stel_input.ti_aux_f=max(ti).*polyval([1 -1 0 0 -1 1],stel_input.ti_aux_s);for j=1:length(ti)
        for j=1:length(vt)
            vt_sig = max([vt_std(j) vt_err(j)]);
            if (vt_sig > vt(j))
                vt_sig = vt(j);
            end
            if (vt(j) == 0.0)
                vt_sig = 1.0E-4;
            elseif (abs(vt_sig/vt(j)) < 0.05)
                vt_sig = vt(j)*0.05;
            end
            fprintf(fid,['  R_VPHI(%3.3d) = %20.10E  PHI_VPHI(%3.3d) = %20.10E' ...
                '  Z_VPHI(%3.3d) = %20.10E  TARGET_VPHI(%3.3d) = %20.10E'...
                '  SIGMA_VPHI(%3.3d) = %20.10E\n'],...
                j,chers_r(j)./100,j,pi*chers_phi(j)./180,...
                j,chers_z(j),j,vt(j),...
                j,vt_sig);
        end
    end
    if (luse_mse)
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        fprintf(fid,'!          MOTIONAL STARK EFFECT (MSE) OPTIMIZATION PARAMETERS\n');
        fprintf(fid,'!-----------------------------------------------------------------------\n');
        disp('  - Requesting MSE');
        coil_data=coil_biot_prep(coil_data);
        mse_phi = pi.*mse_phi./180;
        if (lefield)
            stel_input.phi_aux_s=vmec_input.am_aux_s;
            stel_input.phi_aux_f=1.0*stel_input.phi_aux_f;  % Set to constant =1.0 so normalization works
            dex=find(stel_input.phi_aux_s == max(stel_input.phi_aux_s(4:length(stel_input.phi_aux_s))),1,'first');
            write_namelist_vec(fid,'PHI_AUX_S',stel_input.phi_aux_s(1:dex));
            write_namelist_vec(fid,'PHI_AUX_F',stel_input.phi_aux_f(1:dex));
        end
        gamma = sum(pchip(mse_time,mse_pol,timewindow),2)./(ndt+1);
        sigma = std(pchip(mse_time,mse_pol,timewindow),0,2);
        dex = abs(sigma./gamma) < 0.05;
        sigma(dex) = sigma(dex).*0.0 + 0.05*gamma(dex);
        for j=1:length(mse_r)
            mse_pol_vac = 0.0;
            %x=mse_r(j)*cos(mse_phi(j));
            %y=mse_r(j)*sin(mse_phi(j));
            %z=mse_z(j);
            %[bx,by,bz]=coil_biot(coil_data,x,y,z,vmec_input.extcur);
            %br=bx*cos(mse_phi(j))+by*sin(mse_phi(j));
            %bphi=by*cos(mse_phi(j))-bx*sin(mse_phi(j));
            %mse_pol_vac=atan((mse_a1(j)*bz)/(mse_a2(j)*bphi+mse_a3(j)*br+mse_a4(j)*bz));
            % Write data;
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
            fprintf(fid,'  TARGET_MSE(%3d) = %20.10e ',j,pi.*gamma(j)./180);
            fprintf(fid,'  SIGMA_MSE(%3d) = %20.10e ',j,mse_factor.*pi.*sigma(j)./180);
            fprintf(fid,'  VAC_MSE(%3d) = %20.10e \n',j,mse_pol_vac);
        end
    end
    fprintf(fid,'/\n');
    % Handle the DIAGNO input namelist
    if (luse_mags)
        fprintf(fid,'&DIAGNO_IN\n');
        write_namelist_int(fid,'NU',72);
        write_namelist_int(fid,'NV',72);
        fprintf(fid,'  FLUX_DIAG_FILE = ''/u/slazerso/Sims/NSTX/probes/fluxloop_NSTX_n36.diagno''\n');
        fprintf(fid,'  BPROBES_FILE = ''/u/slazerso/Sims/NSTX/probes/bprobes_NSTX.diagno''\n');
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
    batchname=strtrim([num2str(round(timedex1(i)*1000.),'%4.4d') '_' user '.batch']);
    fid=fopen(batchname,'w+');
    fprintf(fid,'#!/bin/tcsh\n');  % NOTE TCSH
    fprintf(fid,'#Batchname\n');
    fprintf(fid,'#PBS -N %s\n',strtrim(['NSTX' num2str(shotnum,'%d') '_' num2str(round(timedex1(i)*1000.),'%4.4d') '_' user ]));
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
    fprintf(fid,'#PBS -l mem=%-dgb\n',2*numprocs);
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
    if i>1, fprintf(fid,'set ARGS_OLD="%s"\n',short_name_old); end
    fprintf(fid,'# --- echo the command syntax\n');
    fprintf(fid,'echo "The command syntax for this job is:"\n');
    fprintf(fid,'echo mpirun --mca btl ^openib -np $NPROCS $EXEPATH input.$ARGS\n');
    fprintf(fid,'echo " "\n');
    fprintf(fid,'cd $PBS_O_WORKDIR\n');
    fprintf(fid,'echo -n ''Started job at : '' ; date\n');
    fprintf(fid,'mkdir $ARGS\n');
    fprintf(fid,'cd $ARGS\n');
    fprintf(fid,'mv ../input.$ARGS . \n');
    if i>1, fprintf(fid,'cp ../$ARGS_OLD/wout_reset_file.nc . \n'); end
    if i==1 || ~lfull_shot
        fprintf(fid,'time mpirun --mca btl ^openib -np $NPROCS $EXEPATH input.$ARGS >& log.$ARGS\n');
    else
        fprintf(fid,'time mpirun --mca btl ^openib -np $NPROCS $EXEPATH input.$ARGS -restart >& log.$ARGS\n');
    end
    fprintf(fid,'rm *_opt* \n');
    fprintf(fid,' \n');
    fprintf(fid,'echo -n ''Ended job at  : '' ; date\n');
    fprintf(fid,'echo " " \n');
    fprintf(fid,'exit\n');
    fclose(fid);
    if (i==1 && lfull_shot)
        fprintf(fid_sh,'RUNID = $(qsub %s)\n',batchname);
        fprintf(fid_sh,'echo $RUNID\n');
    elseif (i>1 && lfull_shot)
        fprintf(fid_sh,'RUNID = $(qsub -W depend=afterany:$RUNID %s)\n',batchname);
        fprintf(fid_sh,'echo $RUNID\n');
    end
    short_name_old = strtrim(tfile_name(dex+1:length(tfile_name)));
    old_batchname=batchname;
end
if ~isempty(fid_sh)
    fclose(fid_sh);
end

% Close the Connection
mds_server.CloseTree(tree,shotnum);
mds_server.disconnect(server);

end

