function fig = beams3d_fidasimnamelist(varargin)
%BEAMS3D_FIDASIMNAMELIST Generates a BEAMS3D FIDASIM NAMELIST
%   The BEAMS3D_FIDASIMENAMELIST routine outputs the FIDASIM_INPUTS_B3D
%   namelist to a file passed to it.  In it's most general invocation it
%   will append the namelist to a file with default values.  By passing a
%   structure with a field called 'namelist' to the routine the structure
%   will be used to default values of the namelist. The optional argument
%   'plot' causes a plot to be made and the figure handle to be returned by
%   the code.
%
% Example usage
%
%       Basic
%       beams3d_fidasimenamelist('file','input.example'); % default values
%
%       User defined
%       f.namelist=1;
%       f.shot = 20181009043; % shot will be used
%       beams3d_fidasimenamelist('file','input.example',f);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

% Defaults
fid = [];
namelist = struct;
beam_data = [];
lplot=0;
fig = [];

% Handle inputs
if (nargin > 0)
    j=1;
    while j<=numel(varargin)
        if ischar(varargin{j})
            switch(varargin{j})
                case {'file','filename'}
                    j=j+1;
                    filename=varargin{j};
                    fid=fopen(filename,'a');
                case {'plot'}
                    lplot = 1;
            end
        elseif isstruct(varargin{j})
            if isfield(varargin{j},'namelist')
                namelist=varargin{j};
            elseif isfield(varargin{j},'beams3d')
                beam_data=varargin{j};
            end
        end
        j=j+1;
    end
end

% Hanled fileid
if isempty(fid)
    fprintf(fid,' ERROR- No file passed. ');
    return;
end
if fid == -1
    fprintf(fid,[' ERROR- Problem opening file: ' filename]);
    return;
end

% Default namelist
namelist0.shot = 3271980;
namelist0.time = 0.00;
namelist0.runid = 'none';
namelist0.result_dir = '%%FIDASIM_PATH%%';
namelist0.comment = 'Matlab beams3d_fidasimnamelist.m';
namelist0.device = 'BEAMS3D Run';
namelist0.tables_file = '%%FIDASIM_PATH%%/atomic_tables.h5';
namelist0.equilibrium_file = '';
namelist0.geometry_file = '';
namelist0.distribution_file = '';
namelist0.neutrals_file = '';
switch_names={'bes','dcx','halo','cold','brems','fida','npa','pfida',...
        'pnpa','neutron','birth','fida_wght','npa','wght'};
for i=1:length(switch_names)
    namelist0.(['calc_' switch_names{i}]) = 0;
end
namelist0.seed = -1;
namelist0.flr = 2;
namelist0.load_neutrals = 0;
namelist0.verbose = 1;
monte_names={'n_fida','n_npa','n_pfida','n_nbi','n_halo','n_dcx',...
    'n_birth'};
for i=1:length(monte_names)
    namelist0.(monte_names{i}) = 5000000;
end
namelist0.ai = 1;
namelist0.impurity_charge=5;
namelist0.nx = 128;
namelist0.ny = 128;
namelist0.nz = 128;
namelist0.xmin = -120.0;
namelist0.ymin = -120.0;
namelist0.zmin = -120.0;
namelist0.xmax =  120.0;
namelist0.ymax =  120.0;
namelist0.zmax =  120.0;
namelist0.alpha = 0.0;
namelist0.beta  = 0.0;
namelist0.gamma = 0.0;
namelist0.origin = [0,0,0];
namelist0.name = 'NBI';
namelist0.shape = 2;
namelist0.nbi_data_source = 'Matlab';
namelist0.src = [0,0,0];
namelist0.axis_nbi = [100,0,50];
namelist0.widy = 50.0;
namelist0.widz = 100.0;
namelist0.focy = 50.0;
namelist0.focz = 70.0;
namelist0.divy = [0,0,0];
namelist0.divz = [0,0,0];
namelist0.naperature = 1;
namelist0.ashape = 2;
namelist0.awidy = 2;
namelist0.awidz = 2;
namelist0.aoffy = 0;
namelist0.aoffz = 0;
namelist0.aoffz = 0;
namelist0.adist = 0;
namelist0.ab = 1;
namelist0.pinj = 1;
namelist0.einj = 55;
namelist0.current_fractions = [50 45 5];
namelist0.nchan = 1;
namelist0.system = 'SPEC';
namelist0.spec_data_source = 'Matlab';
namelist0.id = {'CHANNEL01'};
namelist0.radius = 100.;
namelist0.lens = zeros(3,1);
namelist0.axis_spec = ones(3,1);
namelist0.spot_size = 20;
namelist0.sigma_pi = 1.0;
namelist0.nlambda = 1024;
namelist0.lambdamin = 647.0;
namelist0.lambdamax = 669.0;
namelist0.ne_wght = 10;
namelist0.np_wght = 10;
namelist0.nphi_wght = 8;
namelist0.emax_wght = 55.0;
namelist0.nlambda_wght = 10;
namelist0.lambdamin_wght = 647.0;
namelist0.lambdamax_wght = 669.0;
namelist0.nchan_npa = 1;
namelist0.system_npa = 'NPA';
namelist0.npa_data_source = 'Matlab';
namelist0.id_npa = {'CHANNEL01'};
namelist0.radius_npa = 100.0;
namelist0.a_shape = 2;
namelist0.d_shape = 2;
namelist0.a_cent = zeros(3,1);
namelist0.d_cent = zeros(3,1);
namelist0.a_redge = ones(3,1);
namelist0.d_redge = ones(3,1);
namelist0.a_tedge = ones(3,1);
namelist0.d_tedge = ones(3,1);


% Default any missing values
names = fields(namelist0);
numfields = length(names);
for i = 1:numfields
    if ~isfield(namelist,names{i})
        namelist.(names{i})=namelist0.(names{i});
    end
end

% Check beam apperature consistency
n = namelist.naperature;
names = {'ashape','awidy','awidz','aoffy','aoffz','adist'};
numfields = length(names);
for i = 1:numfields
    if length(namelist.(names{i})) ~= n
        disp(['ERROR: Length of ' names{i} '='...
            num2str(length(namelist.(names{i})),' %i not naperature=')...
            num2str(n,'%i') ]);
    end
end

% Check spectrometer consistency
n = namelist.nchan;
names = {'id','radius','spot_size','sigma_pi'};
numfields = length(names);
for i = 1:numfields
    if length(namelist.(names{i})) ~= n
        disp(['ERROR: Length of ' names{i} '='...
            num2str(length(namelist.(names{i})),' %i not nchan=')...
            num2str(n,'%i') ]);
    end
end
if size(namelist.lens,2) ~= namelist.nchan
    disp(['ERROR: Size of lens = (3, ' ...
        num2str(size(namelist.lens,2),'%d') ') not nchan = '...
        num2str(namelist.nchan,'%d')]);
end
if size(namelist.axis_spec,2) ~= namelist.nchan
    disp(['ERROR: Size of axis_spec = (3, ' ...
        num2str(size(namelist.axis_spec,2),'%d') ') not nchan = '...
        num2str(namelist.nchan,'%d')]);
end

% Check npa consistency
n = namelist.nchan_npa;
names = {'id_npa','radius_npa','a_shape','d_shape'};
numfields = length(names);
for i = 1:numfields
    if length(namelist.(names{i})) ~= n
        disp(['ERROR: Length of ' names{i} '='...
            num2str(length(namelist.(names{i})),' %i not nchan_npa=')...
            num2str(n,'%i') ]);
    end
end
names = {'a_cent','d_cent','a_redge','d_redge','a_redge','d_redge'};
numfields = length(names);
for i = 1:numfields
    if size(namelist.(names{i}),2) ~= n
        disp(['ERROR: Size of ' names{i} '= (3,'...
            num2str(size(namelist.(names{i}),2),' %i not nchan_npa=')...
            num2str(n,'%i') ]);
    end
end

% Output values to screen
fprintf(fid, '&FIADASIM_INPUTS_B3D\n');
fprintf(fid, '!--------Run Information -------\n');
fprintf(fid,'  SHOT = %d\n',namelist.shot);
fprintf(fid,'  TIME = %5.2d\n',namelist.time);
fprintf(fid,'  RUNID = "%s"\n',namelist.runid);
fprintf(fid,'  COMMENT = "%s"\n',namelist.comment);
fprintf(fid,'  DEVICE = "%s"\n',namelist.device);
fprintf(fid, '!--------File Locations -------\n');
fprintf(fid,'  RESULTS_DIR = "%s"\n',namelist.result_dir);
fprintf(fid,'  TABLES_FILE = "%s"\n',namelist.tables_file);
fprintf(fid,'  EQUILIBRIUM_FILE = "%s"\n',namelist.equilibrium_file);
fprintf(fid,'  GEOMETRY_FILE = "%s"\n',namelist.geometry_file);
fprintf(fid,'  DISTRIBUTION_FILE = "%s"\n',namelist.distribution_file);
fprintf(fid,'  NEUTRALS_FILE = "%s"\n',namelist.neutrals_file);
fprintf(fid,'  LOAD_NEUTRALS = %d\n',namelist.load_neutrals);
fprintf(fid, '!--------Simulation Control -------\n');
for i=1:length(switch_names)
    temp = strcat('CALC_',upper(switch_names{i}));
    fprintf(fid,'  %s = %1d\n',temp,namelist.(['calc_' switch_names{i}]));
end
fprintf(fid, '!--------Debugging Switches -------\n');
fprintf(fid,'  SEED = %d\n',namelist.seed);
fprintf(fid,'  FLR = %d\n',namelist.flr);
fprintf(fid,'  VERBOSE = %d\n',namelist.verbose);
fprintf(fid, '!--------Monte Carlo Settings -------\n');
for i = 1:length(monte_names)
    fprintf(fid,'  %s = %d\n',upper(monte_names{i}),namelist.(monte_names{i}));
end
fprintf(fid, '!--------Wavelength Grid -------\n');
fprintf(fid,'  NLAMBDA = %d\n',namelist.nlambda);
fprintf(fid,'  LAMBDAMIN = %20.10E\n',namelist.lambdamin);
fprintf(fid,'  LAMBDAMAX = %20.10E\n',namelist.lambdamax);
fprintf(fid, '!--------Weight Function Settings -------\n');
fprintf(fid,'  NE_WGHT = %d\n',namelist.ne_wght);
fprintf(fid,'  NP_WGHT = %d\n',namelist.np_wght);
fprintf(fid,'  NPHI_WGHT = %d\n',namelist.nphi_wght);
fprintf(fid,'  NLAMBDA_WGHT = %d\n',namelist.nlambda_wght);
fprintf(fid,'  EMAX_WGHT = %20.10E\n',namelist.emax_wght);
fprintf(fid,'  LAMBDAMIN_WGHT = %20.10E\n',namelist.lambdamin_wght);
fprintf(fid,'  LAMBDAMAX_WGHT = %20.10E\n',namelist.lambdamax_wght);
fprintf(fid, '!--------Plasma Species -------\n');
fprintf(fid,'  AI = %d\n',namelist.ai);
fprintf(fid,'  IMPURITY_CHARGE = %d\n',namelist.impurity_charge);
fprintf(fid, '!--------Neutral Beam Deffinition -------\n');
fprintf(fid,'  NAME = "%s"\n',namelist.name);
fprintf(fid,'  AB = %20.10E\n',namelist.ab);
fprintf(fid,'  PINJ = %20.10E  EINJ = %20.10E\n',namelist.pinj, namelist.einj);
fprintf(fid,'  CURRENT_FRACTIONS = %20.10E %20.10E %20.10E\n',namelist.current_fractions);
fprintf(fid,'  SHAPE = %d\n',namelist.shape);
fprintf(fid,'  NBI_DATA_SOURCE = "%s"\n',namelist.nbi_data_source);
fprintf(fid,'  NX = %d  XMIN = %20.10E  XMAX = %20.10E \n',namelist.nx,namelist.xmin,namelist.xmax);
fprintf(fid,'  NY = %d  YMIN = %20.10E  YMAX = %20.10E \n',namelist.ny,namelist.ymin,namelist.ymax);
fprintf(fid,'  NZ = %d  ZMIN = %20.10E  ZMAX = %20.10E \n',namelist.nz,namelist.zmin,namelist.zmax);
fprintf(fid,'  ALPHA = %20.10E  BETA = %20.10E  GAMMA = %20.10E\n',namelist.alpha,namelist.beta,namelist.gamma);
fprintf(fid,'  ORIGIN = %20.10E %20.10E %20.10E\n',namelist.origin);
fprintf(fid,'  SRC = %20.10E %20.10E %20.10E\n',namelist.src);
fprintf(fid,'  AXIS_NBI = %20.10E %20.10E %20.10E\n',namelist.axis_nbi);
fprintf(fid,'  WIDY = %20.10E  WIDZ = %20.10E\n',namelist.widy,namelist.widz);
fprintf(fid,'  FOCY = %20.10E  FOCZ = %20.10E\n',namelist.focy,namelist.focz);
fprintf(fid,'  DIVY = %20.10E %20.10E %20.10E\n',namelist.divy);
fprintf(fid,'  DIVZ = %20.10E %20.10E %20.10E\n',namelist.divz);
fprintf(fid,'  NAPERATURE = %d\n',namelist.naperature);
for i = 1:namelist.naperature
    fprintf(fid,'    ASHAPE(%3.3d) = %d\n',i,namelist.ashape(i));
    fprintf(fid,'    ADIST(%3.3d) = %20.10E\n',i,namelist.adist(i));
    fprintf(fid,'    AWIDY(%3.3d) = %20.10E\n',i,namelist.awidy(i));
    fprintf(fid,'    AWIDZ(%3.3d) = %20.10E\n',i,namelist.awidz(i));
    fprintf(fid,'    AOFFY(%3.3d) = %20.10E\n',i,namelist.aoffy(i));
    fprintf(fid,'    AOFFZ(%3.3d) = %20.10E\n',i,namelist.aoffz(i));
end
fprintf(fid, '!--------Spectroscopic Deffinition -------\n');
fprintf(fid,'  SYSTEM = "%s"\n',namelist.system);
fprintf(fid,'  SPEC_DATA_SOURCE = "%s"\n',namelist.spec_data_source);
fprintf(fid,'  NCHAN = %d\n',namelist.nchan);
for i = 1:namelist.nchan
    fprintf(fid,'    ID(%3.3d) = "%s"\n',i,namelist.id{i});
    fprintf(fid,'    RADIUS(%3.3d)      = %20.10E\n',i,namelist.radius(i));
    fprintf(fid,'    SPOT_SIZE(%3.3d)   = %20.10E\n',i,namelist.spot_size(i));
    fprintf(fid,'    SIGMA_PI(%3.3d)    = %20.10E\n',i,namelist.sigma_pi(i));
    fprintf(fid,'    LENS(:,%3.3d)      = %20.10E  %20.10E  %20.10E\n',i,namelist.lens(:,i));
    fprintf(fid,'    AXIS_SPEC(:,%3.3d) = %20.10E  %20.10E  %20.10E\n',i,namelist.axis_spec(:,i));
end
fprintf(fid, '!--------NPA Deffinition (Not Yet Implemented -------\n');
fprintf(fid,'!  SYSTEM_NPA = "%s"\n',namelist.system_npa);
fprintf(fid,'!  NPA_DATA_SOURCE = "%s"\n',namelist.npa_data_source);
fprintf(fid,'!  NCHAN_NPA = %d\n',namelist.nchan_npa);
for i = 1:namelist.nchan_npa
    fprintf(fid,'!    ID_NPA(%3.3d) = "%s"\n',j,namelist.id_npa{i});
    fprintf(fid,'!    RADIUS_NPA(%3.3d)      = %20.10E\n',i,namelist.radius_npa(i));
    fprintf(fid,'!    A_SHAPE(%3.3d)   = %1d\n',i,namelist.a_shape(i));
    fprintf(fid,'!    D_SHAPE(%3.3d)    = %1d\n',i,namelist.d_shape(i));
    fprintf(fid,'!    A_CENT(:,%3.3d)      = %20.10E  %20.10E  %20.10E\n',i,namelist.a_cent(:,i));
    fprintf(fid,'!    A_REDGE(:,%3.3d)      = %20.10E  %20.10E  %20.10E\n',i,namelist.a_redge(:,i));
    fprintf(fid,'!    A_TEDGE(:,%3.3d)      = %20.10E  %20.10E  %20.10E\n',i,namelist.a_tedge(:,i));
    fprintf(fid,'!    D_CENT(:,%3.3d)      = %20.10E  %20.10E  %20.10E\n',i,namelist.d_cent(:,i));
    fprintf(fid,'!    D_REDGE(:,%3.3d)      = %20.10E  %20.10E  %20.10E\n',i,namelist.d_redge(:,i));
    fprintf(fid,'!    D_TEDGE(:,%3.3d)      = %20.10E  %20.10E  %20.10E\n',i,namelist.d_tedge(:,i));
end
fprintf(fid,'/\n');


if lplot
    fig = figure('Position',[1,1,1024,768],'Color','white','InvertHardCopy','off');
    % Plot box
    coord = [0 0 0; 1 0 0; 1 1 0; 0 1 0;...
        0 0 1; 1 0 1; 1 1 1; 0 1 1];
    idx = [1 2 3 4 1; 1 2 6 5 1; 2 3 7 6 2; 3 4 8 7 3; 4 1 5 8 4; 5 6 7 8 5]';
    x_0 = coord(:,1).*(namelist.xmax-namelist.xmin)+namelist.xmin;
    y_0 = coord(:,2).*(namelist.ymax-namelist.ymin)+namelist.ymin;
    z_0 = coord(:,3).*(namelist.zmax-namelist.zmin)+namelist.zmin;
    % need to rotate about angles https://d3denergetic.github.io/FIDASIM/page/02_physics/04_neutrals.html
    Rx = [1, 0, 0;
        0, cos(namelist.gamma), -sin(namelist.gamma);
        0, sin(namelist.gamma), cos(namelist.gamma)];
    Ry = [cos(namelist.beta), 0, sin(namelist.beta);
        0, 1, 0;
        -sin(namelist.beta), 0, cos(namelist.beta)];
    Rz = [cos(namelist.alpha), -sin(namelist.alpha), 0;
        sin(namelist.alpha), cos(namelist.alpha), 0;
        0, 0, 1];
    % Combine the rotation matrices
    rotation_matrix = Rz * Ry * Rx;
    % rotate
    points = (rotation_matrix * [x_0 y_0 z_0]')';
    x_plt = points(:,1) + namelist.origin(1);
    y_plt = points(:,2) + namelist.origin(2);
    z_plt = points(:,3) + namelist.origin(3);
    x_plt = x_plt./100; y_plt = y_plt./100; z_plt = z_plt./100; % cm to m
    patch(x_plt(idx),y_plt(idx),z_plt(idx),'r','FaceAlpha',0.33); hold on;
    % Plot source (in machine coordinates)
    x_plt = namelist.src(1);
    y_plt = namelist.src(2);
    z_plt = namelist.src(3);
    nx = namelist.axis_nbi;
    x_plt = x_plt./100; y_plt = y_plt./100; z_plt = z_plt./100; % cm to m
    quiver3(x_plt,y_plt,z_plt,nx(1),nx(2),nx(3),1,'k');
    n = norm(nx);
    nx = nx./n;
    nxz = cross(nx,[0 0 1]);
    nz = cross(nxz,nx);
    ny = cross(nz,nx);
    quiver3(x_plt,y_plt,z_plt,ny(1),ny(2),ny(3),n./200,'g');
    quiver3(x_plt,y_plt,z_plt,nz(1),nz(2),nz(3),n./200,'r');
    coord=[];
    if namelist.shape == 1 % Rectangle
        coord = [0 -1 -1; 0 1 -1; 0 1 1; 0 -1 1];
        idx = [1 2 3 4 1];
        idx_help = [1,2,3,4];
    else
        theta = deg2rad(0:5:360);
        coord(:,1) = zeros(length(theta),1);
        coord(:,2) = cos(theta);
        coord(:,3) = sin(theta);
        idx = 1:length(theta);
        idx_help = [1 19 39 55];
    end
    a = namelist.widy;
    b = namelist.widz;
    x_0 = coord(:,1);
    y_0 = coord(:,2).*a;
    z_0 = coord(:,3).*b;
    x_plt = x_0.*nx(1)+y_0.*ny(1)+z_0.*nz(1) + namelist.src(1);
    y_plt = x_0.*nx(2)+y_0.*ny(2)+z_0.*nz(2) + namelist.src(2);
    z_plt = x_0.*nx(3)+y_0.*ny(3)+z_0.*nz(3) + namelist.src(3);
    x_plt = x_plt./100; y_plt = y_plt./100; z_plt = z_plt./100; % cm to m
    patch(x_plt(idx),y_plt(idx),z_plt(idx),'b','FaceAlpha',0.1,'LineWidth',2);
    % Beam at multiple focy fractions
    x_line=x_plt(idx_help); y_line=y_plt(idx_help); z_line=z_plt(idx_help);
    frac_list = [1,2];
    for frac = frac_list
        frac2 = frac.*namelist.focy./namelist.focz;
        x_0 = coord(:,1);
        y_0 = coord(:,2).*a.*(1-frac);
        z_0 = coord(:,3).*b.*(1-frac2);
        x_plt = x_0.*nx(1)+y_0.*ny(1)+z_0.*nz(1)+ nx(1).*namelist.focy.*frac + namelist.src(1);
        y_plt = x_0.*nx(2)+y_0.*ny(2)+z_0.*nz(2)+ nx(2).*namelist.focy.*frac + namelist.src(2);
        z_plt = x_0.*nx(3)+y_0.*ny(3)+z_0.*nz(3)+ nx(3).*namelist.focy.*frac + namelist.src(3);
        x_plt = x_plt./100; y_plt = y_plt./100; z_plt = z_plt./100; % cm to m
        patch(x_plt(idx),y_plt(idx),z_plt(idx),'b','FaceColor','none');
        x_line = [x_line x_plt(idx_help)];
        y_line = [y_line y_plt(idx_help)];
        z_line = [z_line z_plt(idx_help)];
    end
    plot3(x_line',y_line',z_line','k');
    x_line = x_line(:,1); y_line = y_line(:,1); z_line = z_line(:,1);
    % Beam at multiple focy fractions
    frac_list = [1,2];
    for frac = frac_list
        frac2 = frac.*namelist.focz./namelist.focy;
        x_0 = coord(:,1);
        y_0 = coord(:,2).*a.*(1-frac2);
        z_0 = coord(:,3).*b.*(1-frac);
        x_plt = x_0.*nx(1)+y_0.*ny(1)+z_0.*nz(1)+ nx(1).*namelist.focz.*frac + namelist.src(1);
        y_plt = x_0.*nx(2)+y_0.*ny(2)+z_0.*nz(2)+ nx(2).*namelist.focz.*frac + namelist.src(2);
        z_plt = x_0.*nx(3)+y_0.*ny(3)+z_0.*nz(3)+ nx(3).*namelist.focz.*frac + namelist.src(3);
        x_plt = x_plt./100; y_plt = y_plt./100; z_plt = z_plt./100; % cm to m
        patch(x_plt(idx),y_plt(idx),z_plt(idx),'b','FaceColor','none');
        x_line = [x_line x_plt(idx_help)];
        y_line = [y_line y_plt(idx_help)];
        z_line = [z_line z_plt(idx_help)];
    end
    plot3(x_line',y_line',z_line','k');
    % Beam aperature
    for i = 1:namelist.naperature
        d = namelist.adist(i);
        a = namelist.awidy(i);
        b = namelist.awidz(i);
        if namelist.ashape == 1 % retangle
            u = [-a  a  a -a -a];
            v = [ b  b -b -b  b];
        else % ellipse
            theta = deg2rad(0:5:360);
            u = a.*cos(theta);
            v = b.*sin(theta);
        end
        x_plta = namelist.src(1) + d.*nx(1) + u.*ny(1) + v.*nz(1);
        y_plta = namelist.src(2) + d.*nx(2) + u.*ny(2) + v.*nz(2);
        z_plta = namelist.src(3) + d.*nx(3) + u.*ny(3) + v.*nz(3);
        x_plta = x_plta./100; y_plta = y_plta./100; z_plta = z_plta./100; % cm to m
        plot3(x_plta,y_plta,z_plta,'k');
    end
    % Now do spectrometers
    for i = 1:namelist.nchan
        s = namelist.spot_size(i); nx = namelist.axis_spec(:,i);
        r0 = namelist.axis_spec(:,i);
        n = norm(nx);
        nx = nx./n;
        nxz = cross(nx,[0 0 1]);
        nz = cross(nxz,nx);
        ny = cross(nz,nx);
        theta = linspace(0,2.*pi,18);
        z_0 = s.*sin(theta);
        y_0 = s.*cos(theta);
        x_0 = 0;
        x_plt = x_0.*nx(1)+y_0.*ny(1)+z_0.*nz(1) + namelist.lens(1,i);
        y_plt = x_0.*nx(2)+y_0.*ny(2)+z_0.*nz(2) + namelist.lens(2,i);
        z_plt = x_0.*nx(3)+y_0.*ny(3)+z_0.*nz(3) + namelist.lens(3,i);
        x_plt = [x_plt x_plt+nx(1).*namelist.radius(i)];
        y_plt = [y_plt y_plt+nx(2).*namelist.radius(i)];
        z_plt = [z_plt z_plt+nx(3).*namelist.radius(i)];
        x_plt = x_plt./100; y_plt = y_plt./100; z_plt = z_plt./100; % cm to m
        plot3(x_plt,y_plt,z_plt,'g');
    end

end


fclose(fid);


end