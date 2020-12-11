function fig = ascot5_plotwall(a5file,wallid,runid,varargin)
%ASCOT5_PLOTWALL plots the wall structure of an ASCOT run.
%   The ASCOT5_PLOTWALL function plots the wall of an ASCOT5 run.  If a
%   runid is provided the code will use the result endstate to calculate
%   the heat flux to the wall.  Passing the optional argument 'hits' will
%   cause the code to plot the number of hits per cell rather than the heat
%   flux.  The option 'nfilter' will filter the results so any cell with
%   fewer than this number of hits is zeroed out.
%
%   Example:
%       a5file='ascot5_test.h5';
%       wallid=0838288192;
%       runid=0396210459;
%       fig=ascot5_plotwall(a5file,wallid,[]); % Just wall
%       fig=ascot5_plotwall(a5file,wallid,runid); %Heat flux [W/m^2]
%       fig=ascot5_plotwall(a5file,wallid,runid,'hits'); %Wall strikes
%       fig=ascot5_plotwall(a5file,wallid,runid,'nfilter',5); %Filtering
%       fig=ascot5_plotwall(a5file,wallid,runid,'log'); %LOG10 scale
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.2

nfilter=0;
lhits=0;
llog=0;
lpts=0;
amu = 1.66053906660E-27;
ec = 1.60217662E-19;


% Handle varargin
if nargin > 3
    i=1;
    while i <= nargin-3
        switch varargin{i}
            case{'log'}
                llog=1;
            case{'hits','nhits'}
                lhits = 1;
            case{'nfilter'}
                i=i+1;
                nfilter=varargin{i};
            case{'points'}
                lpts=1;
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end

fig=[];
% Check for file
if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

% Use active
if isempty(wallid)
    wallid=h5readatt(a5file,'/wall','active');
    disp(['  Using wallid: ' wallid]);
end
if isempty(runid)
    try
        runid=h5readatt(a5file,'/results','active');
        disp(['  Using runid: ' runid]);
    catch
        runid=[];
        disp('  No runid found, displaying wall only.');
    end
end

% Pull Wall
path = ['/wall/wall_3D_' num2str(wallid,'%10.10i')];
try
    x1x2x3 = h5read(a5file,[path '/x1x2x3']);
catch
    disp(['ERROR: Could not find wall: ' num2str(wallid,'%10.10i')]);
    return;
end
y1y2y3 = h5read(a5file,[path '/y1y2y3']);
z1z2z3 = h5read(a5file,[path '/z1z2z3']);

wall_load = [];
wall_strikes = [];

if ~isempty(runid)
    if (lhits)
        wall_strikes = ascot5_calcwallload(a5file,wallid,runid,'hits');
    else
        wall_load = ascot5_calcwallload(a5file,wallid,runid);
        if nfilter > 0
            wall_strikes = ascot5_calcwallload(a5file,wallid,runid,'hits');
            dex = wall_strikes <= nfilter;
            wall_load(dex) = 0;
        end
    end
    if lpts
        path = ['/results/run_' num2str(runid,'%10.10i') '/endstate'];
        r_pts = h5read(a5file,[path '/r']);
        p_pts = h5read(a5file,[path '/phi']);
        z_pts = h5read(a5file,[path '/z']);
        x_pts = r_pts.*cosd(p_pts);
        y_pts = r_pts.*sind(p_pts);
    end
end
    
% Handle Log plots
if llog
    if lhits && ~isempty(wall_strikes)
        wall_strikes = log10(wall_strikes);
    elseif ~isempty(wall_load)
        wall_load = log10(wall_load);
    end
    disp('Log10 color contours');
end

% Make plot
fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
if isempty(wall_load) && isempty(wall_strikes)
    hp=patch('XData',x1x2x3,'YData',y1y2y3,'ZData',z1z2z3,'EdgeColor','none','FaceColor',[1 1 1].*0.25);
%set(hp,'AmbientStrength',1.0,'SpecularStrength',0.75,'DiffuseStrength',0.9,'SpecularColorReflectance',0.5);
else
    if lhits
        hp=patch(x1x2x3,y1y2y3,z1z2z3,wall_strikes,'EdgeColor','none');
    else
        hp=patch(x1x2x3,y1y2y3,z1z2z3,wall_load,'EdgeColor','none');
    end
    %set(hp,'AmbientStrength',1.0,'SpecularStrength',0.75,'DiffuseStrength',1);
    set(hp,'AmbientStrength',1.0,'SpecularStrength',0.75,'DiffuseStrength',0.9,'SpecularColorReflectance',0);
end
if lpts
    hold on;
    plot3(x_pts,y_pts,z_pts,'.r','MarkerSize',0.1);
    hold on;
end
set(gca,'Color','blue');

%cmap
%r0 = [138; 138; 138]./255.0;
r0 = [96; 96; 96]./255.0;
r1 = [255; 0; 0]./255.0;
r2 = [255; 255; 0]./255.0;
dr = (r1-r0);
cmap = r0+dr*(0:1.0/127:1);
dr = (r2-r1);
cmap = [cmap r1+dr*(0:1.0/127:1)];
colormap(cmap');

return;

end

