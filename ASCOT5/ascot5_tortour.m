function fig = ascot5_tortour(a5file,vmec_data,varargin)
%ASCOT5_TORTOUR makes a torus tour movie.
%   The ASCOT5_TORTOUR makes a movie of the wall using a VMEC equilibirium
%   for camera path information.  Various options are available the
%   'active' runid and wallid will be used unless passed to the function.
%   The 'movie' option causes the code to create a movie.
%
%   Example:
%       a5file='ascot5_W7X_20180821_012_5100_h8.h5';
%       vmec_data=read_vmec('wout_W7X_20180821_012_5100_H8.nc');
%       wallid=0838288192;
%       runid=0396210459;
%       fig = fig=ascot5_tortour(ascot_file,vmec_data,'log','runid',id{1});
%
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

nfilter=0;
lhits=0;
llog=0;
lmovie=0;
amu = 1.66053906660E-27;
ec = 1.60217662E-19;
wallid=[];
runid=[];


% Handle varargin
if nargin > 2
    i=1;
    while i <= nargin-2
        switch varargin{i}
            case{'wallid'}
                i=i+1;
                wallid=num2str(varargin{i},'%10.10i');
            case{'runid'}
                i=i+1;
                runid=num2str(varargin{i},'%10.10i');
            case{'movie'}
                lmovie=1;
            case{'log'}
                llog=1;
            case{'hits','nhits'}
                lhits = 1;
            case{'nfilter'}
                i=i+1;
                nfilter=varargin{i};
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

% Setup movie name
movname = a5file(1:end-3);

% Use active if not provided
if isempty(wallid)
    wallid=h5readatt(a5file,'/wall','active');
    disp(['  Using wallid: ' wallid]);
end
if isempty(runid)
    try
        runid=h5readatt(a5file,'/results','active');
        disp(['  Using runid: ' runid]);
    catch
        disp('Error: Cannot find active run!');
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


wall_strikes = zeros(1,size(x1x2x3,2));
wall_load = zeros(1,size(x1x2x3,2));
if ~isempty(runid)
    % Pull data
    path = ['/results/run_' num2str(runid,'%10.10i') '/endstate'];
    try
        endcond = h5read(a5file,[path '/endcond']);
    catch
        disp(['ERROR: Could not result: ' num2str(runid,'%10.10i')]);
        return;
    end
    walltile = h5read(a5file,[path '/walltile'])+1;
    % Correct walltile
    dex = endcond ~= 8; % endcond=8 is wall hit
    walltile(dex) = 0;
    % Count hits
    nhits=[];
    mask = unique(walltile);
    for i = mask'
        if (i==0), continue; end
        dex = walltile==i;
        nhits = [nhits; sum(dex)];
    end
    mask2 = mask(mask>0);
    wall_strikes(mask2) = nhits;
    % Filter to zero
    dex = wall_strikes<=nfilter;
    wall_strikes(dex) = 0;
    % Calc heatload if desired
    if ~lhits
        weight = h5read(a5file,[path '/weight']);
        vr = h5read(a5file,[path '/vr']);
        vphi = h5read(a5file,[path '/vphi']);
        vz = h5read(a5file,[path '/vz']);
        mass = h5read(a5file,[path '/mass']).*amu; %in amu
        v2 = vr.*vr+vphi.*vphi+vz.*vz;
        q  = 0.5.*mass.*v2.*weight;
        %vpar = h5read(a5file,[path '/vpar']);
        %mu = h5read(a5file,[path '/mu']).*ec; %in eV/T
        %br = h5read(a5file,[path '/br']);
        %bphi = h5read(a5file,[path '/bphi']);
        %bz = h5read(a5file,[path '/bz']);
        %modb = sqrt(br.*br+bphi.*bphi+bz.*bz);
        %energy = mu.*modb;
        %q = (energy+0.5.*mass.*vpar.*vpar).*weight;
    
        % Calc Area
        V0=[x1x2x3(2,:)-x1x2x3(1,:);y1y2y3(2,:)-y1y2y3(1,:);z1z2z3(2,:)-z1z2z3(1,:)];
        V1=[x1x2x3(3,:)-x1x2x3(1,:);y1y2y3(3,:)-y1y2y3(1,:);z1z2z3(3,:)-z1z2z3(1,:)];
        F = cross(V0,V1);
        A = 0.5.*sqrt(sum(F.*F));
        
        % Convert to heatflux
        qflux=[];
        for i = mask'
            if (i==0), continue; end
            dex = walltile==i;
            qtemp = q(dex);
            qflux = [qflux; sum(qtemp)];
        end
        wall_load(mask2)=qflux;
        wall_load = wall_load./A;
        % Always filter
        dex = wall_strikes==0;
        wall_load(dex) = 0;
    end
end
if llog
    if lhits
        wall_strikes = log10(wall_strikes);
    else
        wall_load = log10(wall_load);
    end
    disp('Log10 color contours');
end

% Make plot
win_res=[1 1 1920 1080];         % Movie resolution
bgcolor='w';   % background color
fig=figure('DoubleBuffer','on','Position',win_res,...
        'Color',bgcolor,'BackingStore','on','MenuBar','none',...
        'Name','Fieldlines','InvertHardcopy','off');
set(gca,'nextplot','replacechildren','XColor',bgcolor,...
        'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
if isempty(wall_load)
    patch('XData',x1x2x3,'YData',y1y2y3,'ZData',z1z2z3,'EdgeColor','none');
else
    if lhits
        patch(x1x2x3,y1y2y3,z1z2z3,wall_strikes,'EdgeColor','none');
    else
        patch(x1x2x3,y1y2y3,z1z2z3,wall_load,'EdgeColor','none');
    end
end

% Create the equilibrium hoops
theta=0:2*pi./255:2*pi;
zeta=0:2*pi./(vmec_data.nfp*8):2*pi;
req=cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
zeq=sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
for i = 1:size(req,3)
    xeq(:,:,i) = req(:,:,i).*cos(zeta(i));
    yeq(:,:,i) = req(:,:,i).*sin(zeta(i));
end
ns=vmec_data.ns;
hold on;
plot3(squeeze(xeq(ns,:,:)),squeeze(yeq(ns,:,:)),squeeze(zeq(ns,:,:)),'r')


% Now set camera position
fps=30;
time=1*60;
nstep = time.*fps;
zeta=0:2*pi./nstep:2*pi;
rax=cfunct(0,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
zax=sfunct(0,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
rax = squeeze(rax(1,:,:));
zax = squeeze(zax(1,:,:));
xax = rax.*cos(zeta');
yax = rax.*sin(zeta');
camproj('perspective');
camva(60);
axis equal;
colormap jet;
if (llog)
    caxis([2 6]);
else
    caxis([0 1E6]);
end

% Now make movie
i=1;
campos([xax(i) yax(i) zax(i)]);
camtarget([xax(i+1) yax(i+1) zax(i+1)]);
ha=camlight(0,0);
drawnow;
if (lmovie)
    mov=VideoWriter([movname '_movie.avi']);
    mov.Quality= 100;
    mov.FrameRate=30;
    open(mov);
    frame=getframe(fig);
    writeVideo(mov,frame);
end
delete(ha);
for i=2:nstep
    campos([xax(i) yax(i) zax(i)]);
    camtarget([xax(i+1) yax(i+1) zax(i+1)]);
    ha=camlight(0,0);
    drawnow;
    if lmovie
        frame=getframe(fig);
        writeVideo(mov,frame);
    end
    delete(ha);
end


if lmovie
    close(mov);
end

return;

end

