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
%       a5file='ascot5_W7X_20180821_012_5100_h8.h5';
%       wallid=0838288192;
%       runid=0396210459;
%       fig=ascot5_plotwall(a5file,wallid,[]); % Just wall
%       fig=ascot5_plotwall(a5file,wallid,runid); %Heat flux [W/m^2]
%       fig=ascot5_plotwall(a5file,wallid,runid,'hits'); %Wall strikes
%       fig=ascot5_plotwall(a5file,wallid,runid,'nfilter',5); %Filtering
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

nfilter=0;
lhits=0;
llog=0;
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
    runid=h5readatt(a5file,'/results','active');
    disp(['  Using runid: ' runid]);
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
        walltile = h5read(a5file,[path '/walltile']);
    catch
        disp(['ERROR: Could not result: ' num2str(runid,'%10.10i')]);
        return;
    end
    % Count hits
    nhits=[];
    mask = unique(walltile)+1;
    for i = mask'
        if (i==0), continue; end
        dex = walltile==i;
        nhits = [nhits; sum(dex)];
    end
    wall_strikes(mask) = nhits;
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
        wall_load(mask)=qflux;
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
fig=figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
if isempty(wall_load)
    patch('XData',x1x2x3,'YData',y1y2y3,'ZData',z1z2z3,'EdgeColor','none');
else
    if lhits
        patch(x1x2x3,y1y2y3,z1z2z3,wall_strikes,'EdgeColor','none');
    else
        patch(x1x2x3,y1y2y3,z1z2z3,wall_load,'EdgeColor','none');
    end
end

return;

end

