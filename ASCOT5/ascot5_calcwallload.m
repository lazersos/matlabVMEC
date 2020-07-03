function wall_load = ascot5_calcwallload(a5file,wallid,runid,varargin)
%ASCOT5_CALCWALLLOAD Calculate wall load array from ASCOT5 run
%   The ASCOT5_CALCWALLLOAD function returns an array of wall load data
%   based on an ASCOT5 runid and wallid. The code takes an ASCOT5
%   filename, wallid, and runid as inputs. If empty arrays are supplied
%   the active ID's are used. An optional argument 'hits' can be supplied
%   in which case the code returns an array of wall hits instead of the
%   wall load.  Wall load is returned in units of [W/m^2], assuming weight
%   is in units of [part/s].
%
%   Example:
%       a5file='ascot5_test.h5';
%       wallid=0838288192;
%       runid=0396210459;
%       wall_load = ascot5_calcwallload(a5file,wallid,runid); %Heat flux [W/m^2]
%       wall_load = ascot5_calcwallload(a5file,[],[]); %Active ID's
%       wall_load = ascot5_calcwallload(a5file,[],[],'hits'); % Strikes
%     
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0  

amu = 1.66053906660E-27;
wall_load = [];
lhits = 0;

% Handle varargin
if nargin > 3
    i=1;
    while i <= nargin-3
        switch varargin{i}
            case{'hits','nhits'}
                lhits = 1;
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end

% Check for file
if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

%Use active
if isempty(wallid)
    wallid=h5readatt(a5file,'/wall','active');
    disp(['  Using wallid: ' wallid]);
end
if isempty(runid)
    runid=h5readatt(a5file,'/results','active');
    disp(['  Using runid: ' runid]);
end

% Pull Wall
path_wall = ['/wall/wall_3D_' num2str(wallid,'%10.10i')];
try
    x1x2x3 = h5read(a5file,[path_wall '/x1x2x3']);
catch
    disp(['ERROR: Could not find wall: ' num2str(wallid,'%10.10i')]);
    return;
end

wall_load = zeros(1,size(x1x2x3,2));

% Pull particle data
path_run = ['/results/run_' num2str(runid,'%10.10i') '/endstate'];
try
    endcond = h5read(a5file,[path_run '/endcond']);
catch
    disp(['ERROR: Could not result: ' num2str(runid,'%10.10i')]);
    return;
end
walltile = h5read(a5file,[path_run '/walltile'])+1; % now in matlab index

% Correct walltile
dex = endcond ~= 8; % endcond=8 is wall hit
walltile(dex) = 0;
% Count hits
nhits=[];
mask = unique(walltile);
mask_final = mask(mask>0);
if ~lhits
    y1y2y3 = h5read(a5file,[path_wall '/y1y2y3']);
    z1z2z3 = h5read(a5file,[path_wall '/z1z2z3']);
    weight = h5read(a5file,[path_run '/weight']);
        mass = h5read(a5file,[path_run '/mass']).*amu; %in amu
    try
        vr = h5read(a5file,[path_run '/vr']);
        vphi = h5read(a5file,[path_run '/vphi']);
        vz = h5read(a5file,[path_run '/vz']);
    catch
        pr = h5read(a5file,[path_run '/prprt']);
        pphi = h5read(a5file,[path_run '/pphiprt']);
        pz = h5read(a5file,[path_run '/pzprt']);
        vr = pr./mass;
        vphi=pphi./mass;
        vz = pz./mass;
    end
    v2 = vr.*vr+vphi.*vphi+vz.*vz;
    q  = 0.5.*mass.*v2.*weight;
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
    wall_load(mask_final)=qflux;
    wall_load = wall_load./A;
else
    for i = mask'
        if (i==0), continue; end
        dex = walltile==i;
        nhits = [nhits; sum(dex)];
    end
    wall_load(mask_final) = nhits;
end

return;

end

