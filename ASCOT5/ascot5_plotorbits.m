function fig = ascot5_plotorbits(a5file,runid,varargin)
%ASCOT5_PLOTORBITS plots the orbit data from an ASCOT5 file
%   The ASCOT5_PLOTORBITS routine plots the 'orbits' field of an ASCOT5
%   HDF5 file. It takes the ASCOT5 filename, and runid as input. If an
%   empty array is passed then the active id for orbits is used. Optional
%   arguments control what type of plot is generated.
%  
%   Options arguments:
%       'xyz':      XYZ orbit plot (default)
%       'flux':     Rho-theta polar plot
%       'pll':      Parallel-Perpendicular Momentum plot
%
%   Example:
%       a5file='ascot5_test.h5';
%       orbitid=0396210459;
%       fig=ascot5_plotorbits(a5file,orbitid);
%       fig=ascot5_plotorbits(a5file,[]); %Use active ID
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

% Defaults
amu = 1.66053906660E-27;
plottype='xyz';

% Handle varargin
if nargin > 2
    i=1;
    while i <= nargin-2
        switch varargin{i}
            case{'xyz','flux','vll','pll'}
                plottype=varargin{i};
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

% Use active run
if isempty(runid)
    try
        runid=h5readatt(a5file,'/results','active');
        disp(['  Using runid: ' runid]);
    catch
        runid=[];
        return;
    end
end

% Set run path
path = ['/results/run_' num2str(runid,'%10.10i') '/orbit'];

% Get IDS and time for sorting
ids = double(h5read(a5file,[path '/ids']));
time = double(h5read(a5file,[path '/mileage']));

% get array sizes
npart = length(unique(ids));
nsteps = round(length(ids)./max(ids));
time = reshape(time,[nsteps npart]);

% create particle sorting array
[time,idex] = sort(time);

fig = figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
switch plottype
    case 'xyz'
        % Extract position
        r = h5read(a5file,[path '/r']);
        z = h5read(a5file,[path '/z']);
        phi = h5read(a5file,[path '/phi']);
        x = r.*cos(phi);
        y = r.*sin(phi);
        
        % Reshape by particle index
        x = reshape(x,[nsteps npart]);
        y = reshape(y,[nsteps npart]);
        z = reshape(z,[nsteps npart]);
        
        % Reorder in time
        for i=1:npart
            x(:,i) = x(idex(:,i),i);
            y(:,i) = y(idex(:,i),i);
            z(:,i) = z(idex(:,i),i);
        end
        
        % Make plot
        plot3(x,y,z);
        set(gca,'FontSize',24);
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Z [m]');
        title('Particle Orbits');
    case 'flux'
        % Extract position
        rho = h5read(a5file,[path '/rho']);
        theta = h5read(a5file,[path '/theta']);
        
        % Reshape by particle index
        rho = reshape(rho,[nsteps npart]);
        theta = reshape(theta,[nsteps npart]);
        
        % Reorder in time
        for i=1:npart
            rho(:,i) = rho(idex(:,i),i);
            theta(:,i) = theta(idex(:,i),i);
        end
        
        % Make plot
        polarplot(theta,rho);
        set(gca,'FontSize',24);
        title('Particle Orbits');
        
    case {'vll','pll'}
        % Extract position
        pr = h5read(a5file,[path '/pr']);
        pz = h5read(a5file,[path '/pz']);
        pphi = h5read(a5file,[path '/pphi']);
        br = h5read(a5file,[path '/br']);
        bz = h5read(a5file,[path '/bz']);
        bphi = h5read(a5file,[path '/bphi']);
        
        % get pll
        b=sqrt(br.*br+bphi.*bphi+bz.*bz);
        pll = (pr.*br+pphi.*bphi+pz.*bz)./b;
        p2  = pr.*pr+pphi.*pphi+pz.*pz;
        pperp = sqrt(p2-pll.*pll);
        
        % Reshape by particle index
        pll = reshape(pll,[nsteps npart]);
        pperp = reshape(pperp,[nsteps npart]);
        
        % Reorder in time
        for i=1:npart
            pll(:,i) = pll(idex(:,i),i);
            pperp(:,i) = pperp(idex(:,i),i);
        end
        
        % Make plot
        plot(pll./amu,pperp./amu);
        set(gca,'FontSize',24);
        xlabel('Parallel Momentum [u*m/s]');
        ylabel('Perpendicular Momentum [u*m/s]')
        title('Phase space');
end


end

