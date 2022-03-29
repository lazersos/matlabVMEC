function fig = ascot5_plotorbits(a5file,runid,varargin)
%ASCOT5_PLOTORBITS plots the orbit data from an ASCOT5 file
%   The ASCOT5_PLOTORBITS routine plots the 'orbits' field of an ASCOT5
%   HDF5 file. It takes the ASCOT5 filename, and runid as input. If an
%   empty array is passed then the active id for orbits is used. Optional
%   arguments control what type of plot is generated.
%
%   Options arguments:
%       'parts':    Select subset of Markers by ID
%       'xyz':      XYZ orbit plot (default)
%       'flux':     Rho-theta polar plot
%       'fluxperp': Rho-theta xy-plot
%       'pll':      Parallel-Perpendicular Momentum plot
%       'rz':       Plot in R-Z Plane
%       'movie':    Make a movie of the plot
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
parts = [];
lmovie = 0;
pspace = 0;
filename = 'default';
startframe = 1; %Number of Start-Frame - ntail
endframe = 20000; %Number of last frame
skipframes = 15; %How many data points to skip between images
ntail = 20;

% Handle varargin
if nargin > 2
    i=1;
    while i <= nargin-2
        switch varargin{i}
            case{'xyz','flux','vll','pll', 'rz', 'fluxperp'}
                plottype=varargin{i};
            case 'parts'
                i = i+1;
                parts=varargin{i};
            case 'movie'
                i = i + 1;
                filename = varargin{i};
                lmovie = 1;
            case 'ntail'
                i = i+1;
                ntail = varargin{i};
            case {'pspace', 'vspace'}
                    pspace = 1;
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
if isempty(parts)
    parts = unique(ids);
end


% get array sizes

nsteps = round(length(ids)./max(ids));
endframe = nsteps;
dex = ismember(ids,parts);
ids = ids(dex);
time = time(dex);

npart = length(unique(parts));

time = reshape(time,[nsteps npart]);

% create particle sorting array
[time,idex] = sort(time);

fig = figure('Position',get(0,'ScreenSize'));%figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
ax1 = gca;
normh = 0.550;
if pspace
    ax1 = axes('Position', [0.03 0.3 0.4 normh]);
    ax2 = axes('Position', [-0.16 0.3 0.4 normh]);
end
v = VideoWriter([filename '.mp4'], 'MPEG-4');
open(v);
endframe = 

switch plottype
    case 'xyz'
        % Extract position
        
        
        modphi = mod(phi,360);
        edges = linspace(0,360,360);
        discphi = discretize(modphi,edges);
        colors = parula(360);
        % Make plot
        %plot3(ax1,x,y,z);
        hold on;
        if lmovie
            
            for i=startframe+ntail:skipframes:endframe
                
                %s1 = scatter3(ax1,x(i-ntail:i,:),y(i-ntail:i,:),z(i-ntail:i,:));
                s2 = plot3(ax1,x(i-ntail:i,:),y(i-ntail:i,:),z(i-ntail:i,:),'Color', colors(discphi(i),:));
                campos([95.6559   -2.0102   32.5105]);
                camtarget([-1.0351   -0.2400    0.1142]);
                camup([0 0 1]);
                set(gca,'FontSize',24);
                xlabel('X [m]');
                ylabel('Y [m]');
                zlabel('Z [m]');
                axis equal
                caxis(ax1,[0 360]);
                zlim(ax1,[-1.25 1.25]);
                ylim(ax1,[-6 6]);
                xlim(ax1,[-6 6]);
                
                title('Particle Orbits');
                
                frame = getframe(fig); %for writing
                writeVideo(v,frame);
                %delete(s1);
                delete(s2);
            end
        else
            scatter3(ax1,x,y,z,0.1,time);
            set(gca,'FontSize',24);
            xlabel('X [m]');
            ylabel('Y [m]');
            zlabel('Z [m]');
            title('Particle Orbits');
            axis equal;
        end
    case 'flux'
        % Extract position
        rho = h5read(a5file,[path '/rho']);
        theta = h5read(a5file,[path '/theta']);
        phi = h5read(a5file,[path '/phi']);
        rho = rho(dex);
        theta = theta(dex);
        phi = phi(dex);
        
        % Reshape by particle index
        rho = reshape(rho,[nsteps npart]);
        theta = reshape(theta,[nsteps npart]);
        phi = reshape(phi,[nsteps npart]);
        
        % Reorder in time
        for i=1:npart
            rho(:,i) = rho(idex(:,i),i);
            theta(:,i) = theta(idex(:,i),i);
            phi(:,i) = phi(idex(:,i),i);
        end
        
        % Make plot
        delete(ax1);
        ax1 = polaraxes;
        modphi = mod(phi,180);
        edges = linspace(0,180,360);
        discphi = discretize(modphi,edges);
        colors = parula(360);
        
        if lmovie
            for i=startframe+ntail:skipframes:endframe
                polarplot(ax1,deg2rad(theta(i-ntail:i,:)),rho(i-ntail:i,:), 'Color', colors(discphi(i),:));
                ax1.RLim = [0 1.5];
                
                %scatter(ax1, deg2rad(theta(i,:),rho(i,:));
                set(gca,'FontSize',24);
                title('Particle Orbits');
                frame = getframe(fig); %for writing
                writeVideo(v,frame);
            end
        else
            polarplot(ax1,deg2rad(theta),rho);
            set(gca,'FontSize',24);
            title('Particle Orbits');
        end
    case 'fluxperp'
        % Extract position
        rho = h5read(a5file,[path '/rho']);
        theta = h5read(a5file,[path '/theta']);
        phi = h5read(a5file,[path '/phi']);
        rho = rho(dex);
        theta = theta(dex);
        phi = phi(dex);
        
        % Reshape by particle index
        rho = reshape(rho,[nsteps npart]);
        theta = reshape(theta,[nsteps npart]);
        phi = reshape(phi,[nsteps npart]);
        
        % Reorder in time
        for i=1:npart
            rho(:,i) = rho(idex(:,i),i);
            theta(:,i) = theta(idex(:,i),i);
            phi(:,i) = phi(idex(:,i),i);
        end
        
        % Make plot
        modphi = mod(phi,180);
        edges = linspace(0,180,360);
        discphi = discretize(modphi,edges);
        colors = parula(360);
        poloidalextent = 360.0;
        if lmovie
            revolutions = floor(theta./ poloidalextent);
            for i=startframe+ntail:skipframes:endframe
                plotx = theta(i-ntail:i,:)- poloidalextent*revolutions(i-ntail,:);
                ploty = rho(i-ntail:i,:);
                plotcolors = colors(discphi(i-ntail:i,:),:);
                plotrevolutions = revolutions(i,:) - revolutions(i-ntail,:); %How many times does the orbit wrap around in current plot for each particle
                
                plot(ax1,plotx,ploty, 'Color', 'k');
                hold on;
                scatter(ax1,reshape(plotx, [numel(plotx) 1]), reshape(ploty, [numel(ploty) 1]), 5, reshape(plotcolors, [size(plotcolors,1) 3]));
                if sum(plotrevolutions)
                for j = 1:numel(plotrevolutions)
                    plot(ax1,plotx(:,j)-poloidalextent*plotrevolutions(j),ploty(:,j), 'Color', 'k');
                    %scatter(ax1,plotx(:,j)-poloidalextent*plotrevolutions(j), ploty(:,j), 5, reshape(plotcolors, [size(plotcolors,1) 3]));
                end
                end
                if pspace
                    
                end
                hold off;
                ax1.YLim = [0 1.5];
                ax1.XLim = [0 poloidalextent];
                %scatter(ax1, deg2rad(theta(i,:),rho(i,:));
                set(gca,'FontSize',24);
                title('Particle Orbits');
                frame = getframe(fig); %for writing
                writeVideo(v,frame);
            end
        else
            plot(ax1,theta,rho);
            set(gca,'FontSize',24);
            title('Particle Orbits');
        end
    case {'vll','pll'}
       
        % Make plot
        
        plot(ax1,pll./amu,pperp./amu);
        set(gca,'FontSize',24);
        xlabel('Parallel Momentum [u*m/s]');
        ylabel('Perpendicular Momentum [u*m/s]')
        title('Phase space');
    case 'rz'
        % Extract position
        r = h5read(a5file,[path '/r']);
        z = h5read(a5file,[path '/z']);
        r = r(dex);
        z = z(dex);
        % Reshape by particle index
        r = reshape(r,[nsteps npart]);
        z = reshape(z,[nsteps npart]);
        
        % Reorder in time
        for i=1:npart
            r(:,i) = r(idex(:,i),i);
            z(:,i) = z(idex(:,i),i);
        end
        % Make plot
        plot(ax1, r,z);
        set(gca,'FontSize',24);
        xlabel('R [m]');
        ylabel('Z [m]')
        title('R-Z Overview');
        
end
end

