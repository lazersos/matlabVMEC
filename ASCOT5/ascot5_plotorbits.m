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
ec = 1.60217662E-19;
plottype='xyz';
parts = [];
lmovie = 0;
pspace = 0;
mirspace = 0;
lgc = 0;
filename = 'default';
startframe = 1; %Number of Start-Frame - ntail
endframe = []; %Number of last frame
ntail = 20;

fig=[];

%Deal with multiple input files (have to have same marker ids)
if ~iscell(a5file)
    files = {a5file};
else
    files = a5file;
end

% Check for file(s)
for i = 1:size(files,2)
    a5file = files{i};
    if ~isfile(a5file)
        disp(['ERROR: ' a5file ' file not found!']);
        return;
    end
end

% Handle varargin
if nargin > 2
    i=1;
    while i <= nargin-2
        switch varargin{i}
            case{'xyz','flux','vll','pll', 'rz'}
                plottype=varargin{i};
            case{'fluxperp'}
                plottype=varargin{i};
                if ~mirspace
                    pspace = 1;
                end
            case{'cartesian'}
                plottype=varargin{i};
                pspace = 1;
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
            case {'bref', 'b_ref', 'bmir', 'b_mir'}
                mirspace =1;
                if pspace
                    pspace = 0;
                end
            case 'startframe'
                i = i+1;
                startframe = varargin{i};
            case 'endframe'
                i = i+1;
                endframe = varargin{i};
                
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end


for f = 1:size(files,2)
    
    % Use active run
    if isempty(runid)
        try
            runid=h5readatt(files{f},'/results','active');
            disp(['  Using runid: ' runid]);
        catch
            runid=[];
            return;
        end
    end
    % Set run path
    path = ['/results/run_' num2str(runid,'%10.10i') '/orbit'];
    %Check if gyro-center
    optionsid = h5readatt(files{f},'/options', 'active');
    path_simmode = ['/options/opt_' num2str(optionsid,'%10.10i')];
    simmode = h5read(files{f}, [path_simmode '/SIM_MODE']);
    
    if simmode == 2
        [time{f}, r{f}, phi{f}, z{f}, br{f}, bphi{f}, bz{f}, pll{f}, mu{f}, rho{f}, theta{f} ] = ascot5_get_properties_time_trace(files{f}, path, parts, 'r', 'phi', 'z', 'br', 'bphi', 'bz', 'ppar', 'mu', 'rho', 'theta');
    else
        [time{f}, r{f}, phi{f}, z{f}, br{f}, bphi{f}, bz{f}, pr{f}, pphi{f}, pz{f}, rho{f}, theta{f} ] = ascot5_get_properties_time_trace(files{f}, path, parts, 'r', 'phi', 'z', 'br', 'bphi', 'bz', 'pr', 'pphi', 'pz', 'rho', 'theta');
        
    end
    
    x{f} = r{f}.*cosd(phi{f}); %Phi is in degrees
    y{f} = r{f}.*sind(phi{f});
    
    xyz_mat = [x{f}, y{f}, z{f}];
    writematrix(xyz_mat, 'xyz_particle_NBI_19.txt', 'Delimiter', ';');
    
    rhocostheta{f} = rho{f}.*cosd(theta{f});
    rhosintheta{f} = rho{f}.*sind(theta{f});
    
    % get pll
    b{f}=sqrt(br{f}.*br{f}+bphi{f}.*bphi{f}+bz{f}.*bz{f});
    diffb{f} = [diff(b{f}); 0];
    if simmode ~= 2
        pll{f} = (pr{f}.*br{f}+pphi{f}.*bphi{f}+pz{f}.*bz{f})./b{f};
        p2{f}  = pr{f}.*pr{f}+pphi{f}.*pphi{f}+pz{f}.*pz{f};
        pperp{f} = sqrt(p2{f}-pll{f}.*pll{f});
        vperp{f} = pperp{f}./amu;
    else
        vperp{f} = sqrt( b{f} .* 2 .* mu{f} ./ amu * ec);
    end
    
    vll{f} = pll{f}./amu; %ONLY CORRECT FOR HYDROGEN!!!
    bmir{f} = b{f}.* (vll{f}.^2 + vperp{f} .^2) ./ vperp{f}.^2;
    
    
    skipframes = ceil((1E-3./24./60)./(time{f}(4,1) - time{f}(3,1)));  %How many data points to skip between images, here 1ms/video minute, assuming 24fps
    if isempty(endframe) || (size(time{f},1)<endframe)
        endframe = size(time{f},1);
    end
    
end
%----------Calculate Orbits

%Prepare Figure
fig = figure('Position',get(0,'ScreenSize'));%figure('Position',[1 1 1024 768],'Color','white','InvertHardCopy','off');
ax1 = gca;
if pspace
    delete(gca);
    ax1 = axes('Position', [0.07 0.18 0.55 0.75]);
    ax2 = axes('Position', [.75 0.18 0.2*2/3.0 0.2]);
    ax3 = axes('Position', [.68 0.58 0.3 0.4]);
end
if mirspace
    delete(gca);
    ax1 = axes('Position', [0.07 0.18 0.55 0.75]);
    ax2 = axes('Position', [.75 0.18 0.2 0.2]);
    ax3 = axes('Position', [.68 0.58 0.3 0.4]);
    for f = 1:size(files,2)
        timediff{f}= repmat(time{f}(1,1), size(time{f},2), 1) - time{f}(1,:)';
        time{f} = time{f} + repmat(timediff{f}, 1, size(time{f},1))';
    end
end
%linecolors = num2cell(lines(size(x,2)),2);
for f = 1:size(files,2)
    linecolors{f} = lines(size(x{f},2));
end
if lmovie
v = VideoWriter([filename '.mp4'], 'MPEG-4');
open(v);
end
switch plottype
    case 'xyz'
        for f = 1:size(files,2)
            modphi = mod(phi,360);
            edges = linspace(0,360,360);
            discphi = discretize(modphi,edges);
        end
        
        colors = parula(360);
        % Make plot
        %plot3(ax1,x,y,z);
        hold on;
        if lmovie
            for i=startframe+ntail:skipframes:endframe
                for f = 1:size(files,2)
                    %s1 = scatter3(ax1,x(i-ntail:i,:),y(i-ntail:i,:),z(i-ntail:i,:));
                    s2 = plot3(ax1,x(i-ntail:i,:),y(i-ntail:i,:),z(i-ntail:i,:),'Color', colors(discphi(i),:));
                end
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
                ylim(ax1,[-6.9 6.9]);
                xlim(ax1,[-6.9 6.9]);
                
                title('Particle Orbits');
                
                frame = getframe(fig); %for writing
                writeVideo(v,frame);
                %delete(s1);
                delete(s2);
            end
        else
            for f = 1:size(files,2)
                scatter3(ax1,x,y,z,0.1,time);
            end
            set(gca,'FontSize',24);
            xlabel('X [m]');
            ylabel('Y [m]');
            zlabel('Z [m]');
            title('Particle Orbits');
            axis equal;
        end
    case 'flux'
        
        % Make plot
        delete(ax1);
        ax1 = polaraxes;
        for f = 1:size(files,2)
            modphi = mod(phi{f},180);
            edges{f} = linspace(0,180,360);
            discphi{f} = discretize(modphi,edges{f});
        end
        colors = parula(360);
        
        if lmovie
            for i=startframe+ntail:skipframes:endframe
                for f = 1:size(files,2)
                    polarplot(ax1,deg2rad(theta{f}(i-ntail:i,:)),rho{f}(i-ntail:i,:), 'Color', colors(discphi{f}(i),:));
                end
                ax1.RLim = [0 1.5];
                
                %scatter(ax1, deg2rad(theta(i,:),rho(i,:));
                set(gca,'FontSize',24);
                title('Particle Orbits');
                frame = getframe(fig); %for writing
                writeVideo(v,frame);
            end
        else
            for f = 1:size(files,2)
                polarplot(ax1,deg2rad(theta{f}),rho{f});
            end
            set(gca,'FontSize',24);
            title('Particle Orbits');
        end
    case 'fluxperp'
        poloidalextent = 360.0;
        toroidalextent = 72.0;
        colors = parula(360);
        for f = 1:size(files,2)
            edges = linspace(min(b{f}, [], 'all'), max(b{f}, [], 'all'), 360);
            discquant = discretize(b{f},edges);
            polrevolutions = floor(theta{f}./ poloidalextent); %How many full revolutions have the particles completed
            torrevolutions = floor(phi{f}./ toroidalextent);       
        end
        
        %phi = phi - torrevolutions(1,:).*toroidalextent; %Reset rotations
        %theta = theta - polrevolutions(1,:).*poloidalextent;
        endframe = size(b{f},1);
        if lmovie   
            tail = ntail;
        else
            tail = endframe - 1;
            startframe = 1;
%             tail = endframe - startframe;
%             startframe = startframe+ntail;
        end
        
        hold(ax1,'on');
        hold(ax2,'on');
        hold(ax3,'on');
        for i=startframe+tail:skipframes:endframe
            plottheta = theta{f}(i-tail:i,:);
            plotrho = rho{f}(i-tail:i,:);
            plotcolors = colors(discquant(i-tail:i,:),:);
            plot_polrevolutions = [polrevolutions(i-tail,:)' polrevolutions(i,:)']; %How many times does the orbit wrap around in current plot for each particle
            plot_torrevolutions = [torrevolutions(i-tail,:)' torrevolutions(i,:)'];
            %                 h1 = plot(ax1,plotx,ploty);
            %                 set(h1, {'color'}, linecolors);
            %                 hold(ax1,'on');
            
            %if sum(plotrevolutions) %replot when part of orbit goes offscreen
            for j = 1:numel(parts) %Plot all frames where particle goes over plot area boundary
                if plot_polrevolutions(j, 1)>plot_polrevolutions(j, 2)
                    plot_polrevolutions = fliplr(plot_polrevolutions);
                end
                if plot_torrevolutions(j,1)>plot_torrevolutions(j,2)
                    plot_torrevolutions = fliplr(plot_torrevolutions);
                end
                for k = plot_polrevolutions(j, 1):plot_polrevolutions(j, 2)%Frames for poloidal crossings
                    plotx = plottheta(:,j)-poloidalextent*k;
                    ploty = plotrho(:,j);
                    scattercolors = plotcolors((j-1)*size(plottheta,1)+1:(j*size(plottheta,1)),:);
                    dex = plotx < 360 & 0 < plotx;
                    plotx = plotx(dex);
                    ploty = ploty(dex);
                    scattercolors = scattercolors(dex,:);
                    %set(h2(j), {'color'}, linecolors);
                    %scatter(ax1,reshape(plotx, [numel(plotx) 1])-poloidalextent*k, reshape(ploty, [numel(ploty) 1]), 5, reshape(plotcolors, [size(plotcolors,1) 3]));
                    
                    plot(ax1,plotx,ploty, 'Color', linecolors{f}(j,:));
                    scatter(ax1,plotx,ploty, 5, scattercolors);
                    
                    if pspace
                        plot(ax2,pll(i-tail:i,j)./amu,pperp(i-tail:i,j)./amu, 'LineWidth', 5, 'Color', linecolors{f}(j,:));
                        for m = plot_torrevolutions(j,1):plot_torrevolutions(j,1) %Frames for toroidal crossings
                            ploty = phi(i-tail:i,j)-toroidalextent*m;
                            ploty = ploty(dex);
                            %ploty = mod(ploty(dex),toroidalextent); %Mod should be removed by previous calculations
                            plot(ax3, plotx, ploty, 'Color', linecolors{f}(j,:));
                        end
                    end
                    
                    if mirspace
                        plot(ax2, time{f}(:,j), bmir{f}(:,j), 'Color', linecolors{f}(j,:));
                        plot(ax2, time{f}(i-tail:i,j), bmir{f}(i-tail:i,j));
                        plot(ax2, time{f}(i,j), bmir{f}(i,j), 'r.', 'MarkerSize', 36);
                        for m = plot_torrevolutions(j,1):plot_torrevolutions(j,1) %Frames for toroidal crossings
                            ploty = phi{f}(i-tail:i,j)-toroidalextent*m;
                            ploty = ploty(dex);
                            %ploty = mod(ploty(dex),toroidalextent); %Mod should be removed by previous calculations
                            plot(ax3, plotx, ploty, 'Color', linecolors{f}(j,:));
                        end
                    end
                    
                    %scatter(ax1,plotx(:,j)-poloidalextent*plotrevolutions(j), ploty(:,j), 5, reshape(plotcolors, [size(plotcolors,1) 3]));
                end
            end
            
        end
        if pspace
            %h3 = plot(ax2,pll(i-tail:i,:)./amu,pperp(i-tail:i,:)./amu, 'LineWidth', 5);
            %set(h3, {'color'}, linecolors);
            ax2.XLim = [-3E6 3E6];
            ax2.YLim = [0 4E6];
            xlabel(ax2,'Parallel Mom. [u*m/s]');
            ylabel(ax2,'Perp. Mom. [u*m/s]')
            %title(ax2,'Phase space');
            set(ax2,'FontSize',20);
            %h5 = plot(ax3, plotx(:,:), mod(phi(i-tail:i,:),toroidalextent));
            ax3.XLim = [0 poloidalextent];
            ax3.YLim = [0 toroidalextent];
            xlabel(ax3, 'Poloidal angle \theta [deg]');
            ylabel(ax3, 'Toroidal angle \phi [deg]');
            set(ax3,'FontSize',20);
            %set(h5, {'color'}, linecolors);
        end
        if mirspace
            %ax2.XLim = [-3E6 3E6];
            %ax2.YLim = [0 4E6];
            xlabel(ax2,'Time [s]');
            ylabel(ax2,'B_{mir} = E/\mu [T]')
            %title(ax2,'Phase space');
            set(ax2,'FontSize',20);
            %h5 = plot(ax3, plotx(:,:), mod(phi(i-tail:i,:),toroidalextent));
            ax3.XLim = [0 poloidalextent];
            ax3.YLim = [0 toroidalextent];
            xlabel(ax3, 'Poloidal angle \theta [deg]');
            ylabel(ax3, 'Toroidal angle \phi [deg]');
            set(ax3,'FontSize',20);
        end
        
        ax1.YLim = [0 1.5];
        ax1.XLim = [0 poloidalextent];
        ax1.CLim = [min(b{f}, [], 'all'), max(b{f}, [], 'all')];
        c = colorbar(ax1,'south');
        c.Label.String = 'B Field [T]';
        %scatter(ax1, deg2rad(theta(i,:),rho(i,:));
        set(ax1,'FontSize',24);
        
        xlabel(ax1, 'Poloidal angle \theta [deg]');
        ylabel(ax1, 'Flux label \rho [a.u.]');
        %title(ax1,'Particle Orbits');
        
        frame = getframe(fig); %for writing
        writeVideo(v,frame);
        if lmovie
            cla(ax1)
            if pspace || mirspace
                cla(ax2)
                cla(ax3)
            end
        end
        
    case 'cartesian'
        colors = parula(360);
        
        edges = linspace(min(b{1}, [], 'all'), max(b{1}, [], 'all'), 360);
        for f = 1:size(files,2)
            discquant{f} = discretize(b{f},edges);
        end
        if lmovie
            %endframe = size(b,1);
            tail = ntail;
        else
            %endframe = size(b,1);
            startframe = endframe-ntail;
            tail = endframe - 1;
            
        end
        
        hold(ax1,'on');
        hold(ax2,'on');
        hold(ax3,'on');
        for i=startframe+ntail:skipframes:endframe
            for f = 1:size(files,2)
                plotcolors{f} = colors(discquant{f}((i-tail):i,:),:);
            end
            for f = 1:size(files,2)
                for j = 1:numel(parts) %Plot all frames where particle goes over plot area boundary
                    plotx = y{f}(i-tail:i,j);
                    ploty =  x{f}(i-tail:i,j);
                    scattercolors = plotcolors{f}((j-1)*size(plotx,1)+1:(j*size(plotx,1)),:);
                    plot(ax1,plotx,ploty, 'Color', linecolors{f}(j,:));
                    scatter(ax1,plotx,ploty, 5, scattercolors);
                    
                    plot(ax2, time{f}(:,j), bmir{f}(:,j), 'Color', linecolors{f}(j,:));
                    plot(ax2, time{f}(i-tail:i,j), bmir{f}(i-tail:i,j));
                    plot(ax2, time{f}(i,j), bmir{f}(i,j), 'r.', 'MarkerSize', 36);
                    
                    %plot(ax3, r(i-tail:i,j), z(i-tail:i,j), 'Color', linecolors(j,:));
                    plot(ax3, rhosintheta{f}(i-tail:i,j), rhocostheta{f}(i-tail:i,j), 'Color', linecolors{f}(j,:));
                end
            end
            
            ax1.YLim = [-8 8];
            ax1.XLim = [-8 7];
            axis equal;
            ax1.CLim = [min(b{1}, [], 'all'), max(b{1}, [], 'all')];
            c = colorbar(ax1,'south');
            c.Label.String = 'B Field [T]';
            set(ax1,'FontSize',24);
            xlabel(ax1, 'X axis [m]');
            ylabel(ax1, 'Y axis [m]');
            
            xlabel(ax2,'Time [s]');
            ylabel(ax2,'B_{mir} = E/\mu [T]')
            set(ax2,'FontSize',20);
            
            %             ax3.XLim = [4.8 6.5];
            %             ax3.YLim = [-1.3 1.3];
            %             xlabel(ax3, 'R axis [m]');
            %             ylabel(ax3, 'Z Axis [m]');
            xlabel(ax3, '\rho*cos(\theta)');
            ylabel(ax3, '\rho*sin(\theta)');
            set(ax3,'FontSize',20);
            ax3.XLim = [-2 2];
            ax3.YLim = [-2 2];
            
            

            if lmovie
                frame = getframe(fig); %for writing
                 writeVideo(v,frame);
                cla(ax1)
                if pspace || mirspace
                    cla(ax2)
                    cla(ax3)
                end
            end
        end
        
    case {'vll','pll'}
        
        % Make plot
        
        plot(ax1,pll./amu,pperp./amu);
        set(gca,'FontSize',24);
        xlabel('Parallel Momentum [u*m/s]');
        ylabel('Perpendicular Momentum [u*m/s]')
        title('Phase space');
    case 'rz'
        % Make plot
        plot(ax1, r,z);
        set(gca,'FontSize',24);
        xlabel('R [m]');
        ylabel('Z [m]')
        title('R-Z Overview');
        
end


end