function make_beam_movie(varargin)
%MAKE_BEAM_MOVIE(data) Plots beams3d movies.
%   The MAKE_BEAM_MOVIE routine creates movies of the various particle
%   trajectories as calculated by the BEAMS3D code.  It takes a beams3d
%   data structure as read by the READ_BEAS3D routine as a mandatory
%   argument.  Optional arguments may be passed which modify what type of
%   movie is created.
%
%   Options:
%       'ntail':    Followed by a numer determines the number of points
%                   defining the 'tail' of a particle. (default 10)
%       'parts':    Followed by an array defining a subset of particles to
%                   follow.
%       'movie':    Output an AVI file. (slows code down)
%       'flux':     Particle orbits are plotted in 'flux' space.  Requires
%                   the code to have been run with the -flux flag, or
%                   produced by STELLOPT.
%       'flux_notail':   Same as 'flux' but with no trailing particles
%                        plotted.
%
%   Usage:
%       beam_data=read_beams3d('beams3d_test.h5');
%       make_beam_movie(beam_data,'ntail',5,'flux');
%
%   See also read_beams3d.
%
%   Created by: S. Lazerson (lazerson@pppl.gov)
%   Version:    1.0
%   Date:       04/10/15

% Set some defaults
ntail=10;  % Number of particles to use in creating fading tail
color=[1 0 0];  % Color of particles (in this case red)
parts=[];    % Subset of particle indicies to plot
bgcolor='k';   % background color
data=[];
vmec_data=[];
ha=[];
lmovie=0;
lflux=0;
ldist=0;
lflux_notail=0;
lbounce_notail=0;
lbounce=0;
movname='beams3d';
dt = 0.01;

% Read the input
if nargin>0
    i=1;
    while i<=nargin
        if isfield(varargin{i},'datatype')
            switch varargin{i}.datatype
                case 'BEAMS3D'
                    data=varargin{i};
                case 'wout'
                    vmec_data=varargin{i};
                case 'coil_data'
                    coil_data=varargin{i};
            end
        else
            switch varargin{i}
                case {'ntail'}
                    i=i+1;
                    ntail = varargin{i};
                case {'parts'}
                    i=i+1;
                    parts=varargin{i};
                case {'movie'}
                    lmovie=1;
                case {'flux'}
                    lflux=1;
                case {'flux_notail'}
                    lflux_notail=1;
                case {'bounce'}
                    lbounce=1;
                case {'bounce_notail'}
                    lbounce_notail=1;
                case {'dist'}
                    ldist=1;
                case {'dt'}
                    i=i+1;
                    dt=varargin{i};
            end
        end
        i=i+1;
    end
else
    disp('ERROR: Specify movie type');
    return
end

% Reset some things
if (isempty(parts))
    parts=1:size(data.X_lines,2);    % Subset of particle indicies to plot
end
total_steps=size(data.X_lines,1);  % Total number of timesteps to output
if (~isempty(vmec_data))
    theta = 0:2*pi./90:2*pi;
    zeta = 0:2*pi./90:2*pi;
    r=cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
    z=sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
end

win_res=[1 1 1920 1080];         % Movie resolution
fig=figure('DoubleBuffer','on','Position',win_res,...
        'Color',bgcolor,'BackingStore','on','MenuBar','none',...
        'Name','Fieldlines','InvertHardcopy','off');
set(gca,'nextplot','replacechildren','XColor',bgcolor,...
        'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
xmax=max(data.raxis);
ymax=xmax;
xmin=-xmax;
ymin=-xmax;
zmax=max(data.zaxis);
zmin=min(data.zaxis);

% Handle movie output
if (lmovie)
    mov=VideoWriter([movname '_movie.avi']);
    mov.Quality= 100;
    mov.FrameRate=30;
    open(mov);
end

% Plot the starting points
if (lflux)
    for i=1:ntail
        d2=i;
        d1=i-ntail;
        if (d1 < 1), d1=1; end
        set(gca,'nextplot','replacechildren','XColor',bgcolor,...
            'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
        hold on
        for j=d1:d2
            shade = double(j-1)/double(d2);
            hline=polar(data.U_lines(j,parts),data.S_lines(j,parts),'.');
            set(hline,'Color',color.*shade,'MarkerSize',18);
        end
        polar(0:2*pi/90:2*pi,ones(1,91).*0.25,'w--')
        polar(0:2*pi/90:2*pi,ones(1,91).*0.5,'w--')
        polar(0:2*pi/90:2*pi,ones(1,91).*0.75,'w--')
        polar(0:2*pi/90:2*pi,ones(1,91),'w')
        hold off
        xlim([-1 1]);
        ylim([-1 1]);
        axis square;
        pause(dt);
        if lmovie
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    end
    
    % Now plot the rest
    for i=ntail+1:total_steps-1
        d2=i;
        d1=d2-ntail;
        cla;
        set(gca,'nextplot','replacechildren','XColor',bgcolor,...
            'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
        hold on;
        for j=d1:d2
            shade = double(j-d1)/double(ntail);
            hline=polar(data.U_lines(j,parts),data.S_lines(j,parts),'.');
            set(hline,'Color',color.*shade,'MarkerSize',18);
        end
        polar(0:2*pi/90:2*pi,ones(1,91).*0.25,'w--')
        polar(0:2*pi/90:2*pi,ones(1,91).*0.5,'w--')
        polar(0:2*pi/90:2*pi,ones(1,91).*0.75,'w--')
        polar(0:2*pi/90:2*pi,ones(1,91),'w')
        hold off
        xlim([-1 1]);
        ylim([-1 1]);
        axis square;
        pause(dt);
        if lmovie
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    end
elseif (lflux_notail)
    % Now plot the rest
    for i=1:total_steps-1
        cla;
        set(gca,'nextplot','replacechildren','XColor',bgcolor,...
            'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
        hold on;
        hline=polar(data.U_lines(i,parts),data.S_lines(i,parts),'.');
        set(hline,'Color',color,'MarkerSize',18);
        polar(0:2*pi/90:2*pi,ones(1,91).*0.25,'w--')
        polar(0:2*pi/90:2*pi,ones(1,91).*0.5,'w--')
        polar(0:2*pi/90:2*pi,ones(1,91).*0.75,'w--')
        polar(0:2*pi/90:2*pi,ones(1,91),'w')
        hold off
        xlim([-1 1]);
        ylim([-1 1]);
        axis square;
        pause(dt);
        if lmovie
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    end
elseif (lbounce)
    vll = data.vll_lines;
    vll = abs(vll./max(abs(vll(1,:))));
    for i=1:ntail
        d2=i;
        d1=i-ntail;
        if (d1 < 1), d1=1; end
        set(gca,'nextplot','replacechildren','XColor',bgcolor,...
            'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
        hold on
        for j=d1:d2
            shade = double(j-1)/double(d2);
            plot(mod(data.U_lines(j,parts),2*pi),vll(j,parts),'.','Color',color.*shade,'MarkerSize',18);
        end
        plot([0 2*pi],[1 1]*0.50,'w--');
        plot([0 2*pi],[1 1]*1.00,'w--');
        plot([0 2*pi],[1 1]*1.50,'w--');
        plot([0 2*pi],[1 1]*2.00,'w-');
        plot([0 2*pi],[0 0],'w-');
        plot([0 0],[0 2],'w-');
        plot([1 1].*2*pi,[0 2],'w-');
        hold off
        xlim([0 2*pi]);
        ylim([0 2]);
        axis square;
        pause(dt);
        if lmovie
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    end
    
    % Now plot the rest
    for i=ntail+1:total_steps-1
        d2=i;
        d1=d2-ntail;
        cla;
        set(gca,'nextplot','replacechildren','XColor',bgcolor,...
            'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
        hold on;
        for j=d1:d2
            shade = double(j-d1)/double(ntail);
            plot(mod(data.U_lines(j,parts),2*pi),vll(j,parts),'.','Color',color.*shade,'MarkerSize',18);
        end
        plot([0 2*pi],[1 1]*0.50,'w--');
        plot([0 2*pi],[1 1]*1.00,'w--');
        plot([0 2*pi],[1 1]*1.50,'w--');
        plot([0 2*pi],[1 1]*2.00,'w-');
        plot([0 2*pi],[0 0],'w-');
        plot([0 0],[0 2],'w-');
        plot([1 1].*2*pi,[0 2],'w-');
        hold off
        xlim([0 2*pi]);
        ylim([0 2]);
        axis square;
        pause(dt);
        if lmovie
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    end
elseif (lbounce_notail)
    % Now plot the rest
    vll = data.vll_lines;
    vll = abs(vll./max(abs(vll(1,:))));
    for i=1:total_steps-1
        cla;
        set(gca,'nextplot','replacechildren','XColor',bgcolor,...
            'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
        hold on;
        %plot(mod(data.V_lines(i,parts),2*pi),data.S_lines(i,parts),'.','Color',color,'MarkerSize',18);
        plot(mod(data.U_lines(i,parts),2*pi),vll(i,parts),'.','Color',color,'MarkerSize',18);
        plot([0 2*pi],[1 1]*0.50,'w--');
        plot([0 2*pi],[1 1]*1.00,'w--');
        plot([0 2*pi],[1 1]*1.50,'w--');
        plot([0 2*pi],[1 1]*2.00,'w-');
        plot([0 2*pi],[0 0],'w-');
        plot([0 0],[0 2],'w-');
        plot([1 1].*2*pi,[0 2],'w-');
        hold off
        xlim([0 2*pi]);
        ylim([0 2]);
        axis square;
        pause(dt);
        if lmovie
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    end
elseif (ldist)
    %for j=1:data.nphi
    %   BX(j,:,:) = data.B_R(:,j,:).*cos(data.phiaxis(j))-data.B_PHI(:,j,:).*sin(data.phiaxis(j));
    %   BY(j,:,:) = data.B_R(:,j,:).*sin(data.phiaxis(j))+data.B_PHI(:,j,:).*cos(data.phiaxis(j));
    %   BZ(j,:,:) = data.B_Z(:,j,:);
    %end
    %B = sqrt(BX.^2 + BY.^2 + BZ.^2);
    %[x1, y1, z1] = meshgrid(min(data.raxis):(max(data.raxis)-min(data.raxis))./(length(data.raxis)-1):max(data.raxis),...
    %     min(data.phiaxis):(max(data.phiaxis)-min(data.phiaxis))./(length(data.phiaxis)-1):max(data.phiaxis),...
    %     min(data.zaxis):(max(data.zaxis)-min(data.zaxis))./(length(data.zaxis)-1):max(data.zaxis));
    mass = 1.6726231E-27;
    ec=1.60217733E-19;
    vperp_lines = sqrt(abs(2.*data.moment_lines.*data.B_lines./mass));
    ymax = max(max(abs(vperp_lines)))./1E6;
    xmax = max(max(abs(data.vll_lines)))./1E6;
    for j=1:total_steps-1
        cla;
        set(gca,'nextplot','replacechildren','XColor','w',...
            'YColor','w','ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
        dex = data.S_lines(j,:) <1;
        hold on;
        plot(data.vll_lines(j,dex)./1E6,vperp_lines(j,dex)./1E6,'.w');
        hold off;
        axis square;
        ylim([0 5]);
        xlim([-5 5]);
        xlabel('v_{||} x10^6 [m/s]');
        ylabel('v_\perp x10^6 [m/s]');
        pause(dt);
        if lmovie
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    end
else
    camva('manual');
    camtarget('manual');
    for i=1:ntail
        d2=i;
        d1=i-ntail;
        if (d1 < 1), d1=1; end
        %plot3(data.X_lines(d1:d2,parts),data.Y_lines(d1:d2,parts),data.Z_lines(d1:d2,parts),'Color',grey);
        set(gca,'nextplot','replacechildren','XColor',bgcolor,...
            'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
        axis off;
        hold on
        for j=d1:d2
            shade = double(j-1)/double(d2);
            plot3(data.X_lines(j,parts),data.Y_lines(j,parts),data.Z_lines(j,parts),'.','Color',color.*shade,'MarkerSize',18);
        end
        if (~isempty(vmec_data))
            delete(ha);
            ha=isotoro(r,z,zeta,vmec_data.ns);
            set(ha,'FaceAlpha',0.2);
        end
        hold off
        xlim([xmin xmax]);
        ylim([ymin ymax]);
        zlim([zmin zmax]);
        axis equal;
        axis off;
        camtarget([0 0 0])
        view(90,30);
        pause(dt);
        if lmovie
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    end
    
    % Now plot the rest
    for i=ntail+1:total_steps-1
        d2=i;
        d1=d2-ntail;
        %plot3(data.X_lines(1:d2,parts),data.Y_lines(1:d2,parts),data.Z_lines(1:d2,parts),'Color',grey);
        cla
        set(gca,'nextplot','replacechildren','XColor',bgcolor,...
            'YColor',bgcolor,'ZColor',bgcolor,'Color',bgcolor,'FontSize',20);
        hold on;
        %disp(num2str([d1 d2],'%i'));
        for j=d1:d2
            shade = double(j-d1)/double(ntail);
            plot3(data.X_lines(j,parts),data.Y_lines(j,parts),data.Z_lines(j,parts),'.','Color',color.*shade,'MarkerSize',18);
        end
        %plot3(x(99,:,1),y(99,:,1),z(99,:,1),'r');
        %plot3(x(99,:,2),y(99,:,2),z(99,:,2),'r');
        %plot3(x(99,:,3),y(99,:,3),z(99,:,3),'r');
        if (~isempty(vmec_data))
            ha=isotoro(r,z,zeta,vmec_data.ns);
            set(ha,'FaceAlpha',0.2);
        end
        hold off
        xlim([xmin xmax]);
        ylim([ymin ymax]);
        zlim([zmin zmax]);
        axis equal;
        axis off;
        %view(90+180.*fieldline.PHI_lines(d2,)./pi,30);
        pause(dt);
        if lmovie
            frame=getframe(fig);
            writeVideo(mov,frame);
        end
    end
end

if lmovie
    close(mov);
end


end

