function make_PSO_mov( xvar ,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Defaults
plottype='scatter';
ptime=0.2;
lmovie = 0;
        
i=1;
if nargin>=4
    n1 = varargin{1};
    n2 = varargin{2};
    n3 = varargin{3};
    i=4;
end
if i <= (nargin-1)
    while i<= (nargin-1)
        switch varargin{i}
            case {'scatter','bar','trail','scatter2','scatter_good'}
                plottype=varargin{i};
            case {'1080','720','480'}
                restype=varargin{i};
            case {'movie'}
                lmovie = 1;
        end
        i=i+1;
    end
end
% Set Color limits
chi_max = xvar.chisq(1);
chi_min = min(xvar.chisq);
if chi_max == chi_min
    chi_min = chi_max - 0.25*chi_max;
    chi_max = chi_max + 0.25*chi_max;
end
if chi_max > 1.5*chi_min
    chi_max = 1.5*chi_max;
end
xmax = max(xvar.x(n1,:));
xmin = min(xvar.x(n1,:));
ymax = max(xvar.x(n2,:));
ymin = min(xvar.x(n2,:));
zmax = max(xvar.x(n3,:));
zmin = min(xvar.x(n3,:));
if xmax == xmin
    xmin = xmax - 0.1*xmax;
    xmax = xmax + 0.1*xmax;
end
if ymax == ymin
    ymin = ymax - 0.1*ymax;
    ymax = ymax + 0.1*ymax;
end
if zmax == zmin
    zmin = zmax - 0.1*zmax;
    zmax = zmax + 0.1*zmax;
end


% Setup Movie File
if (lmovie)
    mov=VideoWriter([plottype '_movie.avi']);
    mov.Quality= 100;
    mov.FrameRate=10;
    open(mov);
end

switch plottype
    case{'scatter'}
        n_max = max(xvar.iter);
        point_size = xvar.chisq.*0.0+300.;
        min_dex = xvar.chisq.*0.0  ~= 0.0;
        min_global = xvar.chisq(1);
        for i=0:n_max
            dex = xvar.iter < i;
            ha=scatter3(xvar.x(n1,dex),xvar.x(n2,dex),xvar.x(n3,dex),point_size(dex),xvar.chisq(dex),'o');
            set(ha,'Linewidth',get(ha,'Linewidth')./i);
            caxis([chi_min chi_max]);
            dex = xvar.iter == i;
            min_local =  min(xvar.chisq(dex));
            if (min_local < min_global)
                min_global = min_local;
            end
            temp = xvar.chisq == min_global;
            min_dex = (min_dex + temp) == 1;
            hold on;
            scatter3(xvar.x(n1,dex),xvar.x(n2,dex),xvar.x(n3,dex),point_size(dex),xvar.chisq(dex),'.');
            plot3(xvar.x(n1,min_dex),xvar.x(n2,min_dex),xvar.x(n3,min_dex),'+','MarkerSize',20);
            hold off;
            pause(ptime);
        end
    case{'scatter2'}
        n_max = max(xvar.iter);
        point_size = xvar.chisq.*0.0+300.;
        min_dex = xvar.chisq.*0.0  ~= 0.0;
        min_global = xvar.chisq(1);
        for i=0:n_max
            dex = xvar.iter < i;
            dex2 = xvar.iter > i-10;
            dex = (dex+dex2)>1;
            ha=scatter3(xvar.x(n1,dex),xvar.x(n2,dex),xvar.x(n3,dex),point_size(dex),xvar.chisq(dex),'o');
            xlim([xmin xmax]);
            ylim([ymin ymax]);
            zlim([zmin zmax]);
            set(ha,'Linewidth',get(ha,'Linewidth')./i);
            caxis([chi_min chi_max]);
            dex = xvar.iter == i;
            min_local =  min(xvar.chisq(dex));
            if (min_local < min_global)
                min_global = min_local;
            end
            temp = xvar.chisq == min_global;
            min_dex = (min_dex + temp) == 1;
            hold on;
            scatter3(xvar.x(n1,dex),xvar.x(n2,dex),xvar.x(n3,dex),point_size(dex),xvar.chisq(dex),'.');
            plot3(xvar.x(n1,min_dex),xvar.x(n2,min_dex),xvar.x(n3,min_dex),'+','MarkerSize',20);
            hold off;
            pause(ptime);
        end
    case{'trail'}
        n_max = max(xvar.iter);
        nproc = sum(xvar.iter == 0);
        for i=0:n_max
            dex = xvar.iter <= i;
            ilast = sum(dex)/nproc;
            x = reshape(xvar.x(n1,dex),[nproc ilast]);
            y = reshape(xvar.x(n2,dex),[nproc ilast]);
            z = reshape(xvar.x(n3,dex),[nproc ilast]);
            plot3(x',y',z','-o');
            xlim([xmin xmax]);
            ylim([ymin ymax]);
            zlim([zmin zmax]);
            %for j=1:nproc
            %    if j==1, hold on; end
            %    plot3(xvar.x(n1,dex(j:nproc:ilast)),xvar.x(n2,dex(j:nproc:ilast)),xvar.x(n3,dex(j:nproc:ilast)),'.-');
            %end
            pause(ptime);
        end
    case{'scatter_good'}
        n_max = max(xvar.iter);
        point_size = xvar.chisq.*0.0+300.;
        min_dex = xvar.chisq.*0.0  ~= 0.0;
        min_global = xvar.chisq(1);
        dex2 = xvar.chisq < 1.0E10;
        chi_max = 100.*min(xvar.chisq(dex2));
        for i=0:n_max
            dex = xvar.iter < i;
            dex = (dex+dex2) > 1;
            ha=scatter3(xvar.x(n1,dex),xvar.x(n2,dex),xvar.x(n3,dex),point_size(dex),xvar.chisq(dex),'o');
            set(ha,'Linewidth',get(ha,'Linewidth')./i);
            caxis([chi_min chi_max]);
            dex = xvar.iter == i;
            min_local =  min(xvar.chisq(dex));
            if (min_local < min_global)
                min_global = min_local;
            end
            temp = xvar.chisq == min_global;
            min_dex = (min_dex + temp) == 1;
            hold on;
            scatter3(xvar.x(n1,dex),xvar.x(n2,dex),xvar.x(n3,dex),point_size(dex),xvar.chisq(dex),'.');
            plot3(xvar.x(n1,min_dex),xvar.x(n2,min_dex),xvar.x(n3,min_dex),'+','MarkerSize',20);
            hold off;
            xlim([xmin xmax]);
            ylim([ymin ymax]);
            zlim([zmin zmax]);
            pause(ptime);
            if (lmovie)
                frame=getframe(gcf);
                writeVideo(mov,frame);
            end
        end
    case{'bar'}
        n_max = max(xvar.iter);
        for i=0:n_max
            dex = xvar.iter == i;
            h=bar3(xvar.x(:,dex),'detached');
            color_data = xvar.chisq(dex);
            for j=1:length(h)
                zdata = get(h(j),'ZData');
                zdata(:) = color_data(j);
                set(h(j),'CData',zdata);
            end
            caxis([chi_min chi_max]);
            title(['Generation  ' num2str(i,'%3d')]);
            pause(ptime);
        end
        
end
if (lmovie)
    close(mov);
    % Make a poster frame
    saveas(gca,[plottype,'_poster.fig'],'fig');
end

end

