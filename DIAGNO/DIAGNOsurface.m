function DIAGNOsurface(varargin)
%DIAGNOsurface Plots DIAGNO surface geometry data
%   This program plots the plasma_surface.geom file produced by DIAGNO.  

% Set defaults
plottype='surf';
filename='plasma_surface.geom';
stepsize=45;
markercolor='blue';
markersize=2;
arrow_linewidth=4;
% Handle varargin
numdefargs=0;   %Number of default arguments
if nargin >numdefargs
    for i=1:nargin-numdefargs
        switch varargin{i}
            case {'surf','phidf','bexn'}
                plottype=varargin{i};
            case 'filename'
                i=i+1;
                filename=varargin{i};
        end
    end
end
% Read file
fid=fopen(filename);
headline=fgetl(fid);                %Read Header
line=fgetl(fid);                    %Read First Line
data=sscanf(line,'%e');
while ~feof(fid)
    line=fgetl(fid);
    if ~strcmp(line,'###');
        data=[data sscanf(line,'%e')];
    end
end
fclose(fid);
% Get array info
nels=size(data,1);
nlines=size(data,2);
% Now plot
fig=figure;
vertex=1:stepsize:nlines;
switch plottype
    case 'surf'
        plot3(data(4,:),data(5,:),data(6,:),'.',...
            'Color',markercolor);
        hold on
        coneplot(data(4,vertex),data(5,vertex),data(6,vertex),...
            data(7,vertex),data(8,vertex),data(9,vertex),...
            'nointerp','quiver');
        hold off
    case 'phidf'
        plot3(data(4,:),data(5,:),data(6,:),'.',...
            'Color',markercolor);
        hold on
        coneplot(data(4,vertex),data(5,vertex),data(6,vertex),...
            data(10,vertex),data(11,vertex),data(12,vertex),...
            'nointerp','quiver');
        hold off
    case 'bexn'
        plot(data(13,:))
end
axis vis3d; 
end

