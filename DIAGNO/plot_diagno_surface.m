function plot_diagno_surface(diagno_data,varargin)
%PLOT_DIAGNO_SURFACE Plots DIAGNO surface geometry data
%   This function plots the DIAGNO surface geometry file.  If an axis file
%   is read and passed to the function the axis is plotted in 3D as a line
%   segment.
%
%   Options:
%       plot_diagno_surface(diagno_data,'surf') Plots all the points
%       making up the DIAGNO surface in x,y,z.  A cone plot is also plotted
%       on the surface in the direction of the surface normals. (default)
%
%       plot_diagno_surface(diagno_data,'phidf') Plots all the points
%       making up the DIAGNO surface in x,y,z.  A cone plot is also plotted
%       on the surface in the direction of the vector Phi*df.
%   
%   Example:
%       dsurf_data=read_diagno_surface('plasma_surf.geom');
%       plot_diagno_surface(dsurf_data);
%       asurf_data=read_diagno_surface('plasma_axis.geom');
%       plot_diagno_surface(asurf_data);
%
%   See also read_diagno_surface.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           1/04/11  

% Set defaults
plottype='surf';
stepsize=45;
markercolor='blue';
markersize=2;
arrow_linewidth=4;
% Handle varargin
numdefargs=1;   %Number of default arguments
if nargin >numdefargs
    for i=1:nargin-numdefargs
        switch varargin{i}
            case {'surf_norm','surf','phidf','bexn'}
                plottype=varargin{i};
        end
    end
end
% Check to see if it's a axis file
if isfield(diagno_data,'datatype')
    switch diagno_data.datatype
        case 'diagno_surf'
            % Get array info
            nels=size(diagno_data.data,1);
            nlines=size(diagno_data.data,2);
            % Now plot
            vertex=1:stepsize:nlines;
            switch plottype
                case 'surf'
                    plot3(diagno_data.data(4,:),diagno_data.data(5,:),diagno_data.data(6,:),'.',...
                        'Color',markercolor);
                case 'surf_norm'
                    plot3(diagno_data.data(4,:),diagno_data.data(5,:),diagno_data.data(6,:),'.',...
                        'Color',markercolor);
                    hold on
                    coneplot(diagno_data.data(4,vertex),diagno_data.data(5,vertex),diagno_data.data(6,vertex),...
                        diagno_data.data(7,vertex),diagno_data.data(8,vertex),diagno_data.data(9,vertex),...
                        'nointerp','quiver');
                    hold off
                case 'phidf'
                    plot3(diagno_data.data(4,:),diagno_data.data(5,:),diagno_data.data(6,:),'.',...
                        'Color',markercolor);
                    hold on
                    coneplot(diagno_data.data(4,vertex),diagno_data.data(5,vertex),diagno_data.data(6,vertex),...
                        diagno_data.data(10,vertex),diagno_data.data(11,vertex),diagno_data.data(12,vertex),...
                        'nointerp','quiver');
                    hold off
                case 'bexn'
                    plot(diagno_data.data(13,:))
            end
        case 'diagno_axis'
            hold on
            plot3(diagno_data.data(2,:),diagno_data.data(3,:),diagno_data.data(4,:),'Color','black')
            hold off
    end
    axis vis3d; 
else
    disp('ERROR: No datatype field found');
end
end

