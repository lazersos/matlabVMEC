function plot_wall_acc(wall_data,varargin)
%PLOT_WALL_ACC Plots an accelerated wall mesh
%    The PLOT_WALL_ACC subroutine plots information regarding the submesh
%    for accelerated wall models.  
%
%   Options:
%
%       PLOT_WALL_ACC(wall_data,'grid')  Plots the entire subgrid as a
%       series of boxes of varying color in cartesian coordinates.
%
%       PLOT_WALL_ACC(wall_data,'subgrid')  Plots only those grids which
%       contain data.
%
%       PLOT_WALL_ACC(wall_data,'subfaces')  Plots the subfaces labeled by
%       color according to their subgrid.
%
%       PLOT_WALL_ACC(wall_data,'block',33)  Plots a specific block with
%       the faces refereced in the block if any data is present in the
%       block.
%
%   See also READ_WALL.
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           2/12/23

% Defaults
lplotgrid = 0;
lplotsubgrid = 0;
lplotsubfaces = 0;
plotblockn = 0;

% Handle Varargin
if nargin > 1
    i = 1;
    while i <= length(varargin)
        switch varargin{i}
            case 'grid'
                lplotgrid = 1;
            case 'subgrid'
                lplotsubgrid = 1;
            case 'subfaces'
                lplotsubfaces = 1;
            case 'block'
                i = i + 1;
                plotblockn=varargin{i};
        end
        i=i+1;
    end
end

% Plot grid
if or(lplotgrid,lplotsubgrid)
    colors = lines(wall_data.nblocks)';
    for i=1:wall_data.nblocks
        if and(wall_data.nblockfaces(i) <=0,lplotsubgrid), continue; end
        xmin = wall_data.blockbounds(1,i);
        xmax = wall_data.blockbounds(2,i);
        ymin = wall_data.blockbounds(3,i);
        ymax = wall_data.blockbounds(4,i);
        zmin = wall_data.blockbounds(5,i);
        zmax = wall_data.blockbounds(6,i);
        coords = [xmin ymin zmin;...
            xmax ymin zmin;...
            xmax ymax zmin;...
            xmin ymax zmin;...
            xmin ymin zmax;...
            xmax ymin zmax;...
            xmax ymax zmax;...
            xmin ymax zmax];
        faces = [1 2 3 4 1;
            1 2 6 5 1;
            2 3 7 6 2;
            3 4 8 7 3;
            4 1 5 8 4;
            5 6 7 8 5];
        ha = patch('vertices',coords,'faces',faces);
        set(ha,'FaceColor',colors(:,i), 'facealpha', 0.1);
    end
end

% Plot subfaces
if lplotsubfaces
    colors = lines(wall_data.nblocks)';
    for i=1:wall_data.nblocks
        if wall_data.nblockfaces(i) <=0, continue; end
        idx = wall_data.blockfaces(i).faces;
        ha = patch('vertices',wall_data.coords','faces',wall_data.faces(:,idx)');
        set(ha,'FaceColor',colors(:,i), 'facealpha', 1.0);
    end
end

% Plot a specific subgrid and model
if any(plotblockn > 0)
    for i = plotblockn
        xmin = wall_data.blockbounds(1,i);
        xmax = wall_data.blockbounds(2,i);
        ymin = wall_data.blockbounds(3,i);
        ymax = wall_data.blockbounds(4,i);
        zmin = wall_data.blockbounds(5,i);
        zmax = wall_data.blockbounds(6,i);
        coords = [xmin ymin zmin;...
            xmax ymin zmin;...
            xmax ymax zmin;...
            xmin ymax zmin;...
            xmin ymin zmax;...
            xmax ymin zmax;...
            xmax ymax zmax;...
            xmin ymax zmax];
        faces = [1 2 3 4 1;
            1 2 6 5 1;
            2 3 7 6 2;
            3 4 8 7 3;
            4 1 5 8 4;
            5 6 7 8 5];
        ha = patch('vertices',coords,'faces',faces);
        set(ha,'FaceColor','red', 'facealpha', 0.1);
        hold on;
        if wall_data.nblockfaces(i) >0
            idx = wall_data.blockfaces(i).faces;
            ha = patch('vertices',wall_data.coords','faces',wall_data.faces(:,idx)');
            set(ha,'FaceColor','b', 'facealpha', 1.0);
        end
    end
end

