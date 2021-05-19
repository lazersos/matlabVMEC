function phandle = plot_wall_accelerated(data,varargin)
%PLOT_WALL_ACCELERATED wrapper of plot_wall
%   For accelerated walls
%
%  See also plot_divertor.
%
%   Written by:     D.J. Engels (d.j.engels@student.tue.nl))
%   Version:        1.0
%   Date:           05/21

%% Figure (disabled due to being slow)
% figure;
% if nargin>1
%     phandle = plot_divertor(data,varargin{:});
% else
%     phandle = plot_divertor(data);
% end
% hold on;
% % Plot boxes
% for i=1:data.nblocks
%     xs = [data.blocks(i).x_min, data.blocks(i).x_min, data.blocks(i).x_max, data.blocks(i).x_max, data.blocks(i).x_min, data.blocks(i).x_min, data.blocks(i).x_min, data.blocks(i).x_min, data.blocks(i).x_min, data.blocks(i).x_max, data.blocks(i).x_max, data.blocks(i).x_max, data.blocks(i).x_max, data.blocks(i).x_max, data.blocks(i).x_max, data.blocks(i).x_min];
%     ys = [data.blocks(i).y_min, data.blocks(i).y_max, data.blocks(i).y_max, data.blocks(i).y_min, data.blocks(i).y_min, data.blocks(i).y_min, data.blocks(i).y_max, data.blocks(i).y_max, data.blocks(i).y_max, data.blocks(i).y_max, data.blocks(i).y_max, data.blocks(i).y_max, data.blocks(i).y_min, data.blocks(i).y_min, data.blocks(i).y_min, data.blocks(i).y_min];
%     zs = [data.blocks(i).z_min, data.blocks(i).z_min, data.blocks(i).z_min, data.blocks(i).z_min, data.blocks(i).z_min, data.blocks(i).z_max, data.blocks(i).z_max, data.blocks(i).z_min, data.blocks(i).z_max, data.blocks(i).z_max, data.blocks(i).z_min, data.blocks(i).z_max, data.blocks(i).z_max, data.blocks(i).z_min, data.blocks(i).z_max, data.blocks(i).z_max];
%     plot3(xs, ys, zs)
% end
% title("Mesh split up in uniform grid");
% hold off;
% % histogram of how faces are spread over boxes
% figure;
% histogram([data.blocks.nfaces], data.nblocks)
% title("Spread of nfaces over the blocks generated");
% xlabel("Number of faces in a box");
% ylabel("Occurance");

%% Plot information
fprintf("----------- Accelerated mesh information -----------\n");
fprintf("Minimum number of faces in a block: %7d\n", min([data.blocks.nfaces]));
fprintf("Mean number of faces in a block: %5.1f\n", mean([data.blocks.nfaces]));
fprintf("Std of faces in a block: %5.1f\n", std([data.blocks.nfaces]));
fprintf("Maximum number of faces in a block: %7d\n", max([data.blocks.nfaces]));
fprintf("Number of blocks: %7d\n", data.nblocks);
fprintf("Total number of faces: %7d\n", data.nfaces);
fprintf("Total number of vertices: %7d\n", data.nvertex);
return
end

