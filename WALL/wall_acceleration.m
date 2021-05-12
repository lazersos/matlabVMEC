function wall_accelerated = wall_acceleration( wall, block_size )
%WALL_ACCELERATION Create an uniform grid acceleration structure for a wall mesh
%   WALL_ACCELERATION reads wall data or takes it from input
%   Wall data has to be input using READ_WALL or in file structure
%   that can be read by READ_WALL
%   Example of file format:
%
%MACHINE:  W-7X_limiter
%DATE:  07-01-14
%07612 14399
%    5.5057190000E+00    -7.4592000000E-02    -4.4169700000E-01
%    5.5057190000E+00     7.4592000000E-02     4.4169700000E-01
%    5.5063620000E+00    -7.2601000000E-02    -4.4126700000E-01
%    5.5063620000E+00     7.2601000000E-02     4.4126700000E-01
%    5.5165410000E+00    -7.5732000000E-02    -4.1455700000E-01
%    5.5165410000E+00     7.5732000000E-02     4.1455700000E-01
%         .                      .                     .
%         .                      .                     .
%         .                      .                     .
%    1 2 3
%    4 5 6
%    . . .
%
%   Here the string after MACHINE is the name of the machine.  The DATE
%   string contains information about the date of the file.  The next line
%   indicates the number of verticies and number of triangular cells.  What
%   follows are a list of vertices (x,y,z) then the list of points
%   composing each cell.
%
%   Written by:     D.J. Engels (d.j.engels@student.tue.nl)
%   Version:        1.0
%   Date:           05/21
if nargin<2
    disp('ERROR: wall_acceleration requires wall structure / filename and a block size');
    return
end
%% check if wall given or location for wall
if isa(wall, 'char')
    wall = read_wall(wall);
end

epsilon = 1e-6;
%% Create uniform grid for mesh
% find maxima
x_max = max(wall.coords(1, :));
x_min = min(wall.coords(1, :));
y_max = max(wall.coords(2, :));
y_min = min(wall.coords(2, :));
z_max = max(wall.coords(3, :));
z_min = min(wall.coords(3, :));

% calculate starting points of blocks using uniform buffer
x_buf = ceil((x_max - x_min)/ block_size) * block_size - (x_max - x_min);
xs = x_min - x_buf/2:block_size:x_max+x_buf/2-epsilon;
y_buf = ceil((y_max - y_min)/ block_size) * block_size - (y_max - y_min);
ys = y_min - y_buf/2:block_size:y_max+y_buf/2-epsilon;
z_buf = ceil((z_max - z_min)/ block_size) * block_size - (z_max - z_min);
zs = z_min - z_buf/2:block_size:z_max+z_buf/2-epsilon;

% create blocks by looping over all x's, y's, and z's
grid_blocks(length(xs)*length(ys)*length(zs)) = struct();
counter=1;
for xi=1:length(xs)
   for yi=1:length(ys)
      for zi=1:length(zs)
          x = xs(xi);
          y = ys(yi);
          z = zs(zi);
          grid_blocks(counter).x_min = x;
          grid_blocks(counter).x_max = x + block_size;
          grid_blocks(counter).y_min = y;
          grid_blocks(counter).y_max = y + block_size;
          grid_blocks(counter).z_min = z;
          grid_blocks(counter).z_max = z + block_size;
          grid_blocks(counter).nfaces = 0;
          grid_blocks(counter).faces = [];
          counter = counter + 1;
      end
   end 
end

%% assign triangles to grid
parfor i=1:length(grid_blocks)
    tmp = grid_blocks(i);
    xmax = tmp.x_max + epsilon;
    ymax = tmp.y_max + epsilon;
    zmax = tmp.z_max + epsilon;
    xmin = tmp.x_min - epsilon;
    ymin = tmp.y_min - epsilon;
    zmin = tmp.z_min - epsilon;
    % check which vertices in this block
    mask = wall.coords(1, :) < xmax & wall.coords(1, :) >= xmin ...
        & wall.coords(2, :) < ymax & wall.coords(2, :) >= ymin ...
        & wall.coords(3, :) < zmax & wall.coords(3, :) >= zmin;
    % get those vertex indexes
    vertices = int32(1:wall.nvertex);
    vertices = vertices(mask);
    
    % check how many faces in block
    counter = int32(CountOccurance(vertices, int32(wall.faces)));
    
    % Get faces
    faces = MakeFaces(vertices, int32(wall.faces), counter);
    % store the unique faces, since you could get double from method above
    if ~isempty(faces) > 0
        % only store the indices of the faces
        grid_blocks(i).faces = unique(faces(4,:)', 'rows')';
        grid_blocks(i).nfaces = size(grid_blocks(i).faces,2);
    end
end
%% finish up new wall
wall_accelerated.machine = wall.machine;
try
    wall_accelerated.date = datestr(today());
catch
    fprintf("You need the financial toolbox for automatic date setting\n");
    wall_accelerated.date = wall.date;
end
wall_accelerated.nvertex = wall.nvertex;
wall_accelerated.nfaces = wall.nfaces;
wall_accelerated.coords = wall.coords;
wall_accelerated.faces = wall.faces;
wall_accelerated.nblocks = length(grid_blocks);
wall_accelerated.nblocks_x = xi;
wall_accelerated.nblocks_y = yi;
wall_accelerated.nblocks_z = zi;
wall_accelerated.block_size = block_size;
wall_accelerated.blocks = grid_blocks;
wall_accelerated.datatype = wall.datatype;
wall_accelerated.zstep = 1;
wall_accelerated.ystep = length(zs);
wall_accelerated.xstep = length(ys)*length(zs);
%% check number of faces
wall_accelerated.nfaces_blocks = sum([wall_accelerated.blocks.nfaces]);
if wall_accelerated.nfaces_blocks ~= wall.nfaces
   fprintf("More faces created by accelerated wall due to faces being in multiple blocks. Not an issue, just a warning\n");
   fprintf("New number of faces: %d\n", wall_accelerated.nfaces_blocks);
   fprintf("Old number of faces: %d\n", wall_accelerated.nfaces);
end
end