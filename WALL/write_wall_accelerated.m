function data = write_wall_accelerated(data,filename )
%WRITE_WALL_ACCELERATED Creates a triangulated mesh of an accelerated wall mesh.
%   WRITE_WALL_ACCELERATED writes a wall part data file.  
%   The wall part data file is a file with the format:
%
%MACHINE:  W-7X_limiter
%DATE:  07-01-14
%0      0   
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
%   is 0 0 to differentiate between an unnaccelerated mesh. The next line
%   indicates the number of verticies and number of triangular cells.  What
%   follows are a list of vertices (x,y,z) and then the blocks with their
%   info and a list of points composing each cell in the block.
%
%   Example:
%       lim_data=read_wall('limiter_trimesh.dat');
%       write_wall(lim_data,'test_file.dat');
%
%   Written by:     D.J. Engels (d.j.engels@student.tue.nl)
%   Version:        1.0
%   Date:           05/21
% Check arguments
if nargin<2
    disp('ERROR: write_wall_accelerated requires wall structure and filename');
    return
end

% Open File
fid=fopen(filename,'w');
% Write Header
fprintf(fid,'MACHINE: %s\n',data.machine);
fprintf(fid,'DATE: %s\n',data.date);
fprintf(fid,'%i %i\n', 0, 0);
fprintf(fid,'%i %i\n',data.nvertex,data.nfaces);
% Write Coords
fprintf(fid,'%20.10E %20.10E %20.10E\n',data.coords);
% Write Faces
fprintf(fid,'%i %i %i\n',data.faces);
% Write each block with its faces
fprintf(fid, '%i %i %i %i\n', data.nblocks, data.xstep, data.ystep, data.zstep);
fprintf(fid, '%20.10E %i %i %i\n', data.block_size, data.nblocks_x, data.nblocks_y, data.nblocks_z);
for i=1:data.nblocks
    tmp = data.blocks(i);
    fprintf(fid, '%20.10E %20.10E %20.10E %20.10E %20.10E %20.10E\n', tmp.x_min, tmp.x_max, tmp.y_min, tmp.y_max, tmp.z_min, tmp.z_max);
    fprintf(fid, '%i\n', tmp.nfaces);
    if tmp.nfaces > 0
        fprintf(fid,'%i\n',tmp.faces);
    end
end
% Close file
fclose(fid);
return;

end