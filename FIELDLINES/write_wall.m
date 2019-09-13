function data = write_wall(data,filename )
%WRITE_WALL Read's a triangulated mesh.
%   WRITE_WALL writes a wall part data file.  The wall part data file is a
%   file with the format:
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
%   Example:
%       lim_data=read_wall('limiter_trimesh.dat');
%       write_wall(lim_data,'test_file.dat');
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           9/29/16
% Check arguments
if nargin<2
    disp('ERROR: read_vessel requires wall structure and filename');
    return
end
% Open File
fid=fopen(filename,'w');
% Write Header
fprintf(fid,'MACHINE: %s\n',data.machine);
fprintf(fid,'DATE: %s\n',data.date);
fprintf(fid,'%i %i\n',data.nvertex,data.nfaces);
% Write Coords
fprintf(fid,'%20.10E %20.10E %20.10E\n',data.coords);
% Write Faces
fprintf(fid,'%i %i %i\n',data.faces);
% Close file
fclose(fid);
return;

end

