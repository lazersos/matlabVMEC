function data = read_wall( filename )
%READ_WALL Read's a triangulated mesh.
%   READ_WALL reads a wall part data file.  The wall part data file is a
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
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           9/29/16
% Check arguments
if nargin<1
    disp('ERROR: read_vessel requires filename');
    return
end
% Open File
fid=fopen(filename,'r');
% Read Header
header_line1=fgetl(fid);
header_line2=fgetl(fid);
data.machine=strtrim(header_line1(strfind(header_line1,':')+1:numel(header_line1)));
data.date=strtrim(header_line2(strfind(header_line2,':')+1:numel(header_line2)));
temp=fscanf(fid,'%d',2);
data.nvertex = temp(1);
data.nfaces = temp(2);
% Read dataset
data.coords=fscanf(fid,'%E %E %E',[3 data.nvertex]);
data.faces=fscanf(fid,'%d %d %d',[3 data.nfaces]);
fclose(fid);
% Fix if 0 index
if min(min(data.faces)) == 0
    data.faces = data.faces+1;
end
% Prep the wall
vertex = data.coords;
face   = data.faces;
dex1 = face(1,:);
dex2 = face(2,:);
dex3 = face(3,:);
data.A  = vertex(:,dex1);
data.V0 = vertex(:,dex3)-vertex(:,dex1);
data.V1 = vertex(:,dex2)-vertex(:,dex1);
data.FN = cross(data.V1,data.V0);
data.DOT00=dot(data.V0,data.V0);
data.DOT01=dot(data.V0,data.V1);
data.DOT11=dot(data.V1,data.V1);
data.invDenom = 1.0./(data.DOT00.*data.DOT11-data.DOT01.*data.DOT01);
data.d = dot(data.FN,data.A);
data.datatype='limiter_trimesh';
return;

end

