function data = read_mumat(filename)
%READ_MUMAT Reads a magnetic materias metrahedral mesh file
%   READ_MUMAT reads a magnetic material tetrahedral mesh in the
%   following format:
%
%MACHINE:  W-7X_limiter
%DATE:  07-01-14
%07612 14399 1
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
%    . . .
%    . . .
%   5
%   1.0 2.0 3.0 4.0 5.0
%   0.0 0.0 0.0 0.0 0.0
%
%   Here the string after MACHINE is the name of the machine.  The DATE
%   string contains information about the date of the file.  The next line
%   indicates the number of verticies, number of tetrahedral cells and
%   the number of state functions. What follows next are the vertices in
%   x,y,z [m] format.  Then comes the array of indices defining each
%   tetrahedral cell.  The first 4 are indices, the last index references
%   the state function index.  Finally, the state functions are listed.
%   First comes the number of points in the function.  For now
%   state-functions are not supported.  If the number is -1 then it is
%   assumed a hard magnet and the magnetization is listed next. If -2 then
%   is is a soft magnet with constant value.
%
%   Example:
%       sphere_data=read_mumag('sphere.dat');
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           9/12/23

% Defaults
data=[];

% Check arguments
if nargin<1
    disp('ERROR: read_mumat requires filename');
    return
end

% Open File
fid=fopen(filename,'r');
% Read Header
header_line1=fgetl(fid);
header_line2=fgetl(fid);
data.machine=strtrim(header_line1(strfind(header_line1,':')+1:numel(header_line1)));
data.date=strtrim(header_line2(strfind(header_line2,':')+1:numel(header_line2)));
temp=fscanf(fid,'%d',3);
data.nvertex = temp(1);
data.ntet = temp(2);
data.nstate = temp(3);
data.coords=fscanf(fid,'%E %E %E',[3 data.nvertex]);
temp=fscanf(fid,'%d %d %d %d',[5 data.ntet]);
data.tet = temp(1:4,:);
data.func_dex = temp(5,:);
for i = 1:data.nstate
    data.state_func(i).type = fscanf(fid,'%i',1);
    if data.state_func(i).type == 1 %hard
        data.state_func(i).M = fscanf(fid,'%E',3);
    elseif data.state_func(i).type == 2 %Soft
        temp = fscanf(fid,'%i',2);
        data.state_func(i).nh = temp(1);
        data.state_func(i).nt = temp(2);
        data.state_func(i).H = fscanf(fid,'%E',data.state_func(i).nt);
        data.state_func(i).M = fscanf(fid,'%E',data.state_func(i).nh);
    elseif data.state_func(i).type == 3 %Soft Constant
        data.state_func(i).mu = fscanf(fid,'%E',1);
end

fclose(fid);
return;

end