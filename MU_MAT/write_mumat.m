function write_mumat(data,filename)
%WRITE_MUMAT Write as a magnetic material tetrahedral mesh
%   WRITE_MUMAT writes a magnetic material tetrahedral mesh in the
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
%       write_mumag(Example,'sphere2.dat');
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           9/12/23

% Check arguments
if nargin<2
    disp('ERROR: write_mumat requires wall structure and filename');
    return
end
% Open File
fid=fopen(filename,'w');
% Write Header
fprintf(fid,'MACHINE: %s\n',data.machine);
fprintf(fid,'DATE: %s\n',data.date);
fprintf(fid,'%i %i %i\n',data.nvertex,data.ntet,data.nstate);
% Write Coords
fprintf(fid,'%20.10E %20.10E %20.10E\n',data.coords);
% Write Faces
fprintf(fid,'%i %i %i %i %i\n',[data.tet; data.func_dex]);
% Write functions
for i = 1:data.nstate
    fprintf(fid,'%i\n',data.state_func(i).type);
    if data.state_func(i).type == 1 %hard
        fprintf(fid,'%20.10E %20.10E %20.10E\n',data.state_func(i).M);
    elseif data.state_func(i).type == 2 % Soft
        fprintf(fid,'%i %i\n',data.state_func(i).nh,data.state_func(i).nt);
        fprintf(fid,' %20.10E',data.state_func(i).H);
        fprintf(fid,'\n');
        fprintf(fid,' %20.10E',data.state_func(i).M);
        fprintf(fid,'\n');
    elseif data.state_func(i).type == 3 % Soft constant
        fprintf(fid,'%20.10E\n',data.state_func(i).mu);
    end
% Close file
fclose(fid);
return;

end