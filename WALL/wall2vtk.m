function wall2vtk(wall_data)
%WALL2VTK Routine outputs wall data on VTK format.
%   The WALL2VTK routine outputs wall data in the VTK format.  The routine
%   takes a WALL data structure as returned by READ_WALL as an argument.
%
% Example usage
%      wall_data=read_wall('wall_element.dat'); %.dat or .stl
%      wall2vtk(wall_data);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0


filename='wall_data.vtk';
fid = fopen(filename,'w');
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,['WALL DATA by wall2vtk\n']);
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET POLYDATA\n');
n = wall_data.nvertex;
n2 = wall_data.nfaces;
fprintf(fid,'POINTS %i float\n',n);
fprintf(fid,'%7.4f %7.4f %7.4f\n',wall_data.coords);
fprintf(fid,'POLYGONS %i %i\n',n2,4*n2);
fprintf(fid,'3 %i %i %i\n',wall_data.faces-1);
fclose(fid);

end

