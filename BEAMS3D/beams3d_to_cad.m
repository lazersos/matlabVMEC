function beams3d_to_cad(beam_data,filename)
%BEAMS3D_TO_CAD Outputs wall heatflux data in ASCII format for CAD
%   The BEAMS3D_TO_CAD subroutine outputs the BEAMS3D wall load data to a
%   ASCII text file (.asc) where each line is composed of:
%   X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3 HEATFLUX
%   where there is a line for each 'face' in the model.  The default
%   behavior is to convert from m to mm and output the heatflux in units of
%   W/m^2.
%
% Example usage
%      beam_data=read_beams3d('beams3d_test.h5');
%      beams3d_to_cad(beam_data,'beams3d_heatflux_test.asc');
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

factor=1000; %m to mm

faces = beam_data.wall_faces;
vertex = beam_data.wall_vertex.*factor;
heatflux = sum(beam_data.wall_load,1);

heatflux(heatflux>99E6) = 99E6; % limit to 99 MW/m^2

x1 = vertex(faces(:,1),1);
y1 = vertex(faces(:,1),2);
z1 = vertex(faces(:,1),3);
x2 = vertex(faces(:,2),1);
y2 = vertex(faces(:,2),2);
z2 = vertex(faces(:,2),3);
x3 = vertex(faces(:,3),1);
y3 = vertex(faces(:,3),2);
z3 = vertex(faces(:,3),3);

outdata=[x1,y1,z1,x2,y2,z2,x3,y3,z3,heatflux'];
fid = fopen(filename,'w');
fprintf(fid,'%+7.1f\t%+7.1f\t%+7.1f\t%+7.1f\t%+7.1f\t%+7.1f\t%+7.1f\t%+7.1f\t%+7.1f\t%+9.1f\n',outdata');
fclose(fid);
end