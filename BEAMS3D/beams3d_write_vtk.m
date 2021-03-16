function beams3d_write_vtk(beam_data,quantity)
%BEAMS3D_WRITE_VTK Creates a VTK legacy text file from BEAMS3D data
%   The BEAMS3D_WRITE_VTK routine takes a BEAMS3D data strcuture as read by
%   the READ_BEAMS3D routine and a string specifying a quantity as
%   arguments.  It outputs a text file which VTK can read in legacy format.
%
%   Quantities
%       orbits:     Outputs the orbits of the particles as stored in _lines
%       beam:       Optputs the beam deposition trace
%       wall_loss:  Not implemented
%       wall_power: Not implemented
%
% Example usage
%      data=read_beams3d('beams3d_test.h5');  % Reads BEAMS3D HDF5 file
%      beams3d_write_vtk(data,'orbits');
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       0.1

switch quantity
    case {'orbits'}
        filename='beams3d_vtk_orbits.vtk';
        fid = fopen(filename,'w');
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,['BEAMS3D DATA by beams3d_write_vtk :' quantity '\n']);
        fprintf(fid,'ASCII\n');
        fprintf(fid,'DATASET POLYDATA\n');
        n = size(beam_data.X_lines,1)-1;
        n2 = size(beam_data.X_lines,2);
        fprintf(fid,'POINTS %i float\n',n*n2);
        for i=1:n2
            out = [beam_data.X_lines(1:n,i) beam_data.Y_lines(1:n,i) beam_data.Z_lines(1:n,i)]';
            fprintf(fid,'%7.4f %7.4f %7.4f\n',out);
        end
        fprintf(fid,'LINES %i %i\n',n2, (n+1)*n2);
        d = (1:n)-1;
        for i=1:size(beam_data.X_lines,2)
            fprintf(fid,' %i ',[n d]);
            fprintf(fid,'\n');
            d = d+n;
        end
        fclose(fid);
    case {'beam'}
        filename='beams3d_vtk_beam.vtk';
        fid = fopen(filename,'w');
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,['BEAMS3D DATA by beams3d_write_vtk :' quantity '\n']);
        fprintf(fid,'ASCII\n');
        fprintf(fid,'DATASET POLYDATA\n');
        n = 2;
        n2 = size(beam_data.X_lines,2);
        fprintf(fid,'POINTS %i float\n',n*n2);
        for i=1:n2
            out = [beam_data.X_lines(1:n,i) beam_data.Y_lines(1:n,i) beam_data.Z_lines(1:n,i)]';
            fprintf(fid,'%7.4f %7.4f %7.4f\n',out);
        end
        fprintf(fid,'LINES %i %i\n',n2, (n+1)*n2);
        d = (1:n)-1;
        for i=1:size(beam_data.X_lines,2)
            fprintf(fid,' %i ',[n d]);
            fprintf(fid,'\n');
            d = d+n;
        end
        fprintf(fid,'POINT_DATA %i\n',n*n2);
        fprintf(fid,'SCALARS sample_scalars float 1\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        temp = double([beam_data.end_state beam_data.end_state]);
        fprintf(fid,'%7.4f\n',temp(:));
        fclose(fid);
        
end


end

