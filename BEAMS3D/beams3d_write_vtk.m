function beams3d_write_vtk(beam_data,quantity)
%BEAMS3D_WRITE_VTK Creates a VTK legacy text file from BEAMS3D data
%   The BEAMS3D_WRITE_VTK routine takes a BEAMS3D data strcuture as read by
%   the READ_BEAMS3D routine and a string specifying a quantity as
%   arguments.  It outputs a text file which VTK can read in legacy format.
%
%   Quantities
%       orbits:     Outputs the orbits of the particles as stored in _lines
%       beam:       Optputs the beam deposition trace
%       wall:       Basic Wall Output
%       wall_loss:  Lost particle depostion pattern
%       wall_power: Wall power deposition
%       wall_shine: Shinethrough wall power
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
        fprintf(fid,'POINT_DATA %i\n',n*n2);
        fprintf(fid,'SCALARS mod(B) float 1\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        for i=1:n2
            fprintf(fid,'%7.4f\n',beam_data.B_lines(1:n,i));
        end
        fclose(fid);
    case {'wall'}
        filename='beams3d_vtk_wall.vtk';
        fid = fopen(filename,'w');
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,['BEAMS3D DATA by beams3d_write_vtk :' quantity '\n']);
        fprintf(fid,'ASCII\n');
        fprintf(fid,'DATASET POLYDATA\n');
        n = size(beam_data.wall_vertex,1);
        n2 = size(beam_data.wall_faces,1);
        fprintf(fid,'POINTS %i float\n',n);
        fprintf(fid,'%7.4f %7.4f %7.4f\n',beam_data.wall_vertex');
        fprintf(fid,'POLYGONS %i %i\n',n2,4*n2);
        fprintf(fid,'3 %i %i %i\n',beam_data.wall_faces'-1);
        fclose(fid);
    case {'wall_loss'}
        filename='beams3d_vtk_wall_loss.vtk';
        fid = fopen(filename,'w');
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,['BEAMS3D DATA by beams3d_write_vtk :' quantity '\n']);
        fprintf(fid,'ASCII\n');
        fprintf(fid,'DATASET POLYDATA\n');
        n = size(beam_data.wall_vertex,1);
        n2 = size(beam_data.wall_faces,1);
        fprintf(fid,'POINTS %i float\n',n);
        fprintf(fid,'%7.4f %7.4f %7.4f\n',beam_data.wall_vertex');
        fprintf(fid,'POLYGONS %i %i\n',n2,4*n2);
        fprintf(fid,'3 %i %i %i\n',beam_data.wall_faces'-1);
        fprintf(fid,'CELL_DATA %i\n',n2);
        fprintf(fid,'SCALARS wall_hits float 1\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        fprintf(fid,'%7.4f\n',beam_data.wall_strikes);
        fclose(fid);
    case {'wall_power'}
        filename='beams3d_vtk_wall_load.vtk';
        fid = fopen(filename,'w');
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,['BEAMS3D DATA by beams3d_write_vtk :' quantity '\n']);
        fprintf(fid,'ASCII\n');
        fprintf(fid,'DATASET POLYDATA\n');
        n = size(beam_data.wall_vertex,1);
        n2 = size(beam_data.wall_faces,1);
        fprintf(fid,'POINTS %i float\n',n);
        fprintf(fid,'%7.4f %7.4f %7.4f\n',beam_data.wall_vertex');
        fprintf(fid,'POLYGONS %i %i\n',n2,4*n2);
        fprintf(fid,'3 %i %i %i\n',beam_data.wall_faces'-1);
        fprintf(fid,'CELL_DATA %i\n',n2);
        fprintf(fid,'SCALARS Heatflux float 1\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        fprintf(fid,'%7.4f\n',sum(beam_data.wall_load));
        fclose(fid);
    case {'wall_shine'}
        filename='beams3d_vtk_wall_shine.vtk';
        fid = fopen(filename,'w');
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,['BEAMS3D DATA by beams3d_write_vtk :' quantity '\n']);
        fprintf(fid,'ASCII\n');
        fprintf(fid,'DATASET POLYDATA\n');
        n = size(beam_data.wall_vertex,1);
        n2 = size(beam_data.wall_faces,1);
        fprintf(fid,'POINTS %i float\n',n);
        fprintf(fid,'%7.4f %7.4f %7.4f\n',beam_data.wall_vertex');
        fprintf(fid,'POLYGONS %i %i\n',n2,4*n2);
        fprintf(fid,'3 %i %i %i\n',beam_data.wall_faces'-1);
        fprintf(fid,'CELL_DATA %i\n',n2);
        fprintf(fid,'SCALARS Heatflux float 1\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        fprintf(fid,'%7.4f\n',sum(beam_data.wall_shine));
        fclose(fid);
    case {'beam'}
        n = 2;
        filenames={'depo','therm','wall','shine','port'};
        for j=[0 3 4]
            dex = beam_data.end_state==j;
            n2 = sum(dex);
            if n2>0
                filename=['beams3d_vtk_beam_' filenames{j+1} '.vtk'];
                fid = fopen(filename,'w');
                fprintf(fid,'# vtk DataFile Version 2.0\n');
                fprintf(fid,['BEAMS3D DATA by beams3d_write_vtk :' quantity ' port\n']);
                fprintf(fid,'ASCII\n');
                fprintf(fid,'DATASET POLYDATA\n');
                fprintf(fid,'POINTS %i float\n',n*n2);
                x = beam_data.X_lines(1:n,dex);
                y = beam_data.Y_lines(1:n,dex);
                z = beam_data.Z_lines(1:n,dex);
                for i=1:n2
                    out = [x(1:n,i) y(1:n,i) z(1:n,i)]';
                    fprintf(fid,'%7.4f %7.4f %7.4f\n',out);
                end
                fprintf(fid,'LINES %i %i\n',n2, (n+1)*n2);
                d = (1:n)-1;
                for i=1:n2
                    fprintf(fid,' %i ',[n d]);
                    fprintf(fid,'\n');
                    d = d+n;
                end
                fprintf(fid,'POINT_DATA %i\n',n*n2);
                fprintf(fid,'SCALARS BEAM_INDEX float 1\n');
                fprintf(fid,'LOOKUP_TABLE default\n');
                temp = double([beam_data.Beam(dex) beam_data.Beam(dex)]);
                fprintf(fid,'%7.4f\n',temp(:));
                fclose(fid);
            end
        end
        
end


end

