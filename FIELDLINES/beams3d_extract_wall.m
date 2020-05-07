function wall_data = beams3d_extract_wall(beam_data)
%BEAMS3D_EXTRACT_WALL Extracts a wall structure from a BEAMS3D structure
%   The BEAMS3D_EXTRACT_WALL function is design to extract a wall strucutre
%   as returned by read_wall from a BEAMS3D data structure as read by
%   read_beams3d.
%
%   Example:
%       beam_data = read_beams3d('test.h5');
%       wall_data = beams3d_extract_wall(beam_data);
%       plot_wall(wall_data,'solid');
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           5/05/20

% Defaults
wall_data=[];

% check input
if isstruct(beam_data)
    if isfield(beam_data,'datatype')
        if ~strcmp(beam_data.datatype,'BEAMS3D')
            disp('Error: beam_data not a BEAMS3D data strucutre');
            return;
        end
    else
        disp('Error: beam_data does not posses a datatype field.');
        return;
    end
    if ~isfield(beam_data,'wall_faces')
        disp('Error: wall_faces not found in beam_data.');
        return;
    end
    if ~isfield(beam_data,'wall_vertex')
        disp('Error: wall_vertex not found in beam_data.');
        return;
    end
else
    disp('Error: beam_data not a structure.');
    return;
end

% Extract Wall
wall_data.machine=['BEAMS3D ' num2str(beam_data.VERSION,'v%3.2f') ' wall'];
wall_data.date=datestr(now,'mm-dd-yyyy');
wall_data.coords=beam_data.wall_vertex';
wall_data.faces=beam_data.wall_faces';
wall_data.nvertex = size(wall_data.coords,2);
wall_data.nfaces = size(wall_data.faces,2);
wall_data.datatype='limiter_trimesh';



end

