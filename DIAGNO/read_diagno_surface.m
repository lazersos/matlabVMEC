function surf_data=read_diagno_surface(filename)
%READ_DIAGNO_SURFACE(filename) Reads the DIAGNO surface file.
%   This function reads the DIAGNO surface.geom or axis.geom file and 
%   returns the contents as a data structure.  This structure has the
%   following elements:
%
%   surf_data:
%       data:       Array of values.
%       header:     Header of file.
%       version:    DIAGNO Version (based on number of elements)
%       datatype:   'diagno_surf' or 'diagno_axis'
%
%   The code returns -1 if an error is found.
%
%   Example:
%       dsurf_data=read_diagno_surface('plasma_surface.geom');
%
%   See also plot_diagno_surface.
%
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        2.0
%   Date:           3/20/14

% Check to see if it's an plasma_axis or plasma_surface file
if strcmp(filename,'plasma_axis.geom')
    surf_data.datatype='diagno_axis';
    % Read file
    fid=fopen(filename);
    line=fgetl(fid);                              %Read First Line
    surf_data.data=sscanf(line,'%e');
    while ~feof(fid)
        line=fgetl(fid);
        if ~strcmp(line,'###');
            surf_data.data=[surf_data.data sscanf(line,'%e')];
        end
    end
    fclose(fid);
    surf_data.version=1.0;
    surf_data.header='Element x y z';
elseif strcmp(filename,'plasma_surface.geom')
    surf_data.datatype='diagno_surf';
    % Read file
    fid=fopen(filename);
    surf_data.header=fgetl(fid);                %Read Header
    line=fgetl(fid);                              %Read First Line
    surf_data.data=sscanf(line,'%e');
    while ~feof(fid)
        line=fgetl(fid);
        if ~strcmp(line,'###');
            surf_data.data=[surf_data.data sscanf(line,'%e')];
        end
    end
    fclose(fid);
    % We should add some version information
    switch size(surf_data.data,1)
        case 13
            surf_data.version=1.0;
        case 16
            surf_data.version=1.5;
        case 12
            surf_data.version=2.0;
        otherwise
            disp(['ERROR: Unknown version, numels=' num2str(size(surf_data.data,2)) '!']);
            surf_data.version=-1.0;
    end
else
    try
        fid=fopen(filename);
        surf_data.header=fgetl(fid);                %Read Header
        surf_data.data = fscanf(fid,'%e',[14 inf]);
        fclose(fid);
        surf_data.datatype='diagno_surf';
    catch
        disp(['ERROR: Unknown filetype ' filename]);
        surf_data=-1;
    end
end
return
end