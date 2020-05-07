function ascot5_writewall(filename,wall_data,varargin)
%ASCOT5_WRITEWALL Writes a wall structure to an ASCOT5 style HDF5 file.
%   The ASCOT5_WRITEWALL function writes a wall structure as read by
%   read_wall to an HDF5 file in a format which is readable by the ASCOT5
%   code.  A description is automatically created if none is provided.  If
%   the file exists the /wall/wall_3D_ dataset is added to the file and set
%   as active.
%
%   Example:
%       wall_data = read_wall('test.dat');
%       ascot5_writewall('ascot5_test.h5',wall_data); % Auto Description
%       desc_text = 'This is a test wall, have fun.';
%       ascot5_writewall('ascot5_test.h5',wall_data,'description',desc_text);
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00

% defaults
desc_text=[];

% Handle varargin
if nargin > 2
    i = 1;
    while (i <= nargin-2)
        switch varargin{i}
            case{'desciption'}
                i=i+1;
                desc_text=varargin{i};
        end
        i=i+1;
    end
end

% Check wall
if isstruct(wall_data)
    if isfield(wall_data,'datatype')
        if ~strcmp(wall_data.datatype,'limiter_trimesh')
            disp(['  ERROR: wall_data not a limiter_trimesh strucutre.']);
            return;
        end
    else
        disp(['  ERROR: wall_data missing datatype field.']);
        return;
    end
else
    disp(['  ERROR: wall_data not a structure.']);
    return;
end


% Default descirption if not provided
if isempty(desc_text)
    desc_text = ['Generated in Matlab from limiter_trimesh: ' wall_data.machine];
end

% Divide up into ASCOT5 format
d1 = wall_data.faces(1,:);
d2 = wall_data.faces(2,:);
d3 = wall_data.faces(3,:);
x123 = [wall_data.coords(1,d1); wall_data.coords(1,d2); wall_data.coords(1,d3)];
y123 = [wall_data.coords(2,d1); wall_data.coords(2,d2); wall_data.coords(2,d3)];
z123 = [wall_data.coords(3,d1); wall_data.coords(3,d2); wall_data.coords(3,d3)];

% Create random id string
id = num2str(round(rand*1E10),'%10.10i');

% Handle existing file
if isfile(filename)
    disp(['  ' filename ' exists, adding wall_3D element ' id ' to file.']);
end

% Create Dataset
h5create(filename,['/wall/wall_3D_' id '/nelements'],1,'Datatype','int32');
h5create(filename,['/wall/wall_3D_' id '/x1x2x3'],size(x123));
h5create(filename,['/wall/wall_3D_' id '/y1y2y3'],size(y123));
h5create(filename,['/wall/wall_3D_' id '/z1z2z3'],size(z123));

% Write Attributes
h5writeatt(filename,['/wall'],'active',id);
h5writeatt(filename,['/wall/wall_3D_' id ],'date',datestr(now,'yyyy-mm-dd hh:MM:ss'));
h5writeatt(filename,['/wall/wall_3D_' id ],'description',desc_text);


% Write Variables
h5write(filename,['/wall/wall_3D_' id '/nelements'],int32(wall_data.nfaces));
h5write(filename,['/wall/wall_3D_' id '/x1x2x3'],x123);
h5write(filename,['/wall/wall_3D_' id '/y1y2y3'],y123);
h5write(filename,['/wall/wall_3D_' id '/z1z2z3'],z123);

return;

end



