function ascot5_writemarker_ml(filename,r,phi,z,pitch,weight,time,idnum)
%ASCOT5_WRITEMARKER_ML Writes marker_ml data to an HDF5 file for ASCOT5
%   The ASCOT5_WRITEMARKER_ML routine writes marker_ml data to an HDF5 file
%   in the style read by the ASCOT5 code.  If the file exists the
%   /marker/fl_ dataset is added to the file and set as active.
%
%   Example:
%       r      = [5.5 6.5 7.5]; % [m]
%       phi    = [0.0 0.0 0.0]; % [deg]
%       z      = [0.0 0.0 0.0]; % [m]
%       pitch  = [1.0 1.0 1.0];   % [+/- direction]
%       weight = [1.0 1.0 1.0];   % [particles/particle]
%       time   = [0.0 0.0 0.0];   % [s]
%       id     = [  1   2   3];   % [none]
%       ascot5_writemarker_fl('test.h5',r,phi,z,pitch,weight,time,id);
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00

%Constants

% Create random id string
id = num2str(round(rand*1E10),'%10.10i');

% Handle existing file
if isfile(filename)
    disp(['  ' filename ' exists, adding /marker/fl_' id ' to file.']);
end

% Create datasets
h5create(filename,['/marker/fl_' id '/n'],1,'Datatype','int32');
h5create(filename,['/marker/fl_' id '/r'],size(r));
h5create(filename,['/marker/fl_' id '/phi'],size(phi));
h5create(filename,['/marker/fl_' id '/z'],size(z));
h5create(filename,['/marker/fl_' id '/pitch'],size(pitch));
h5create(filename,['/marker/fl_' id '/weight'],size(weight));
h5create(filename,['/marker/fl_' id '/time'],size(time));
h5create(filename,['/marker/fl_' id '/id'],size(idnum),'Datatype','int32');


% Write Attributes
h5writeatt(filename,'/marker','active',id,'TextEncoding','system');
h5writeatt(filename,['/marker/fl_' id ],'date',datestr(now,'yyyy-mm-dd hh:MM:ss'),'TextEncoding','system');
h5writeatt(filename,['/marker/fl_' id ],'description','Written by MATLAB','TextEncoding','system');

% Write Variables
h5write(filename,['/marker/fl_' id '/n'],int32(length(r)));
h5write(filename,['/marker/fl_' id '/r'],r);
h5write(filename,['/marker/fl_' id '/phi'],phi);
h5write(filename,['/marker/fl_' id '/z'],z);
h5write(filename,['/marker/fl_' id '/pitch'],pitch);
h5write(filename,['/marker/fl_' id '/weight'],weight);
h5write(filename,['/marker/fl_' id '/time'],time);
h5write(filename,['/marker/fl_' id '/id'],int32(idnum));


end

