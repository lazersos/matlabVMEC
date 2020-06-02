function ascot5_writemarker_gc(filename,r,phi,z,energy,pitch,zeta,mass,charge,anum,znum,weight,time,idnum)
%ASCOT5_WRITEMARKER_GC Writes marker_gc data to an HDF5 file for ASCOT5
%   The ASCOT5_WRITEMARKER_GC routine writes marker_gc data to an HDF5 file
%   in the style read by the ASCOT5 code.  If the file exists the
%   /marker/marker_gc_ dataset is added to the file and set as active.
%
%   Example:
%       r      = [5.5 6.5 7.5]; % [m]
%       phi    = [0.0 0.0 0.0]; % [deg]
%       z      = [0.0 0.0 0.0]; % [m]
%       energy = [55E3 55E3 55E3]; % [eV]
%       pitch  = [0.1 0.1 0.1];   % [Vll/Vtotal]
%       zeta   = [0.0 0.0 0.0];   % [???]
%       mass   = [1.0 1.0 1.0];   % [amu]
%       charge = [  1   1   1];   % [ec]
%       anum   = [  1   1   1];   % A [norm]
%       znum   = [  1   1   1];   % Z [norm]
%       weight = [1.0 1.0 1.0];   % [particles/particle]
%       time   = [0.0 0.0 0.0];   % [s]
%       id     = [  1   2   3];   % [none]
%       ascot5_writemarker_gc('test.h5',r,phi,z,energy,pitch,zeta,mass,...
%                               charge,anum,znum,weight,time,id);
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00

%Constants
amu = 1.66053906660E-27;
ec  = 1.60217662E-19;
cspeed = 2.99792458E+08;

% Create random id string
id = num2str(round(rand*1E10),'%10.10i');

% Handle existing file
if isfile(filename)
    disp(['  ' filename ' exists, adding /marker/gc_' id ' to file.']);
end

% Create datasets
h5create(filename,['/marker/gc_' id '/n'],1,'Datatype','int32');
h5create(filename,['/marker/gc_' id '/r'],size(r));
h5create(filename,['/marker/gc_' id '/phi'],size(phi));
h5create(filename,['/marker/gc_' id '/z'],size(z));
h5create(filename,['/marker/gc_' id '/energy'],size(energy));
h5create(filename,['/marker/gc_' id '/pitch'],size(pitch));
h5create(filename,['/marker/gc_' id '/zeta'],size(zeta));
h5create(filename,['/marker/gc_' id '/mass'],size(mass));
h5create(filename,['/marker/gc_' id '/charge'],size(charge),'Datatype','int32');
h5create(filename,['/marker/gc_' id '/anum'],size(anum),'Datatype','int32');
h5create(filename,['/marker/gc_' id '/znum'],size(znum),'Datatype','int32');
h5create(filename,['/marker/gc_' id '/weight'],size(weight));
h5create(filename,['/marker/gc_' id '/time'],size(time));
h5create(filename,['/marker/gc_' id '/id'],size(idnum),'Datatype','int32');


% Write Attributes
h5writeatt(filename,'/marker','active',id,'TextEncoding','system');
h5writeatt(filename,['/marker/gc_' id ],'date',datestr(now,'yyyy-mm-dd hh:MM:ss'),'TextEncoding','system');
h5writeatt(filename,['/marker/gc_' id ],'description','Written by MATLAB','TextEncoding','system');

% Relativistic adjustments
vtotal = sqrt(2.*energy.*ec./(mass.*amu));
gammarel = sqrt(1.0 ./ ( (1.0+vtotal./cspeed) .* (1.0-vtotal./cspeed)));
energycorr = (gammarel - 1.0).*mass.*amu.*cspeed.*cspeed./ec;

% Write Variables
h5write(filename,['/marker/gc_' id '/n'],int32(length(r)));
h5write(filename,['/marker/gc_' id '/r'],r);
h5write(filename,['/marker/gc_' id '/phi'],phi);
h5write(filename,['/marker/gc_' id '/z'],z);
h5write(filename,['/marker/gc_' id '/energy'],energycorr);
h5write(filename,['/marker/gc_' id '/pitch'],pitch);
h5write(filename,['/marker/gc_' id '/zeta'],zeta);
h5write(filename,['/marker/gc_' id '/mass'],mass);
h5write(filename,['/marker/gc_' id '/charge'],int32(charge));
h5write(filename,['/marker/gc_' id '/anum'],int32(anum));
h5write(filename,['/marker/gc_' id '/znum'],int32(znum));
h5write(filename,['/marker/gc_' id '/weight'],weight);
h5write(filename,['/marker/gc_' id '/time'],time);
h5write(filename,['/marker/gc_' id '/id'],int32(idnum));


end

