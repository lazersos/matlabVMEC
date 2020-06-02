function ascot5_1dto1ds(a5file,plasmaid)
%ASCOT5_1DSTO1D Copies 1D data to 1DS
%   The ASCOT5_1DTO1DS function copies a 1D plasma ID to a 1DS plasma id
%   in an ASCOT5 file setting 1DS as active.  It takes a ASCOT5 filename and
%   a plasma id number as input.  If an empty array is passed to plasmaid
%   then the currently active ID is used.
%   
%   Example:
%       ascot5_1dto1ds('ascot5.h5',[]);
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00


% Check for file
if ~isfile(a5file)
    disp(['ERROR: ' a5file ' file not found!']);
    return;
end

% Use active
if isempty(plasmaid)
    plasmaid=h5readatt(a5file,'/plasma','active');
    disp(['  Using plasma: ' plasmaid]);
end

%Pull Profile information
path = ['/plasma/plasma_1D_' num2str(plasmaid,'%10.10i')];
try
    rho = h5read(a5file,[path '/rho']);
    nion = h5read(a5file,[path '/nion']);
    nrho = h5read(a5file,[path '/nrho']);
    znum = h5read(a5file,[path '/znum']);
    anum = h5read(a5file,[path '/anum']);
    charge = h5read(a5file,[path '/charge']);
    mass = h5read(a5file,[path '/mass']);
    ne = h5read(a5file,[path '/edensity']);
    te = h5read(a5file,[path '/etemperature']);
    ni = h5read(a5file,[path '/idensity']);
    ti = h5read(a5file,[path '/itemperature']);
catch
    disp(['ERROR: Could not find plasma: ' num2str(plasmaid,'%10.10i')]);
    return;
end

% Write to new file
ascot5_writeplasma(a5file,znum,anum,charge,mass,rho,ne,ni,te,ti,'1DS','fix_profs')

return;


end

