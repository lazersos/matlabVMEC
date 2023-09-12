function data = gmsh2mumat(filename)
%GMSH2MUMAT Converts a GMSH output .m file into a mumat structure
%   This subroutine takes a .m file as produced by the GMSH code and
%   outputs a data structure as returned by READ_MUMAT.  Note that the
%   state function for this material is defaulted to type=3 (soft magnet
%   with constant permeability), with a value of 500.
%
%   Example:
%       sphere_data=read_mumag('sphere.dat');
%
%   Written by:     S.Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:        1.0
%   Date:           9/12/23

data = [];

run(filename);

% Check file
if ~isvarname('msh')
    disp(['ERROR: Was ' filename ' produced by GMSH?']);
    return;
end

if ~isfield(msh,'POS')
    disp('ERROR: Nos POS field found in msh.');
end

if ~isfield(msh,'TETS')
    disp('ERROR: No TETS field found in msh.');
end

data.machine = filename;
data.date = datestr(today);
data.coords = msh.POS';
data.tet = msh.TETS(:,1:4)';
data.nvertex = size(data.coords,2);
data.ntet = size(data.tet,2);
data.func_dex = ones(1,data.ntet);
data.nstate = 1;
data.state_func(1).type = 3;
data.state_func(1).mu = 500;

return;
end