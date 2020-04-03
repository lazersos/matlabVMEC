function beams3d_zavg(narr,zarr,marr)
%BEAMS3D_ZAVG Calculates the ZAVG and ZMEAN BEAMS3D parameters
%   The BEAMS3D_ZAVG function takes three array parameters the density
%   (narr), charge number (zarr), and mass (marr).  The plasma mass is
%   assumed to be the density weigted average mass.  Note that Z and n are
%   unit agnostic.  Mass can be input in kg or amu.
%
%   Example (D, H, C)
%       narr = [0.9E20, 0.09E20,       0.01E20];
%       zarr = [1.0,    1.0,           6.0];
%       marr = [1.008,  2.01410177811, 12.011];
%       beams3d_zavg(narr,zarr,marr);
%       
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.00

amu = 1.66053906660E-27;
% Assume plasma mass is first mass
marr2 = marr;
if marr>0.1, marr2 = marr2*amu; end
mp = sum(marr2.*narr)./sum(narr);
nzz = sum(narr.*zarr.*zarr);
nz  = sum(narr.*zarr);
nzzm = sum(narr.*zarr.*zarr.*marr2./mp);
Zavg = nzz./nz;
Zmean = nzzm./nz;
disp(['  PLASMA_MASS = ' num2str(mp,'%20.12E')]);
disp(['  ZAVG = ' num2str(Zavg,'%20.12E')]);
disp(['  ZMEAN = ' num2str(Zmean,'%20.12E')]);
return;

end

