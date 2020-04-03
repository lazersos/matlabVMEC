function [plasma_mass, plasma_zavg, plasma_zmean]=beams3d_zavg(narr,zarr,marr,varargin)
%BEAMS3D_ZAVG Calculates the ZAVG and ZMEAN BEAMS3D parameters
%   The BEAMS3D_ZAVG function takes three array parameters the density
%   (narr), charge number (zarr), and mass (marr).  The plasma mass is
%   assumed to be the density weigted average mass.  Note that Z and n are
%   unit agnostic.  Mass can be input in kg or amu.  The function accepts
%   'quiet' as an optional argument to supress screen output.
%
%   Example (D, H, C)
%       narr = [0.9E20, 0.09E20,       0.01E20];
%       zarr = [1.0,    1.0,           6.0];
%       marr = [1.008,  2.01410177811, 12.011];
%       [plasma_mass, plasma_zavg, plasma_zmean] = beams3d_zavg(narr,zarr,marr);
%       
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.00

% Defaults
amu = 1.66053906660E-27;
plasma_mass=[];
plasma_zavg=[];
plasma_zmean=[];
lverb=1;

% Handle varargin
if nargin > 3
    i=1;
    while i <= length(varargin)
        switch varargin{i}
            case{'quiet','QUIET'}
                lverb = 0;
            otherwise
                disp(['Unrecognized Option: ' varargin{i}]);
                return
        end
        i = i + 1;
    end
end

% Assume plasma mass is first mass
marr2 = marr;
if marr>0.1, marr2 = marr2*amu; end
plasma_mass = sum(marr2.*marr2.*narr)./sum(marr2.*narr);
nzz = sum(narr.*zarr.*zarr);
nz  = sum(narr.*zarr);
nzzm = sum(narr.*zarr.*zarr.*marr2./plasma_mass);
plasma_zavg = nzz./nz;
plasma_zmean = nzzm./nz;
if lverb
    disp(['  PLASMA_MASS = ' num2str(plasma_mass,'%20.12E')]);
    disp(['  PLASMA_ZAVG = ' num2str(plasma_zavg,'%20.12E')]);
    disp(['  PLASMA_ZMEAN = ' num2str(plasma_zmean,'%20.12E')]);
end
return;

end

