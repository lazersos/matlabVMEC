function [s, pot, dpotds] = beams3d_er(beam_data)
%BEASM3D_ER Extracts potential and E profile data from BEAMS3D Data
%   The BEAMS3D_ER function calculates the radial potential profiles
%   (in flux) from the background grid data.  The returned profiles are on
%   a 128 point grid.
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [s, pot, dpotds] = beams3d_er(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.5
%        Note earlier versions had incorrectly added a negative sign.


ns = 128;
ns1 = ns - 1;
ds = 1./(ns1);
s = (0:ns1).*ds;
pot=[];
dpotds=[];

% Extract
[C, IA, ~] = unique(beam_data.S_ARR);
pot_temp = beam_data.POT_ARR(IA);

% Make mirror
C=[-flipud(C); C];
pot_temp = [flipud(pot_temp); pot_temp];

% spline
[C, IA, ~] = unique(C);
pot = pchip(C,pot_temp(IA),s);
shalf = 0.5.*(s(1:ns1)+s(2:ns));
dEds = (pot(2:ns)-pot(1:ns1))./(s(2:ns)-s(1:ns1));
dpotds = [0 pchip(shalf,dEds,s(2:end))];

end

