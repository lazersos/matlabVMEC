function [B0] = beams3d_b0(beam_data)
%BEASM3D_B0 Extracts the field on axis
%   The BEASM3D_B0 function returns the magnetic field on axis.
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       b0 = beams3d_b0(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


B0=[];
% Assume we want zeta=0 plane
phidex=1;

% Make 2D
S2D = squeeze(beam_data.S_ARR(:,phidex,:));
B_PHI2D = squeeze(beam_data.B_PHI(:,phidex,:));
B_R2D = squeeze(beam_data.B_R(:,phidex,:));
B_Z2D = squeeze(beam_data.B_Z(:,phidex,:));
B2D = sqrt(B_R2D.*B_R2D+B_PHI2D.*B_PHI2D+B_Z2D.*B_Z2D);

% Find minimum
smin = min(min(S2D));
dex = S2D == smin;

% Return Value
B0 = mean(B2D(dex));

return;

end

