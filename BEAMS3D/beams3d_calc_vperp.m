function vperp = beams3d_calc_vperp(beam_data)
%BEAMS3D_CALC_VPERP Returns perpendicular velocity
%   The BEAMS3D_CALC_VPERP routine returns the perpendicular velocity based
%   on the magnetic moment from a run of the BEAMS3D code.
%
% Example usage
%      data=read_beams3d('fieldline_test.h5');  % Reads BEAMS3D file
%      vperp = beams3d_calc_vperp(data);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

vperp = [];

nsteps = beam_data.npoinc+1;
mass = repmat(beam_data.mass',[nsteps 1]);
mu = beam_data.moment_lines;
B  = beam_data.B_lines;
vperp = sqrt(2.*mu.*B./mass);
vperp(B<=0) = 0;
return
end

