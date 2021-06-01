function U_BEAM = beams3d_fixUlines(beam_data)
%BEAMS3D_FIXULINES Corrects errors in U_lines by recalculating
%   The BEAMS3D_FIXULINES code corrects possible errors in U_lines
%   attributed to the fact that U and a discontinuity in it when
%   represented on an R/Z grid.  
%
%   Example
%       beam_data=read_beams3d('beams3d_test.h5');
%       beam_data.U_lines=beams3d_fixUlines(beam_data);
%
%   Maintained by: Samuel Lazerson (sameul.lazerson@ipp.mpg.de)
%   Version:       1.00

phimax = beam_data.phiaxis(end);
R_BEAM = beam_data.R_lines;
P_BEAM = mod(beam_data.PHI_lines,phimax);
Z_BEAM = beam_data.Z_lines;
X = permute(cos(beam_data.U_ARR),[2 1 3]);
Y = permute(sin(beam_data.U_ARR),[2 1 3]);
X_BEAM= interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
    X,R_BEAM,P_BEAM,Z_BEAM);
Y_BEAM= interp3(beam_data.raxis,beam_data.phiaxis,beam_data.zaxis,...
    Y,R_BEAM,P_BEAM,Z_BEAM);
U_BEAM = atan2(Y_BEAM,X_BEAM);

return;

end

