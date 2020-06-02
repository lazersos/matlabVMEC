function beams3d_output_particles(beam_data,tdex,pdex)
%beams3d_output_particles Outputs select particles to screen
%   The BEAMS3D_OUTPUT_PARTICLES subroutine displays a set of particles to
%   the screen from a beams3d run. The user can specify the time indices
%   (tdex) and particle indicies (pdex) to print to screen. The time
%   indicies can be either a single value or the same size as pdex.
%
% Example usage
%      data=read_beams3d('fieldline_test.h5');  % Reads BEAMS3D HDF5 file
%      beams3d_output_particles(data,2,[5 25 150 170]);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

% Check inputs
if (tdex > beam_data.npoinc)
    disp('  ERROR: tdex > NPOINC!');
    disp(['     NPOINC: ' num2str(beam_data.npoinc,'%i')]);
    return;
end
if (any(pdex>beam_data.nparticles))
    disp('  ERROR: pdex > nparticles');
    disp(['      max(pdex): ' num2str(max(pdex),'%i')]);
    disp(['     nparticles: ' num2str(beam_data.nparticles,'%i')]);
    return;
end

% Handle tdex
if length(tdex)>1
    timedex=tdex;
else
    timedex = ones(1,length(pdex)).*tdex;
end
j = length(timedex);
k = length(pdex);
if (j~=k)
    disp('  ERROR: length(tdex) .ne. length(pdex)');
    return;
end

% Create Strings
R_STR=[];
Z_STR=[];
P_STR=[];
MU_STR = [];
VLL_STR = [];
for i = 1:length(pdex)
    j=timedex(i);
    k=pdex(i);
    R_STR = [R_STR '  ' num2str(beam_data.R_lines(j,k),'%20.10E')];
    Z_STR = [Z_STR '  ' num2str(beam_data.Z_lines(j,k),'%20.10E')];
    P_STR = [P_STR '  ' num2str(beam_data.PHI_lines(j,k),'%20.10E')];
    MU_STR = [MU_STR '  ' num2str(beam_data.moment_lines(j,k),'%20.10E')];
    VLL_STR = [VLL_STR '  ' num2str(beam_data.vll_lines(j,k),'%20.10E')];
end

% Print to Screen
j=length(pdex);
k=mean(tdex)./double(beam_data.npoinc-1);
disp(['  R_START_IN = ' R_STR]);
disp(['  Z_START_IN = ' Z_STR]);
disp(['  PHI_START_IN = ' P_STR]);
disp(['  MU_START_IN = ' MU_STR]);
disp(['  VLL_START_IN = ' VLL_STR]);
disp(['  CHARGE_IN = ' num2str(j,'%i*') num2str(mean(beam_data.charge(pdex)),'%20.10E')]);
disp(['  MASS_IN = ' num2str(j,'%i*') num2str(mean(beam_data.mass(pdex)),'%20.10E')]);
disp(['  ZATOM_IN = ' num2str(j,'%i*') num2str(mean(beam_data.Zatom(pdex)),'%20.10E')]);
disp(['  T_END_IN = ' num2str(j,'%i*') num2str(k.*mean(beam_data.t_end(pdex)),'%20.10E')]);
end

