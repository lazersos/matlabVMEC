function [s, deposition] = beams3d_deposition(beam_data)
%BEAMS3D_BIRTHS Calculates neutral beam deposition profile
%   The BEAMS3D_DEPOSITION routine calculates the neutral beam deposition
%   given a BEAMS3D data structure as returned by the routine READ_BEAMS3D.
%   The routine returns the normalized toroidal flux grid and an array for
%   the number of deposited particles in units of part/s.  The array is
%   broken down by beam population.
%
% Example usage
%      beam_data=read_beams3d('beams3d_test.h5');  
%      [s, deposition]=beams3d_deposition(beam_data);
%
% Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
% Version:       1.0

s=[];
deposition=[];

% Check for a neutral beam run
if ~all(beam_data.neut_lines(1,:))
    disp(' No neutral data found!');
    return;
end

%Initialize
s = zeros(1,beam_data.ns_prof);
deposition = zeros(beam_data.nbeams,beam_data.ns_prof);

% Get volume elements
[sv, ~, vp] = beams3d_volume(beam_data);
vp = 2.*sqrt(sv).*vp; % dV/ds => dV/drho

% Downselect to deposited particles
dex = beam_data.neut_lines(3,:)==0;
f   = sqrt(beam_data.S_lines(2,dex)); %rho
w   = beam_data.Weight(dex);
b   = beam_data.Beam(dex);

% Bin with weight
dh = beam_data.rho(1);
for i=1:beam_data.nbeams
    for j=1:beam_data.ns_prof
        dex = and(f>=beam_data.rho(i)-dh,f<beam_data.rho(i)+dh);
        dex = and(b==i,dex);
        deposition(i,j) = sum(w(dex)).*vp(j);
    end
end

% Calc s
s = rho.*rho;


end

