function [s, deposition] = beams3d_deposition(beam_data)
%BEAMS3D_BIRTHS Calculates neutral beam deposition profile
%   The BEAMS3D_DEPOSITION routine calculates the neutral beam deposition
%   given a BEAMS3D data structure as returned by the routine READ_BEAMS3D.
%   The routine returns the normalized toroidal flux grid and an array for
%   the number of deposited particles in units of part*m^-3/s.  The array is
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
bdex = 2;

% Check for a neutral beam run
if ~all(beam_data.neut_lines(1,:))
    disp(' No neutral data found!');
    bdex=1;
%    return;
end

%Initialize
s = zeros(1,beam_data.ns_prof);
deposition = zeros(beam_data.nbeams,beam_data.ns_prof);

% Get volume elements
[sv, ~, vp] = beams3d_volume(beam_data);
drho = 1./beam_data.ns_prof;
vp = pchip(sqrt(sv),2.*sqrt(sv).*vp,beam_data.rho)'.*drho; % dV/ds => dV=dV/drho*drho

% Downselect to deposited particles
if bdex==2
    dex = and(beam_data.neut_lines(3,:)==0,beam_data.S_lines(2,:)<1);
else
    dex = 1:double(beam_data.nparticles);
end
f   = sqrt(beam_data.S_lines(2,dex)); %rho
w   = beam_data.Weight(dex);
b   = beam_data.Beam(dex);

% Bin with weight
for i=1:beam_data.nbeams
    dex = b==i;
    w2 = w(dex)';
    f2 = f(dex)';
    deposition(i,:) = hist1d_weighted(f2,w2,beam_data.rho).*vp;
end

% Calc s
rho = beam_data.rho;
s = rho.*rho;


end

