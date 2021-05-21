function [] = fieldlines_generate_trim_coils(Icoil, phase)
% fieldlines_generate_trim_coils 
% Small function to print out trim coil current that need to be set for a
% certain phase and current sweep
%
% Inputs: Icoil: the current desired
% Inputs: phase: array of phases desired
% Output: printed in terminal settings for trim coils, ready to be copied
%
%   Created by: D.J. Engels (d.j.engels@student.tue.nl)
%   Version:    1
%   Date:       May 2021

func = @(Icoil, phase) fprintf("Phase %d: %2.8E %2.8E %2.8E %2.8E %2.8E\n",phase, [72 48 48 48 48].*Icoil.*cosd(72.*[1 2 3 4 0]-phase));

for i=1:length(phase)
   func(Icoil, phase(i));
end
end