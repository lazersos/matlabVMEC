function [] = generate_coils_for_sweep(type, Icoil, phase)
% fieldlines_generate_coils 
% Small function to print out trim or sweep coil current that need to be set for a
% certain phase and current sweep
%
% Inputs: type: Type of coil (trim/sweep)
% Inputs: Icoil: the current desired
% Inputs: phase: array of phases desired
% Output: printed in terminal settings for coils, ready to be copied
%
%   Created by: D.J. Engels (d.j.engels@student.tue.nl)
%   Version:    1.1
%   Date:       May 2021


switch type
    case 'trim'
        func = @(Icoil, phase) fprintf("Phase %d: %2.8E %2.8E %2.8E %2.8E %2.8E\n",phase, [72 48 48 48 48].*Icoil.*cosd(72.*[1 2 3 4 0]-phase));
    case 'sweep'
        phid_icc=[14.64049     -14.64049     -57.35951     -86.64049     -129.3595     -158.6405      158.6405      129.3595      86.64049      57.35951];
        func = @(Icoil, phase) fprintf("Phase %d: %2.8E %2.8E %2.8E %2.8E %2.8E %2.8E %2.8E %2.8E %2.8E %2.8E\n", phase, Icoil.*8.*cosd(2.*phid_icc+phase));
        
    otherwise
        fprintf("Unknown type of coils given. Use trim or sweep\n");
        return
end

for i=1:length(phase)
   func(Icoil, phase(i));
end
end