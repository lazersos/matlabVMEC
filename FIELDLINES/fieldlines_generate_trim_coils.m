func = @(Icoil, phase) fprintf("Phase %d: %2.8E %2.8E %2.8E %2.8E %2.8E\n",phase, [72 48 48 48 48].*Icoil.*cosd(72.*[1 2 3 4 0]-phase));

Icoil = 200;

phase = linspace(0, 360, 11);

for i=1:length(phase)
   func(Icoil, phase(i));
end