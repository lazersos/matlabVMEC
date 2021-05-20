func = @(Icoil, phase) fprintf("Phase %d: %2.8E %2.8E %2.8E %2.8E %2.8E\n",phase, [72 48 48 48 48].*Icoil.*cosd(72.*[1 2 3 4 0]-phase));

Icoil = 116.75;

phase = linspace(150, 150, 1);

for i=1:length(phase)
   func(Icoil, phase(i));
end