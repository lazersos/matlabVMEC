function [time,loss] = beams3d_calc_loss(beam_data)
%BEASM3D_CALC_LOSS Calculates losses as a function of time
%   The BEAMS3D_CALC_LOSS routine calculates the cumulative losses of
%   particles (not markers) for each beam line. It returns a vector of time
%   and a loss matrix of size (nbeam,ntime).
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       [time,loss] = beams3d_calc_loss(beam_data);
%       plot(time,loss'); set(gca,'XScale','log');
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0

% Calculate time array
time = [linspace(0,0.999E-6,1001) linspace(1E-6,0.999E-3,1001) linspace(1E-3,1,10001)];

% Initialize loss array
loss = zeros(beam_data.nbeams,length(time));

for i = 1:beam_data.nbeams
    lostdex = and(beam_data.end_state==2,beam_data.Beam==i);
    t = beam_data.t_end(lostdex);
    w = beam_data.Weight(lostdex);
    counts = hist1d_weighted(t,w,time);
    nlost = cumsum(counts);
    loss(i,:) = nlost;
end

end