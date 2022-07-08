function Aminor = beams3d_calc_aminor(beam_data)
%BEAMS3D_CALC_AMINOR Calculates the minor radius from a BEAMS3D run
%   The BEAMS3D_CALC_AMINOR routine calculates the minor radius of an
%   equilbirium from the background grid value S_ARR stored in a beams3d
%   output data struture (as produced by READ_BEAMS3D).  
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       Aminor = beams3d_calc_aminor(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


Aminor = 0;
for j = 1:beam_data.nphi
    S2D = squeeze(beam_data.S_ARR(:,j,:));
    S2D(S2D>1.2) = 1.2;
    [i,k]=find(S2D == min(S2D,[],'all'),1);
    r0 = beam_data.raxis(i);
    z0 = beam_data.zaxis(k);
    cc=contourc(beam_data.raxis,beam_data.zaxis,S2D',1.0);
    cc = cc(:,2:end-1);
    cc(1,:) = cc(1,:)-r0;
    cc(2,:) = cc(2,:)-z0;
    Aminor = Aminor + mean(sqrt(sum(cc.^2,1)));
end
Aminor = Aminor./double(beam_data.nphi);



end