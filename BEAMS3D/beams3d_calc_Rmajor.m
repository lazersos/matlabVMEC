function Rmajor = beams3d_calc_Rmajor(beam_data)
%BEAMS3D_CALC_RMAJOR Calculates the major radius from a BEAMS3D run
%   The BEAMS3D_CALC_RMAJOR routine calculates the major radius of an
%   equilbirium from the background grid value S_ARR stored in a beams3d
%   output data struture (as produced by READ_BEAMS3D).  
%
%   Example:
%       beam_data=read_beams3d('beams3d_test.h5');
%       Rmajor = beams3d_calc_Rmajor(beam_data);
%
%   Maintained by: Samuel Lazerson (samuel.lazerson@ipp.mpg.de)
%   Version:       1.0


Rmajor = zeros(1,beam_data.ns_prof1);
s = beam_data.rho.*beam_data.rho;
for j = 1:beam_data.nphi
    S2D = squeeze(beam_data.S_ARR(:,j,:));
    S2D(S2D>1.2) = 1.2;
    for k = 1:beam_data.ns_prof1
        cc=contourc(beam_data.raxis,beam_data.zaxis,S2D',[s(k) s(k)]);
        cc = cc(:,2:end-1);
        if isempty(cc) || size(cc,2) < 16
            cc=[];
            [i,~]=find(S2D == min(S2D,[],'all'),1);
            r0 = beam_data.raxis(i);
            cc(1,1) = r0;
            %cc(1,:) = cc(1,:)-r0;
            %cc(2,:) = cc(2,:)-z0;
        end
        Rmajor(k) = Rmajor(k) + mean(cc(1,:));
    end
end
Rmajor = Rmajor./double(beam_data.nphi);



end